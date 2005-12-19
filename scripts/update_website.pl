#!/usr/local/bin/perl5.8.0 -w
#
# update_website.pl
# 
# by Keith Bradnam aged 12 and half (is this the reincarnation of Peter Pan?)
#
# A script to finish the last part of the weekly build by updating all of the
# relevant WormBase and Wormpep web pages.
#
# Last updated by: $Author: mh6 $     
# Last updated on: $Date: 2005-12-19 16:01:01 $      


#################################################################################
# load modules etc.                                                             #
#################################################################################
use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Symbol 'gensym';
use Ace;
use Carp;

##############################
# command-line options       #
##############################

my $debug;              # Verbose debug mode careful: only chanegs source dirs
my $help;               # Help/Usage page
my $all;                # option all - do everything
my $header;             # create release_header.shtml
my $dna;                # create DNA.shtml
my $finished;           # create finished.shtml
my $dbcomp;             # create dbcomp.shtml
my $wormpep_diff;       # diff log for wormpep releases
my $copyGFF;            # copy GFF files to WWW site 
my $test;               # test run (requires a release number) careful: it doesnt change anything
my $directories;        # create website directories
my $overlap;            # run copy_overlapcheck_files
my $EST_files;          # run copy_EST_files
my $create_GFF;         # run create_GFF_intron_files
my $update_wormpep;
my $store;		# to specifiy a file containing stored arguments

# check for command-line options if none given then you do everything
unless (defined $ARGV[0]) {
    $all = 1;
}

GetOptions ("all"            => \$all,
	    "header"         => \$header,
	    "dna"            => \$dna,
	    "finished"       => \$finished,
            "dbcomp"         => \$dbcomp,
            "wormpepdiff"    => \$wormpep_diff,
	    "copygff"        => \$copyGFF,
            "debug"          => \$debug,
	    "test=s"         => \$test,
            "help"           => \$help,
            "h"              => \$help,
	    "dirs"           => \$directories,
	    "overlap"        => \$overlap,
	    "est"            => \$EST_files,
	    "create_gff"     => \$create_GFF,
	    "wormpep"        => \$update_wormpep,
	    'store=s'	=> \$store
	   );


############################
# recreate configuration   #
############################
my $wb;
if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug, -test => $test, ) }

###########################################
# Variables Part II (depending on $wb)    #
###########################################
$test  = $wb->test  if $wb->test;     # Test mode
$debug = $wb->debug if $wb->debug;    # Debug mode, output only goes to one user



##############################
# Script variables           #
##############################
my $day;
my $month;
my $year;
($day, $month, $year) = (localtime)[3,4,5];
my $maintainers      = "All";
my $rundate          = `date +%y%m%d`; chomp $rundate;
my $runtime          = `date +%H:%M:%S`; chomp $runtime;
my @chrom            = qw ( I II III IV V X ); 
my $WS_current       = $wb->get_wormbase_version;
my $release_date2    = `date +%d/%m/%y`; chomp $release_date2;
my $release_date     = sprintf ("%02d %02d %04d", $day, $month+1, $year+1900);
my $WS_previous      = $WS_current - 1;
my $WS_name          = $wb->get_wormbase_version_name;
my $WS_previous_name = "WS".$WS_previous;


# file path info
my $www              = "/nfs/WWWdev/SANGER_docs/htdocs/Projects/C_elegans/WORMBASE";
my $wwwdata          = "/nfs/WWWdev/SANGER_docs/data/Projects/C_elegans";
my $gff              = $wb->gff_splits;
my $dbpath           = $wb->autoace;
my $basedir	     = $wb->basedir;
my $chromdir	     = $wb->chromosomes;

#make new log
my $log = Log_files->make_build_log($wb);

#############
# TEST MODE
############

# Use debug mode?
if ($debug) {
    print "DEBUG = \"$debug\"\n\n";
     $maintainers = $debug . '\@sanger.ac.uk';
}


###########################################################################################################
# Main subroutine calls    
###########################################################################################################

&create_log_file;

&create_web_directories        if ($all || $directories);
&create_release_header         if ($all || $header);    # makes release_header.shtml
&create_DNA_table              if ($all || $dna);       # makes DNA.shtml
&create_finished_seq           if ($all || $finished);  # makes finished_seqs.shtml
&create_dbcomp                 if ($all || $dbcomp);    # makes dbcomp.shtml
&create_wormpep_page           if ($all || $wormpep_diff); # the wormpep page accessible from the main wormbase page
&copy_overlapcheck_files       if ($all || $overlap); # copies check files/wormsrv2/autoace/CHECKS and converts to HTML etc.
&copy_EST_files                if ($all || $EST_files); # copies EST_analysis.html and then edits this file then copies other EST files
&copy_GFF_files                if ($all || $copyGFF); # copies some of the GFF files in /wormsrv2/autoace/GFF_SPLITS
&create_GFF_intron_files       if ($all || $create_GFF); # copies the GFF_introns_confirmed_CDS* files from previous WS release, updates totals
&update_wormpep_pages          if ($all || $update_wormpep); # update the main wormpep web pages



#########################
# final bit of tidying up
#########################

# update 'current' symlink on dev site
$log->write_to("\nChanging 'current symbolic link to point to new release\n"); 	
chdir("$www") || $log->write_to("Couldn't chdir to $www\n");
system("rm -f $www/current") && croak "Couldn't remove 'current' symlink\n"; 	 
system("ln -s $WS_name current") && croak "Couldn't create new symlink\n";


$log->write_to("\n\nC'est finis\n\naus is\n\n");

# mail log
$log->mail( "$maintainers", "BUILD REPORT: $0" );


exit(0);


###########################################################################################################
# Create log file          
###########################################################################################################
sub create_log_file{
  $log->write_to("# update_website.pl\n\n");     
  $log->write_to("# run details    : $rundate $runtime\n");
  $log->write_to("# WormBase version : $WS_name\n\n\n");
}

###########################################################################################################
# Create new directories on website     
###########################################################################################################
sub create_web_directories{

  if (! -e "$www/$WS_name"){
    system("mkdir $www/$WS_name") && warn "Cannot create $WS_name directory $!\n";
    system("mkdir $www/$WS_name/Checks") && warn "Cannot create $WS_name/Checks directory $!\n";
    system("mkdir $www/$WS_name/GFF") && warn "Cannot create $WS_name/GFF directory $!\n";
  }

}



###########################################################################################################
# Copy overlapcheck files to website    
###########################################################################################################

sub copy_overlapcheck_files{

  $log->write_to("\ncopy_overlapcheck_files\n");
  $log->write_to("-----------------------\n");

  # list of files to be copied
  my @filenames = qw( overlapping_TSL_cam overlapping_TSL_stl overlapping_genes_cam overlapping_genes_stl EST_in_intron_cam EST_in_intron_stl repeat_in_exon_cam repeat_in_exon_stl );

  $log->write_to("copying from $www/$WS_previous_name/Checks/index.shtml to $www/$WS_name/Checks\n"); 
  system("cp -f $www/$WS_previous_name/Checks/index.shtml $www/$WS_name/Checks/") && warn "Cannot copy index.shtml $!\n";

  $log->write_to("copying three ace files from $www/$WS_previous_name/Checks/*.ace to $www/$WS_name/Checks\n"); 
  system("cp -f $www/$WS_previous_name/Checks/*.ace $www/$WS_name/Checks/") && warn "Cannot copy *.ace files $!\n";

  $log->write_to("copying short genes file $dbpath/CHECKS/short_spurious_genes.$WS_name.csv.html to $www/$WS_name/Checks/short_genes.html\n"); 
  system("cp -f $dbpath/CHECKS/short_spurious_genes.$WS_name.csv.html $www/$WS_name/Checks/short_genes.html") && warn "Cannot copy short_genes file $!\n";

  $log->write_to("Creating symbolic link for header.ini file\n");
  chdir("$www/$WS_name/Checks") || $log->write_to("Couldn't chdir to $www/$WS_name/Checks\n");
  system("ln -sf ../../header.ini header.ini") && croak "Couldn't create new symlink\n";

  foreach my $file (@filenames) {
    my $line_total;
    my @line_stats;


    foreach my $chrom (@chrom) {
      $log->write_to("Copying file CHROMOSOME_$chrom.$file to $www/$WS_name/Checks\n");
      next unless (-s "$dbpath/CHECKS/CHROMOSOME_$chrom.$file");
      system ("cp -f $dbpath/CHECKS/CHROMOSOME_$chrom.$file $www/$WS_name/Checks/")
	&& croak "Could not copy $dbpath/CHECKS/CHROMOSOME_$chrom.$file $!\n";
      
      $log->write_to("Calculating line numbers for CHROMOSOME_$chrom.$file\n");
      my $line = `wc -l $dbpath/CHECKS/CHROMOSOME_$chrom.$file`;
      # take line count and add to array
      my ($new) = ($line =~ /(\d+)/);
      push @line_stats, $new;
      # add to line total for all chromosomes
      $line_total += $new;
    }
    
    
    # make new html files
    $log->write_to("Generating new html files in $www/$WS_name/Checks/\n");
    my $fh     = gensym();
    my $newfh  = gensym();
    my $count = 0;
    
    open ($fh, "$www/$WS_previous_name/Checks/$file.html") || croak "Cannot open old html file $!\n";
    open ($newfh, ">$www/$WS_name/Checks/$file.html") || croak "Cannot open new html file $!\n";
    
    $log->write_to("Generating $www/$WS_name/Checks/$file.html\n");
    
    while (<$fh>) {
      if ((/<TD ALIGN=\"center\" COLSPAN=\"2\"><B>\s*(\d+)\s+\[\s*\d+\]<\/B><\/TD>/) && ($count < 6)) {
	my $old = $1;
	print $newfh "<TD ALIGN=\"center\" COLSPAN=\"2\"><B> ".$line_stats[$count]." [$old]</B></TD>\n";
	$count++;
      }
      elsif ((/<TD ALIGN=\"center\" COLSPAN=\"2\"><B>\s*(\d+)\s+\[\s*\d+\]<\/B><\/TD>/) && ($count > 5)) {
	my $old = $1;
	$count++;
	print $newfh "<TD ALIGN=\"center\" COLSPAN=\"2\"><B> ".$line_total." [$old]</B></TD>\n";
      }
      else {
	print $newfh $_;
      }
    }    
    close($fh)    || croak "Cannot close old html file $!\n";
    close($newfh) || croak "Cannot close new html file $!\n";
  }
}


###########################################################################################################
# Creates: release_header.shtml, finished_seqs.shtml, DNA_table.shtml,
# dbcomp.shtml, and wormpep.shtml
###########################################################################################################

sub create_release_header {

  # create release_header.shtml  
  $log->write_to("\ncreate_top_level_web_pages\n");
  $log->write_to("--------------------------\n");
  
  $log->write_to("Creating $www/$WS_name/release_header.shtml\n");
  open  (RELEASE_HEADER, ">$www/$WS_name/release_header.shtml") || croak "Couldn't create release_header.shtml\n";
  print  RELEASE_HEADER  "<b>$WS_name release, updated $release_date.</b>\n";
  close (RELEASE_HEADER);

}

sub create_DNA_table {

  # create DNA.table.shtml
  $log->write_to("Creating $www/$WS_name/DNA.table.shtml\n");

  my $DNA_tot;   # total DNA length

  my $DNA_d  = 0;
  my $DNA_dp = 0;

  my $DNA_n  = 0;
  my $DNA_np = 0;

  my $DNA_g  = 0;
  my $DNA_gp = 0;

  my $DNA_a  = 0;
  my $DNA_ap = 0;

  my $DNA_t  = 0;
  my $DNA_tp = 0;

  my $DNA_c  = 0;
  my $DNA_cp = 0;

  open (DNA_table, ">$www/$WS_name/DNA_table.shtml") || croak "Couldn't create '$www/$WS_name/DNA_table.shtml'\n\n";
  print DNA_table "<P>\n";
  
  open (COMPOSITION, "<$chromdir/composition.all") || croak "Failed to open composition.all\n";
  while (<COMPOSITION>) {
    $DNA_tot = $1 if (/6 sequences >= 0, (\d+) total/);
    ($DNA_d,$DNA_dp) = ($1,$2) if (/  - (\d+)\s+(\d+) %/);
    ($DNA_n,$DNA_np) = ($1,$2) if (/  n (\d+)\s+(\d+) %/);
    ($DNA_g,$DNA_gp) = ($1,$2) if (/  g (\d+)\s+(\d+) %/);
    ($DNA_a,$DNA_ap) = ($1,$2) if (/  a (\d+)\s+(\d+) %/);
    ($DNA_t,$DNA_tp) = ($1,$2) if (/  t (\d+)\s+(\d+) %/);
    ($DNA_c,$DNA_cp) = ($1,$2) if (/  c (\d+)\s+(\d+) %/);
  }
  close (COMPOSITION);
    

  print DNA_table "<p>\n";
  print DNA_table "<!--#include virtual=\"/SSI/tabletop.shtml\"-->\n";
  print DNA_table "<table WIDTH=\"350\" CELLSPACING=\"0\" CELLWIDTH=\"0\" BORDER=\"0\">\n";

  print DNA_table "<tr class=\"h2bg\">\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <th align=\"center\" class=\"barialw\">&nbsp;</th>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <th align=\"center\" class=\"barialw\">No. bases</th>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <th align=\"center\" class=\"barialw\">%</th>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "</tr>\n";

  print DNA_table "<tr class=\"violet1\">\n";
  print DNA_table "  <td colspan=\"7\"><img src=\"/icons/blank.gif\" width=\"10\" height=\"10\"></td>\n";
  print DNA_table "</tr>\n";

  print DNA_table "<tr class=\"violet3\">\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"15\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"center\">Total</td>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"15\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"right\" >$DNA_tot</td>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"15\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"right\" >100</td>\n";   
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"15\" height=\"22\"></td>\n";
  print DNA_table "</tr>\n";

  print DNA_table "<tr class=\"violet2\">\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"center\">-</td>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"right\" >$DNA_d</td>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"right\" >$DNA_dp</td>\n";   
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "</tr>\n";

  print DNA_table "<tr class=\"violet3\">\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"center\">G</td>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"right\" >$DNA_g</td>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"right\" >$DNA_gp</td>\n";   
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "</tr>\n";

  print DNA_table "<tr class=\"violet2\">\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"center\">A</td>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"right\" >$DNA_a</td>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"right\" >$DNA_ap</td>\n";   
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "</tr>\n";

  print DNA_table "<tr class=\"violet3\">\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"center\">T</td>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"right\" >$DNA_t</td>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"right\" >$DNA_tp</td>\n";   
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "</tr>\n";

  print DNA_table "<tr class=\"violet2\">\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"center\">C</td>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"right\" >$DNA_c</td>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"right\" >$DNA_cp</td>\n";   
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "</tr>\n";

  print DNA_table "<tr class=\"violet3\">\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"center\">N</td>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"right\" >$DNA_n</td>\n";
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "  <td align=\"right\" >$DNA_np</td>\n";   
  print DNA_table "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DNA_table "</tr>\n";

  print DNA_table "</table>\n";
  print DNA_table "<!--#include virtual=\"/SSI/tablebottom.shtml\"-->\n";

  print DNA_table "</p>\n";
  close (DNA_table);
}


sub create_finished_seq {

  # create finished_seqs.shtml
  $log->write_to("Creating $www/$WS_name/finished_seqs.shtml\n");
  my $db = Ace->connect(-path=>$dbpath) || do { print "Connection failure: ",Ace->error; croak();};
  my $count = $db->fetch(-query=> 'find genome_sequence where Finished');

  open (FINISHED_SEQS, ">$www/$WS_name/finished_seqs.shtml") || croak "Couldn't create finished_seqs.shtml\n";
  print FINISHED_SEQS  "<P>\n";
  print FINISHED_SEQS "$count projects\n";
  print FINISHED_SEQS  "</P>\n";
  close (FINISHED_SEQS);
  $db->close;

}

sub create_dbcomp {

  # create dbcomp.shtml
  $log->write_to("Creating $www/$WS_name/dbcomp.shtml\n");

  system ("cp $dbpath/COMPARE/current.out $wwwdata/WS.dbcomp_output") && croak "Couldn't copy current.out\n";
  my $dbcomplen;
  
  open (DBCOMP, "wc -l $dbpath/COMPARE/current.dbcomp |");
  while (<DBCOMP>) {
      ($dbcomplen) = (/^\s+(\d+)/);
  }
  close DBCOMP;
  
  $dbcomplen = $dbcomplen -2;     # six lines of rubbish in the file but need 4 in the html


  open (DB_comp, ">$www/$WS_name/dbcomp.shtml") || croak "Couldn't create dbcomp.shtml\n\n";
  print DB_comp "<P>\n";
  print DB_comp "<SPACER TYPE=\"horizontal\" SIZE=\"50\">\n";
  print DB_comp "<TABLE WIDTH=\"40%\" CELLSPACING=\"0\" CELLWIDTH=\"0\" BORDER=\"0\">\n";

  print DB_comp "\n<tr VALIGN=\"top\">\n";
  print DB_comp "  <td rowspan=\"$dbcomplen\" class=\"grey1\"> <img src=\"/icons/blank.gif\" width=\"1\"   height=\"1\"></td>\n";
  print DB_comp "  <td colspan=\"9\" class=\"grey1\"> <img src=\"/icons/blank.gif\" width=\"250\" height=\"1\"></td>\n";
  print DB_comp "  <td rowspan=\"$dbcomplen\" class=\"grey1\"> <img src=\"/icons/blank.gif\" width=\"1\"   height=\"1\"></td>\n";
  print DB_comp "</tr>\n\n";

  print DB_comp "<tr class=\"violet3\">\n";
  print DB_comp "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DB_comp "  <th>Class</th>\n";
  print DB_comp "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DB_comp "  <th>$WS_previous_name</th>\n";
  print DB_comp "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DB_comp "  <th>$WS_name</th>\n";
  print DB_comp "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DB_comp "  <th>Diff</th>\n";
  print DB_comp "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print DB_comp "</tr>\n";

  print DB_comp "<tr valign=\"top\">\n";
  print DB_comp "  <td colspan=\"11\" class=\"grey1\"><img src=\"/icons/blank.gif\" width=\"1\" height=\"1\"></td>\n";
  print DB_comp "</tr>\n";

  
  my $dbcompcount = 0;

  open (COMP, "<$dbpath/COMPARE/current.dbcomp") || croak "Failed to open current.dbcomp\n\n";
  while (<COMP>) {
      next if (/\+/);
      next if (/ACEDB/);
      next if (/change/);
      
      my ($nowt,$class,$count_1,$count_2,$diff) = split (/\|/,$_);
      $class =~ s/^\s+//g;
      
      if ($dbcompcount % 2 == 0) {
	  print DB_comp "<tr class=\"blue\">\n";
      }
      else {
	  print DB_comp "<tr class=\"violet2\">\n";
      }
      
      print DB_comp "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
      print DB_comp "  <td align=\"RIGHT\"><A  href=\"/cgi-bin/Projects/C_elegans/wormbase_dbcomp.pl?class=$class\" TARGET=\"_blank\">$class</A></td>\n";
      print DB_comp "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
      print DB_comp "  <td ALIGN=\"RIGHT\">$count_1</td>\n";
      print DB_comp "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
      print DB_comp "  <td ALIGN=\"RIGHT\">$count_2</td>\n";
      print DB_comp "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
      print DB_comp "  <td ALIGN=\"RIGHT\">$diff</td>\n";
      print DB_comp "</tr>\n";
      print DB_comp "\n\n";

      $dbcompcount++;

  }
  close (COMP);

  print DB_comp "<tr valign=\"top\">\n";
  print DB_comp "  <td colspan=\"11\" class=\"grey1\"><img src=\"/icons/blank.gif\" width=\"1\" height=\"1\"></td>\n";
  print DB_comp "</tr>\n";

  print DB_comp "</TABLE>\n";
  print DB_comp "</P>\n";
  close (DB_comp);
}



###########################################################################################################
# Creates the wormpep web page that is part of the WORMBASE web pages
###########################################################################################################

sub create_wormpep_page{

  # create wormpep.shtml  
  $log->write_to("\ncreate_wormpep_diff.shtml\n");
  $log->write_to("-------------------\n");


  $log->write_to("Creating $www/$WS_name/wormpep_diff.shtml\n");

  open (WORMPEP, ">$www/$WS_name/wormpep_diff.shtml") || croak "Failed to create wormpep.shtml\n\n";
  print WORMPEP "<P>\n";
  print WORMPEP "<SPACER TYPE=\"horizontal\" SIZE=\"50\">\n";

  print WORMPEP "<TABLE WIDTH=\"40%\" CELLSPACING=\"0\" CELLWIDTH=\"0\" BORDER=\"0\">\n";

  print WORMPEP "<tr class=\"h2bg\">\n";
  print WORMPEP "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print WORMPEP "  <th align=\"center\" class=\"barialw\">Release</th>\n";
  print WORMPEP "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print WORMPEP "  <th align=\"center\" class=\"barialw\">Sequences</th>\n";
  print WORMPEP "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print WORMPEP "  <th align=\"center\" class=\"barialw\">Letters</th>\n";
  print WORMPEP "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print WORMPEP "  <th align=\"center\" class=\"barialw\">Isoforms</th>\n";
  print WORMPEP "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print WORMPEP "</tr>\n";

  print WORMPEP "<tr class=\"violet1\">\n";
  print WORMPEP "  <td colspan=\"9\"><img src=\"/icons/blank.gif\" width=\"10\" height=\"10\"></td>\n";
  print WORMPEP "</tr>\n";
  
  $log->write_to("Opening log file '$basedir/WORMPEP/wormpep$WS_current/wormpep_current.log'\n");

  my ($wp_seq,$wp_let);
 open (WP_1, "</$basedir/WORMPEP/wormpep$WS_current/wormpep_current.log") || croak "Failed to open wormpep.log\n";
  while (<WP_1>) {
    # No. of sequences (letters) written:  22,221  (9,696,145)
    if (/No\. of sequences \(letters\) written:  (\d+,\d+)  \((\d+,\d+,\d+)\)/) {
      ($wp_seq,$wp_let) = ($1,$2);
    }
  }
  close (WP_1);

  
  # get details from last wormpep log file <= argh
  my $wp_alt;
  my @possible_logs = `ls -t $dbpath/logs/make_wormpep.final.WS${WS_current}*`; # added new logfile name 
  my $wormpeplog  = "$possible_logs[0]";
  open (WP_2, "<$wormpeplog") || croak "Failed to open wormpep.log\n";
  while (<WP_2>) {
    if (/\((\d+)\) alternate splice_forms/) {
      ($wp_alt) = ($1);
    }
  }
  close (WP_2);

  print WORMPEP "<tr class=\"violet3\">\n";
  print WORMPEP "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print WORMPEP "  <td align=\"center\">wormpep${WS_current}</td>\n";
  print WORMPEP "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print WORMPEP "  <td align=\"center\">$wp_seq</td>\n";
  print WORMPEP "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print WORMPEP "  <td align=\"center\" >$wp_let</td>\n";
  print WORMPEP "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print WORMPEP "  <td align=\"center\" >$wp_alt</td>\n";
  print WORMPEP "  <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";

  print WORMPEP "</table>\n";
  print WORMPEP "</P><PRE>\n";
  

  my (@changed, @lost, @new, @reappeared);
  open (WP_3, "<$basedir/WORMPEP/wormpep${WS_current}/wormpep.diff${WS_current}") || croak "Failed to open wormpep.history\n";

  while (<WP_3>) {
    (push (@changed,$_)) if (/changed:/);
    (push (@lost,$_)) if (/lost:/);
    (push (@new,$_)) if (/new:/);
    (push (@reappeared,$_)) if (/reappeared:/);
  }
  close (WP_3);
  
  my $wp_changed = scalar @changed;
  my $wp_lost = scalar @lost;
  my $wp_new = scalar @new;
  my $wp_reappeared = scalar @reappeared;    
  
  print WORMPEP "<FONT SIZE=\"+1\"><B>Modified entries: $wp_changed</B></FONT>\n";
  print WORMPEP "<FONT SIZE=\"-1\">\n";
  foreach (@changed) {
    next if ($_ eq "");
    print WORMPEP $_;
  }
  print WORMPEP "</FONT>\n";  

  print WORMPEP "<FONT SIZE=\"+1\"><B>Deleted entries: $wp_lost</B></FONT>\n";
  print WORMPEP "<FONT SIZE=\"-1\">\n";
  foreach (@lost) {
    next if ($_ eq "");
    print WORMPEP $_;
  }
  print WORMPEP "</FONT>\n";

  print WORMPEP "<FONT SIZE=\"+1\"><B>New entries: $wp_new</B></FONT>\n";
  print WORMPEP "<FONT SIZE=\"-1\">\n";
  foreach (@new) {
    next if ($_ eq "");
    print WORMPEP $_;
  }
  print WORMPEP "</FONT>\n";
  
  print WORMPEP "<FONT SIZE=\"+1\"><B>Reappeared entries: $wp_reappeared</B></FONT>\n";
  print WORMPEP "<FONT SIZE=\"-1\">\n";
  foreach (@reappeared) {
    next if ($_ eq "");
    print WORMPEP $_;
  }
  print WORMPEP "</FONT>\n";
  print WORMPEP "</PRE>\n";
  print WORMPEP "</P>\n";

  close (WORMPEP);
}


###########################################################################################################
# Makes the pages of loci designations for each letter of the alphabet
###########################################################################################################

sub copy_EST_files {

  # make EST_*.shtml using EST_*.txt in /wormsrv2/autoace/CHECKS  
  $log->write_to("\ncopy_EST_files\n");
  $log->write_to("--------------\n");


  my %files = ("EST_no_accession" => "No Accession number", 
               "EST_unassigned" => "Unassigned reads",
               "EST_mismatched" => "mismatched CDS assignments");

  my @line_counts;
  
  foreach my $file (keys(%files)){
    $log->write_to("creating $file.shtml in $www/$WS_name/Checks\n"); 
    open(EST_HTML,">$www/$WS_name/Checks/$file.shtml") || croak "Couldn't create $file.shtml\n";
    
#    print EST_HTML &SangerWeb::virtual_header();
print EST_HTML qq(<!--#include virtual="/perl/header" -->\n);
   
    print EST_HTML "<TABLE WIDTH=\"100%\" CELLSPACING=\"0\" CELLPADDING=\"0\"><TR VALIGN=\"top\" BGCOLOR=\"darkblue\" ><TD WIDTH=\"100%\"><BR><H2 align=\"center\"><FONT COLOR=\"white\">";
    print EST_HTML "EST analysis: \"$files{$file}\"";
    print EST_HTML "</FONT></H2></TD></TR></TABLE>\n";

    print EST_HTML "<P><PRE>";
    open (EST_TXT, "<$dbpath/CHECKS/$file.txt") || croak "Couldn't open EST_total.txt\n";

    my $line = `wc -l $dbpath/CHECKS/$file.txt`;
    # take line count and add to array
    my ($new) = ($line =~ /(\d+)/);
    push(@line_counts,$new);

    while (<EST_TXT>) {
     print EST_HTML "$_";
    }
    print EST_HTML "</PRE></P>";
#    print EST_HTML &SangerWeb::virtual_footer();
    print EST_HTML qq(<!--#include virtual="/perl/footer" -->\n);
    close(EST_TXT);
    close(EST_HTML);
  }

  # need to add 'n/a' into the line_counts array to avoid problems later on
  splice(@line_counts,2,0,"n/a");


  $log->write_to("updating EST_analysis.shtml from $www/$WS_previous/Checks to $www/$WS_name/Checks, adding new info\n"); 

  open(EST_OLD, "<$www/$WS_previous_name/Checks/EST_analysis.html") || croak "Couldn't read old version of EST_analysis.html\n";
  open(EST_NEW, ">$www/$WS_name/Checks/EST_analysis.html") || croak "Couldn't create new version of EST_analysis.html\n";

  my $count=0;
  while (<EST_OLD>) {
      if (/<TD ALIGN=\"center\" COLSPAN=\"2\"><B>.*<\/B><\/TD>/) {
	print EST_NEW "<TD ALIGN=\"center\" COLSPAN=\"2\"><B> ".$line_counts[$count]."</B></TD>\n";
	$count++;
      }
      else {
	print EST_NEW $_;
      }
    }    

  close(EST_NEW);
  close(EST_OLD);
}


###########################################################################################################
# Copy some of the GFF files from $gff and some from CHECKS/
###########################################################################################################

sub copy_GFF_files{
  
  $log->write_to("\ncopy_GFF_files\n");
  $log->write_to("--------------\n");

  $log->write_to("Copying across GFF files from $gff\n");

  #simple double foreach loop to loop through each chromosome and file name
  my @gff_files = ("clone_ends", "clone_path", "exon", "clone_acc", "CDS", "repeats", "intron", "rna");
  foreach my $chrom (@chrom) {
    foreach my $file (@gff_files){
      system("sort -u $gff/CHROMOSOME_$chrom.$file.gff | gff_sort > $www/$WS_name/GFF/CHROMOSOME_$chrom.$file.gff")
	&& croak "Couldn't copy CHROMOSOME_$chrom.$file.gff\n";
    }
  }

  $log->write_to("Copying across GFF files from $dbpath/CHECKS/\n");
  system ("cp $dbpath/CHECKS/*.gff $www/$WS_name/GFF/") && croak "Could not copy GFF files from autoace/CHECKS $!\n";
}


###########################################################################################################
# Creates the two GFF_introns_confirmed_CDS* files using format of previous WS release, and wc -l of the
# GFF files in the GFF directory
###########################################################################################################

sub create_GFF_intron_files{
  
  $log->write_to("\ncreate_GFF_intron_files\n");
  $log->write_to("---------------------\n");

  $log->write_to("Creating new GFF_introns_confirmed_CDS... files in $www/$WS_name/Checks/\n");

 foreach my $lab ("cam", "stl"){
   my @line_counts = ();
   my $line_total = 0;

   foreach my $chrom (@chrom) {
     my $line = `wc -l $www/$WS_name/GFF/CHROMOSOME_$chrom.check_intron_$lab.gff`;
     # take line count and add to array
     my ($new) = ($line =~ /(\d+)/);
     push(@line_counts,$new);
     $line_total += $new;
   }

   # open old and new files
   open(OLDFILE, "<$www/$WS_previous_name/Checks/GFF_introns_confirmed_CDS_$lab.html") || croak "Couldn't read old $lab GFF intron file\n";
   open(NEWFILE, ">$www/$WS_name/Checks/GFF_introns_confirmed_CDS_$lab.html") || croak "Couldn't create new $lab GFF intron file\n";

   my $count = 0;
   # loop through old file replacing old info with new info from @line_counts 
   while(<OLDFILE>){
     if ((/<TD ALIGN=\"center\" COLSPAN=\"2\"><B>\s*(\d+)\s+\[\s*\d+\]<\/B><\/TD>/) && ($count <6)) {
	print NEWFILE "<TD ALIGN=\"center\" COLSPAN=\"2\"><B> ".$line_counts[$count]." [$1]</B></TD>\n";
	$count++;

     }
     elsif ((/<TD ALIGN=\"center\" COLSPAN=\"2\"><B>\s*(\d+)\s+\[\s*\d+\]<\/B><\/TD>/) && ($count > 5)) {
       print NEWFILE "<TD ALIGN=\"center\" COLSPAN=\"2\"><B> ".$line_total." [$1]</B></TD>\n";
       $count++;
     }
     else{
       print NEWFILE $_;
     }
   }
   close(OLDFILE);
   close(NEWFILE);
 }

}


###########################################################################################################
# Update the main wormpep pages at ~C_elegans/WORMBASE/WSXXX/                 
###########################################################################################################

sub update_wormpep_pages{
  
  $log->write_to("\nupdate_wormpep_pages\n");
  $log->write_to("--------------------\n");

  $log->write_to("Updating wormpep pages at $www/$WS_name\n");


  # first copy/create relevant files from previous release

  # make symbolic
  $log->write_to("Creating symbolic link for header.ini file\n");
  chdir("$www/$WS_name") || $log->write_to("Couldn't chdir to $www/$WS_name\n");
  system("ln -sf ../header.ini header.ini") && croak "Couldn't create new symlink\n";    
  
  # copy over static files that don't change
  $log->write_to("copying static wormpep files\n"); 
  system("cp -f $www/$WS_previous_name/wormpep.shtml $www/$WS_name/") && warn "Cannot copy wormpep.shtml $!\n";
  system("cp -f $www/$WS_previous_name/wormpep_changes.shtml $www/$WS_name/") && warn "Cannot copy wormpep_changes.shtml $!\n";
  system("cp -f $www/$WS_previous_name/wormpep_download.shtml $www/$WS_name/") && warn "Cannot copy wormpep_download.shtml $!\n";
  system("cp -f $www/$WS_previous_name/wormpep_format.shtml $www/$WS_name/") && warn "Cannot copy wormpep_format.shtml $!\n";
  system("cp -f $www/$WS_previous_name/wormpep_moreinfo.shtml $www/$WS_name/") && warn "Cannot copy wormpep_moreinfo.shtml $!\n";



  # write a new paragraph for the index.shtml page 
  undef $/;       
  open(LOGFILE,"<$basedir/WORMPEP/wormpep${WS_current}/wormpep_current.log") || croak "Couldn't open wormpep log file\n";
  my $text = <LOGFILE>;
  close(LOGFILE);
  $/ = "\n";
  $text =~ /No\. of sequences \(letters\) written:  (\d+,\d+)  \((\d+,\d+,\d+)\)/;

  # grab details of wormpep release, number of sequences etc
  my $count = $1;
  my $letters = $2;
  # calculate number of splice variants by looking for proteins ending in 'A' in wormpep.table file
  my $alt_spliced = `cut -f 1 $basedir/WORMPEP/wormpep${WS_current}/wormpep.table${WS_current} | sed 's/.*[0-9b-z]\$//' | grep \".\"| wc -l`; 
  $alt_spliced =~ s/\s+//g;

  # create release_paragraph.shtml
  open (PARAGRAPH, ">$www/$WS_name/release_paragraph.shtml") || croak "Can't create the file: $www/$WS_name/release_paragraph.shtml\n\n";
  print PARAGRAPH "The current Wormpep database, wormpep${WS_current} (released $release_date), contains $letters residues in $count protein sequences (of which $alt_spliced have splice variants) - wormpep${WS_current} is based on the <A href=\"ftp://ftp.sanger.ac.uk/pub/wormbase/WS$WS_current\">current WS$WS_current release</A> of the <I>C. elegans</I> AceDB database.\n";
  close (PARAGRAPH);

  # update the 'current_release.shtml' file
  open(RELEASE,">$www/$WS_name/current_release.shtml") || croak "Coudn't write to current_release.shtml\n";
  print RELEASE "The current release is wormpep${WS_current}\n";
  close(RELEASE);

  # update the history of wormpep releases, i.e. wormpep_release.txt
  # need to use file from previous release 
  open (TABLE, "<$www/$WS_previous_name/wormpep_release.txt") || croak "Can't read from the file: $www/$WS_name/wormpep_release.txt\n\n";
  open (TABLE2, ">$www/$WS_name/tmp_table.txt") || croak "Can't create the file: $www/$WS_name/tmp_table.txt\n\n";
  print TABLE2 "${WS_current}\t$release_date2\t$count\n";
  while(<TABLE>){
    print TABLE2 $_;
  }
  close(TABLE); 
  close(TABLE2);
  system("mv $www/$WS_name/tmp_table.txt $www/$WS_name/wormpep_release.txt") && croak "Couldn't rename file\n";


  # write a new table for the wormpep_download.shtml page    
  $log->write_to("Updating wormpep_download.shtml page\n");
  my $rows = $WS_current + 5;

  open (LIST, ">$www/$WS_name/releases.shtml") || croak "Can't open the file: $www/$WS_name/releases.shtml\n\n";
  print LIST "<CENTER><TABLE CELLPADDING=\"0\" CELLSPACING=\"0\" BORDER=\"0\">\n";
  print LIST "<TR VALIGN=\"top\">\n";
  print LIST "  <TD rowspan=\"$rows\" class=\"grey1\"><img src=\"/icons/blank.gif\" width=\"1\"   height=\"1\"></td>\n";
  print LIST "  <TD colspan=\"9\"  class=\"grey1\"><img src=\"/icons/blank.gif\" width=\"250\" height=\"1\"></td>\n";
  print LIST "  <TD rowspan=\"$rows\" class=\"grey1\"><img src=\"/icons/blank.gif\" width=\"1\"   height=\"1\"></td>\n";
  print LIST "</tr>\n";
  print LIST "<tr class=\"violet3\">\n";
  print LIST "    <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print LIST "    <TH>Release</TH>\n";
  print LIST "    <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print LIST "    <TH>Date</TH>\n";
  print LIST "    <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print LIST "    <TH>Proteins</TH>\n";
  print LIST "    <td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
  print LIST "</tr>\n";
  print LIST "<tr valign=\"top\">\n";
  print LIST "    <td colspan=\"9\" class=\"grey1\"><img src=\"/icons/blank.gif\" width=\"1\" height=\"1\"></td>\n";
  print LIST "</tr>\n";


  open (RETABLE, "<$www/$WS_name/wormpep_release.txt");
  my ($wp_rel,$wp_dat,$wp_CDS);
  while (<RETABLE>) {
    ($wp_rel,$wp_dat,$wp_CDS) = (/^(\S+)\s+(\S+)\s+(\S+)\s*/);
    print LIST "<!-- wormpep$wp_rel -->\n";
    
    if (($wp_rel % 2) == 0) { 
        print LIST "<tr class=\"violet2\">\n";
    }
    else {
        print LIST "<tr class=\"blue\">\n";
    }

    print LIST "<td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";

    if (($wp_rel < 8) || ($wp_rel == $WS_current)) {
        print LIST "<td>wormpep${wp_rel}</td>\n";
    }
    else {
        print LIST "<td><A href=\"ftp://ftp.sanger.ac.uk/pub/databases/wormpep/wormpep${wp_rel}\">wormpep${wp_rel}</A></td>\n";
    }

    print LIST "<td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
    print LIST "<td>$wp_dat</td>\n";
    print LIST "<td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
    print LIST "<td>$wp_CDS</td>\n";
    print LIST "<td><img src=\"/icons/blank.gif\" width=\"16\" height=\"22\"></td>\n";
    print LIST "</tr> \n\n";
  }
  close(RETABLE); 
  print LIST "<tr valign=\"top\">\n";
  print LIST "  <td colspan=\"9\" class=\"grey1\"><img src=\"/icons/blank.gif\" width=\"250\" height=\"1\"></td>\n";
  print LIST "</tr>\n";
  print LIST "</TABLE></CENTER>\n";
  close (LIST);


}


