#!/usr/local/bin/perl5.8.0 -w
#
# update_website.pl
# 
# by Keith Bradnam aged 12 and half (is this the reincarnation of Peter Pan?)
#
# A script to finish the last part of the weekly build by updating all of the
# relevant WormBase and Wormpep web pages.
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2004-06-18 14:49:15 $      


#################################################################################
# load modules etc.                                                             #
#################################################################################

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Cwd;
use Symbol 'gensym';
use Ace;
use lib "/nfs/WWW/perl";
use SangerWeb;
use Carp;

##############################
# command-line options       #
##############################


my $debug;              # Verbose debug mode
my $help;               # Help/Usage page
my $all;                # option all - do everything
my $header;             # create release_header.shtml
my $dna;                # create DNA.shtml
my $finished;           # create finished.shtml
my $dbcomp;             # create dbcomp.shtml
my $wormpep_diff;       # diff log for wormpep releases
my $copyGFF;            # copy GFF files to WWW site 
my $test;               # test run (requires a release number)
my $directories;        # create website directories
my $overlap;            # run copy_overlapcheck_files
my $EST_files;          # run copy_EST_files
my $create_GFF;         # run create_GFF_intron_files
my $update_wormpep;


GetOptions (
            "all"            => \$all,
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
	    "wormpep"        => \$update_wormpep
	   );


##############################
# Script variables           #
##############################

my $maintainers      = "All";
my $rundate          = `date +%y%m%d`; chomp $rundate;
my $runtime          = `date +%H:%M:%S`; chomp $runtime;
my @chrom            = qw ( I II III IV V X ); 
my $WS_current       = &get_wormbase_version;
my $release_date     = &get_wormbase_release_date("long");
my $release_date2    = &get_wormbase_release_date("short");
my $WS_previous      = $WS_current - 1;
my $WS_name          = &get_wormbase_version_name;
my $WS_previous_name = "WS".$WS_previous;


# file path info
my $www_root         = "/nfs/WWWdev/htdocs/Projects/C_elegans";
my $www              = "/nfs/WWWdev/htdocs/Projects/C_elegans/WORMBASE";
my $wwwdata          = "/nfs/WWWdev/SANGER_docs/data/Projects/C_elegans";
my $wwwlive          = "/nfs/WWW/htdocs/Projects/C_elegans/WORMBASE";
my $gff              = "/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS";
my $dbpath           = "/wormsrv2/autoace";
my $log             = "/wormsrv2/logs/update_website.WS${WS_current}.$rundate.$$";


#############
# TEST MODE
############

if (defined $test) {
    $WS_current       = $test;
    $WS_previous      = $test - 1;
    $WS_name          = "WS".$WS_current;
    $WS_previous_name = "WS".$WS_previous;
}

# debug mode modifies $maintainers to reduce e-mail load
($maintainers = $debug . "\@sanger.ac.uk") if ($debug);

# check for command-line options if none given then you do everything
unless (defined $ARGV[0]) {
    $all = 1;
}

###########################################################################################################
# Main subroutine calls    
###########################################################################################################

&create_log_file               if ($all);

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
print LOG "\nChanging 'current symbolic link to point to new release\n"; 	
chdir("$www") || print LOG "Couldn't chdir to $www\n";
system("rm -f $www/current") && croak "Couldn't remove 'current' symlink\n"; 	 
system("ln -s $WS_name current") && croak "Couldn't create new symlink\n";

print LOG "\n\nC'est finis\n\n";

# mail log
&mail_maintainer("WormBase Report: update_website.pl",$maintainers,$log);

close(LOG);

exit(0);


###########################################################################################################
# Create log file          
###########################################################################################################

sub create_log_file{
    
  open (LOG,">$log") || croak "Cannot open logfile $!\n";
  LOG->autoflush();
  
  print LOG "# update_website.pl\n\n";     
  print LOG "# run details    : $rundate $runtime\n";
  print LOG "# WormBase version : $WS_name\n";
  print LOG "\n\n";

}




###########################################################################################################
# Create new directories on website     
###########################################################################################################

sub create_web_directories{

  if (! -e "$www".$WS_name){
    system("mkdir $www/$WS_name") && warn "Cannot create $WS_name directory $!\n";
    system("mkdir $www/$WS_name/Checks") && warn "Cannot create $WS_name/Checks directory $!\n";
    system("mkdir $www/$WS_name/GFF") && warn "Cannot create $WS_name/GFF directory $!\n";
  }

}



###########################################################################################################
# Copy overlapcheck files to website    
###########################################################################################################

sub copy_overlapcheck_files{

  print LOG "\ncopy_overlapcheck_files\n";
  print LOG "-----------------------\n";

  # list of files to be copied
  my @filenames = qw( overlapping_genes_cam overlapping_genes_stl EST_in_intron_cam EST_in_intron_stl repeat_in_exon_cam repeat_in_exon_stl );

  print LOG "copying files from /wormsrv2/autoace/CHECKS/ to $www/$WS_name/Checks\n"; 
  system("cp -f $www/$WS_previous_name/Checks/index.shtml $www/$WS_name/Checks/") && warn "Cannot copy index.shtml $!\n";

  print LOG "copying three ace files from /wormsrv2/autoace/CHECKS/ to $www/$WS_name/Checks\n"; 
  system("cp -f $www/$WS_previous_name/Checks/*.ace $www/$WS_name/Checks/") && warn "Cannot copy *.ace files $!\n";

  print LOG "copying short genes file /wormsrv2/autoace/CHECKS/ to $www/$WS_name/Checks\n"; 
  system("cp -f /wormsrv2/autoace/CHECKS/short_spurious_genes.$WS_name.csv.html $www/$WS_name/Checks/short_genes.html") && warn "Cannot copy short_genes file $!\n";


  foreach my $file (@filenames) {
    my $line_total;
    my @line_stats;


    foreach my $chrom (@chrom) {
      print LOG "Copying file CHROMOSOME_$chrom.$file to $www/$WS_name/Checks\n";
      system ("cp -f /wormsrv2/autoace/CHECKS/CHROMOSOME_$chrom.$file $www/$WS_name/Checks/")
	&& croak "Could not copy /wormsrv2/autoace/CHECKS/CHROMOSOME_$chrom.$file $!\n";
      
      print LOG "Calculating line numbers for CHROMOSOME_$chrom.$file\n";
      my $line = `wc -l /wormsrv2/autoace/CHECKS/CHROMOSOME_$chrom.$file`;
      # take line count and add to array
      my ($new) = ($line =~ /(\d+)/);
      push @line_stats, $new;
      # add to line total for all chromosomes
      $line_total += $new;
    }
    
    
    # make new html files
    print LOG "Generating new html files in $www/$WS_name/Checks/\n";
    my $fh     = gensym();
    my $newfh  = gensym();
    my $count = 0;
    
    open ($fh, "$www/$WS_previous_name/Checks/$file.html") || croak "Cannot open old html file $!\n";
    open ($newfh, ">$www/$WS_name/Checks/$file.html") || croak "Cannot open new html file $!\n";
    
    print LOG "Generating $www/$WS_name/Checks/$file.html\n";
    
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
  print LOG "\ncreate_top_level_web_pages\n";
  print LOG "--------------------------\n";
  
  print LOG "Creating $www/$WS_name/release_header.shtml\n";
  open  (RELEASE_HEADER, ">$www/$WS_name/release_header.shtml") || croak "Couldn't create release_header.shtml\n";
  print  RELEASE_HEADER  "<b>$WS_name release, updated $release_date.</b>\n";
  close (RELEASE_HEADER);

}

sub create_DNA_table {

  # create DNA.table.shtml
  print LOG "Creating $www/$WS_name/DNA.table.shtml\n";

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
  
  open (COMPOSITION, "</wormsrv2/autoace/CHROMOSOMES/composition.all") || croak "Failed to open composition.all\n";
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
  print LOG "Creating $www/$WS_name/finished_seqs.shtml\n";
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
  print LOG "Creating $www/$WS_name/dbcomp.shtml\n";

  system ("cp /wormsrv2/autoace/COMPARE/current.out $wwwdata/WS.dbcomp_output") && croak "Couldn't copy current.out\n";
  my $dbcomplen;
  
  open (DBCOMP, "wc -l /wormsrv2/autoace/COMPARE/current.dbcomp |");
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

  open (COMP, "</wormsrv2/autoace/COMPARE/current.dbcomp") || croak "Failed to open composition.all\n\n";
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
  print LOG "\ncreate_wormpep_page\n";
  print LOG "-------------------\n";


  print LOG "Creating $www/$WS_name/wormpep.shtml\n";

  open (WORMPEP, ">$www/$WS_name/wormpep.shtml") || croak "Failed to create wormpep.shtml\n\n";
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

  print LOG "Opening log file '/wormsrv2/WORMPEP/wormpep$WS_current/wormpep_current.log'\n";

  my ($wp_seq,$wp_let);
  open (WP_1, "</wormsrv2/WORMPEP/wormpep$WS_current/wormpep_current.log") || croak "Failed to open wormpep.log\n";
  while (<WP_1>) {
    # No. of sequences (letters) written:  22,221  (9,696,145)
    if (/No\. of sequences \(letters\) written:  (\d+,\d+)  \((\d+,\d+,\d+)\)/) {
      ($wp_seq,$wp_let) = ($1,$2);
    }
  }
  close (WP_1);

  
  # get details from last wormpep log file
  my $wp_alt;
  my @possible_logs = `ls -t /wormsrv2/logs/make_wormpep.final.WS${WS_current}*`; # added new logfile name 
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
  open (WP_3, "</wormsrv2/WORMPEP/wormpep${WS_current}/wormpep.diff${WS_current}") || croak "Failed to open wormpep.history\n";
  while (<WP_3>) {
    (push (@changed,$_)) if (/changed:/);
    (push (@lost,$_)) if (/lost:/);
    (push (@new,$_)) if (/new:/);
    (push (@reappeared,$_)) if (/reappeared:/);
    #       print WORMPEP "$_";
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

sub copy_EST_files{

  # make EST_*.shtml using EST_*.txt in /wormsrv2/autoace/CHECKS  
  print LOG "\ncopy_EST_files\n";
  print LOG "--------------\n";


  my %files = ("EST_total" => "Main report",
               "EST_no_accession" => "No Accession number", 
               "EST_unassigned" => "Unassigned reads",
               "EST_mismatched" => "mismatched CDS assignments");

  my @line_counts;
  
  foreach my $file (keys(%files)){
    print LOG "creating $file.shtml in $www/$WS_name/Checks\n"; 
    open(EST_HTML,">$www/$WS_name/Checks/$file.shtml") || croak "Couldn't create $file.shtml\n";
    
    print EST_HTML &SangerWeb::virtual_header();
    
    print EST_HTML "<TABLE WIDTH=\"100%\" CELLSPACING=\"0\" CELLPADDING=\"0\"><TR VALIGN=\"top\" BGCOLOR=\"darkblue\" ><TD WIDTH=\"100%\"><BR><H2 align=\"center\"><FONT COLOR=\"white\">";
    print EST_HTML "EST analysis: \"$files{$file}\"";
    print EST_HTML "</FONT></H2></TD></TR></TABLE>\n";

    print EST_HTML "<P><PRE>";
    open (EST_TXT, "</wormsrv2/autoace/CHECKS/$file.txt") || croak "Couldn't open EST_total.txt\n";

    my $line = `wc -l /wormsrv2/autoace/CHECKS/$file.txt`;
    # take line count and add to array
    my ($new) = ($line =~ /(\d+)/);
    push(@line_counts,$new);

    while (<EST_TXT>) {
     print EST_HTML "$_";
    }
    print EST_HTML "</PRE></P>";
    print EST_HTML &SangerWeb::virtual_footer();
    close(EST_TXT);
    close(EST_HTML);
  }

  # need to add 'n/a' into the line_counts array to avoid problems later on
  splice(@line_counts,2,0,"n/a");


  print LOG "updating EST_analysis.shtml from $www/$WS_previous/Checks to $www/$WS_name/Checks, adding new info\n"; 

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
# Copy some of the GFF files from /wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS/ and some from CHECKS/
###########################################################################################################

sub copy_GFF_files{
  
  print LOG "\ncopy_GFF_files\n";
  print LOG "--------------\n";

  print LOG "Copying across GFF files from /wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS\n";

  #simple double foreach loop to loop through each chromosome and file name
  my @gff_files = ("clone_ends", "clone_path", "exon","intron_all", "clone_acc", "genes", "repeats", "intron", "rna", "UTR");
  foreach my $chrom (@chrom) {
    foreach my $file (@gff_files){
      system("sort -u $gff/CHROMOSOME_$chrom.$file.gff | gff_sort > $www/$WS_name/GFF/CHROMOSOME_$chrom.$file.gff")
	&& croak "Couldn't copy CHROMOSOME_$chrom.$file.gff\n";
    }
  }

  print LOG "Copying across GFF files from /wormsrv2/autoace/CHECKS/\n";
  system ("cp /wormsrv2/autoace/CHECKS/*.gff $www/$WS_name/GFF/") && croak "Could not copy GFF files from autoace/CHECKS $!\n";
}


###########################################################################################################
# Creates the two GFF_introns_confirmed_CDS* files using format of previous WS release, and wc -l of the
# GFF files in the GFF directory
###########################################################################################################

sub create_GFF_intron_files{
  
  print LOG "\ncreate_GFF_intron_files\n";
  print LOG "---------------------\n";

  print LOG "Creating new GFF_introns_confirmed_CDS... files in $www/$WS_name/Checks/\n";

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
   open(NEWFILE, ">$www/$WS_name/Checks/GFF_introns_$lab.html") || croak "Couldn't create new $lab GFF intron file\n";

   my $count = 0;
   # loop through old file replacing old info with new info from @line_counts 
   while(<OLDFILE>){
     if ((/<TD ALIGN=\"center\" COLSPAN=\"2\"><B>\s*(\d+)\s+\[\s*\d+\]<\/B><\/TD>/) && ($count <6)) {
	my $old = $1;
	print NEWFILE "<TD ALIGN=\"center\" COLSPAN=\"2\"><B> ".$line_counts[$count]." [$old]</B></TD>\n";
	$count++;

     }
     elsif ((/<TD ALIGN=\"center\" COLSPAN=\"2\"><B>\s*(\d+)\s+\[\s*\d+\]<\/B><\/TD>/) && ($count > 5)) {
       my $old = $1;       
       print NEWFILE "<TD ALIGN=\"center\" COLSPAN=\"2\"><B> ".$line_total." [$old]</B></TD>\n";
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
# Update the main wormpep pages at ~C_elegans/wormpep                 
###########################################################################################################

sub update_wormpep_pages{
  
  print LOG "\nupdate_wormpep_pages\n";
  print LOG "--------------------\n";

  print LOG "Updating wormpep pages at $www_root/wormpep\n";
  
  # write a new paragraph for the index.shtml page 
  undef $/;       
  open(LOGFILE,"</wormsrv2/WORMPEP/wormpep${WS_current}/wormpep_current.log") || croak "Couldn't open wormpep log file\n";
  my $text = <LOGFILE>;
  close(LOGFILE);
  $/ = "\n";
  $text =~ /No\. of sequences \(letters\) written:  (\d+,\d+)  \((\d+,\d+,\d+)\)/;

  # grab details of wormpep release, number of sequences etc
  my $count = $1;
  my $letters = $2;
  # calculate number of splice variants by looking for proteins ending in 'A' in wormpep.table file
  my $alt_spliced = `cut -f 1 /wormsrv2/WORMPEP/wormpep${WS_current}/wormpep.table${WS_current} | sed 's/.*[0-9b-z]\$//' | grep \".\"| wc -l`; 
  $alt_spliced =~ s/\s+//g;

  # create release_paragraph.shtml
  system("rm -f $www_root/wormpep/release_paragraph.shtml") && croak "Couldn't remove old release_paragraph.shtml\n";
  open (PARAGRAPH, ">$www_root/wormpep/release_paragraph.shtml") || croak "Can't create the file: $www_root/wormpep/release_paragraph.shtml\n\n";
  print PARAGRAPH "The current Wormpep database, wormpep${WS_current} (released $release_date), contains $letters residues in $count protein sequences (of which $alt_spliced have splice variants) - wormpep${WS_current} is based on the <A href=\"ftp://ftp.sanger.ac.uk/pub/wormbase/WS$WS_current\">current WS$WS_current release</A> of the <I>C. elegans</I> AceDB database.\n";
  close (PARAGRAPH);

  # update the 'current_release.shtml' file
  system("rm -f $www_root/wormpep/current_release.shtml") && croak "Couldn't remove current_release.shtml\n";
  open(RELEASE,">$www_root/wormpep/current_release.shtml") || croak "Coudn't write to current_release.shtml\n";
  print RELEASE "The current release is wormpep${WS_current}\n";
  close(RELEASE);

  # update the history of wormpep releases, i.e. wormpep_release.txt

  open (TABLE, "<$www_root/wormpep/wormpep_release.txt") || croak "Can't read from the file: $www_root/wormpep/wormpep_release.txt\n\n";
  open (TABLE2, ">$www_root/wormpep/tmp_table.txt") || croak "Can't create the file: $www_root/wormpep/tmp_table.txt\n\n";
  print TABLE2 "${WS_current}\t$release_date2\t$count\n";
  while(<TABLE>){
    print TABLE2 $_;
  }
  close(TABLE); 
  close(TABLE2);
  system("mv $www_root/wormpep/tmp_table.txt $www_root/wormpep/wormpep_release.txt") && croak "Couldn't rename file\n";


  # write a new table for the wormpep_download.shtml page    
  print LOG "Updating wormpep_download.shtml page\n";
  my $rows = $WS_current + 5;

  open (LIST, ">$www_root/wormpep/releases.shtml") || croak "Can't open the file: $www_root/wormpep/releases.shtml\n\n";
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


  open (RETABLE, "<$www_root/wormpep/wormpep_release.txt");
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
        print LIST "<td><A href=\"ftp://ftp.sanger.ac.uk/pub/databases/wormpep/old_wormpep${wp_rel}\">wormpep${wp_rel}</A></td>\n";
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


