#!/usr/local/bin/perl5.8.0 -w
#
# update_website.pl
# 
# by Keith Bradnam aged 12 and half (is this the reincarnation of Peter Pan?)
#
# A script to finish the last part of the weekly build by updating all of the
# relevant WormBase and Wormpep web pages.
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2009-12-22 13:31:50 $      


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
use File::Basename;
use IO::File;
use Switch;
use threads;
use Thread::Queue;

##############################
# command-line options       #
##############################

my $debug;              # Verbose debug mode careful: only chanegs source dirs
my $help;               # Help/Usage page
my $all;                # option all - do everything
my $header;             # create release_header.shtml
my $dna;                # create DNA.shtml
my $finished;           # create finished.shtml
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
            "wormpepdiff"    => \$wormpep_diff,
	    "copygff"        => \$copyGFF,
            "debug=s"        => \$debug,
	    "test"           => \$test,
            "help"           => \$help,
            "h"              => \$help,	# 
	    "dirs"           => \$directories,
	    "overlap"        => \$overlap,
	    "est"            => \$EST_files,
	    "create_gff"     => \$create_GFF,
	    "wormpep"        => \$update_wormpep,
	    'store=s'	     => \$store
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

###########################################################################################################
# Main subroutine calls    
###########################################################################################################


&create_web_directories        if ($all || $directories);
&create_release_header         if ($all || $header);    # makes release_header.shtml
&create_DNA_table              if ($all || $dna);       # makes DNA.shtml
&create_wormpep_page           if ($all || $wormpep_diff); # the wormpep page accessible from the main wormbase page
#&copy_overlapcheck_files       if ($all || $overlap); # copies check files/wormsrv2/autoace/CHECKS and converts to HTML etc.
#&copy_EST_files                if ($all || $EST_files); # copies EST_analysis.html and then edits this file then copies other EST files
#&copy_GFF_files                if ($all || $copyGFF); # copies some of the GFF files in /wormsrv2/autoace/GFF_SPLITS
#&create_GFF_intron_files       if ($all || $create_GFF); # copies the GFF_introns_confirmed_CDS* files from previous WS release, updates totals
&update_wormpep_pages          if ($all || $update_wormpep); # update the main wormpep web pages



#########################
# final bit of tidying up
#########################

##################
# Check the files
##################

if ($all || $directories) {
  $wb->check_file("$www/$WS_name", $log);
  $wb->check_file("$www/$WS_name/GFF", $log);
}
if ($all || $header) {
  $wb->check_file("$www/$WS_name/release_header.shtml", $log, 
		  minsize => 42, 
		  maxsize => 42);
}
if ($all || $dna) {
  $wb->check_file("$www/$WS_name/DNA_table.shtml", $log, 
		  minsize => 3200,
		  maxsize => 3300);
}


if ($all || $wormpep_diff) {
  $wb->check_file("$www/$WS_name/wormpep_diff.shtml", $log, 
		  minsize => 1000,
		  maxsize => 10000);
}

if ($all || $update_wormpep) {
  $wb->check_file("$www/$WS_name/release_paragraph.shtml", $log,
		  minsize => 10,);		
  $wb->check_file("$www/$WS_name/current_release.shtml", $log,
		  minsize => 33,		
		  maxsize => 36,);		
  $wb->check_file("$www/$WS_name/wormpep_release.txt", $log,
		  minsize => 3490,);		
  $wb->check_file("$www/$WS_name/releases.shtml", $log,
		  minsize => 70000,);		
  $wb->check_file("$www/$WS_name/wormpep.shtml", $log,
		  samesize => "$www/$WS_previous_name/wormpep.shtml");
  $wb->check_file("$www/$WS_name/wormpep_changes.shtml", $log,
		  samesize => "$www/$WS_previous_name/wormpep_changes.shtml");
  $wb->check_file("$www/$WS_name/wormpep_download.shtml", $log,
		  samesize => "$www/$WS_previous_name/wormpep_download.shtml");
  $wb->check_file("$www/$WS_name/wormpep_format.shtml", $log,
		  samesize => "$www/$WS_previous_name/wormpep_format.shtml");
  $wb->check_file("$www/$WS_name/wormpep_moreinfo.shtml", $log,
		  samesize => "$www/$WS_previous_name/wormpep_moreinfo.shtml");
  
}


$log->write_to("\n\nFinished\n");

# mail log
$log->mail( "$maintainers", "BUILD REPORT: $0" );
exit(0);


###########################################################################################################
# Create new directories on website     
###########################################################################################################
sub create_web_directories{

  if (! -e "$www/$WS_name"){
    system("mkdir $www/$WS_name") && warn "Cannot create $WS_name directory $!\n";
    system("mkdir $www/$WS_name/GFF") && warn "Cannot create $WS_name/GFF directory $!\n";
  }

}




###########################################################################################################
# Creates: release_header.shtml, DNA_table.shtml and wormpep.shtml
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
    $DNA_tot = $1 if (/(\d+) total/);
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
  
  $log->write_to("Opening log file ".$wb->wormpep."/wormpep_current.log'\n");

  my ($wp_seq,$wp_let);
 open (WP_1, "<".$wb->wormpep."/wormpep_current.log") || croak "Failed to open WORMPEP/wormpep_current.log\n";
  while (<WP_1>) {
    # No. of sequences (letters) written:  22,221  (9,696,145)
    if (/No\. of sequences \(letters\) written:  (\d+,\d+)  \((\d+,\d+,\d+)\)/) {
      ($wp_seq,$wp_let) = ($1,$2);
    }
  }
  close (WP_1);

  
  # get details from file prepared for reports.
  my $wp_alt;
  my $wormpeplog  = $wb->reports."/wormpep";
  open (WP_2, "<$wormpeplog") || croak "Failed to open REPORTS/wormpep\n:$!";
  while (<WP_2>) {
    if (/(\d+) alternate splice forms/) {
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
  open (WP_3, "<".$wb->wormpep."/wormpep.diff${WS_current}") or $log->write_to("Failed to open wormpep.diff\n");

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
  print PARAGRAPH "The current Wormpep database, wormpep${WS_current} (released $release_date), contains $letters residues in $count protein sequences (of which $alt_spliced have splice variants) - wormpep${WS_current} is based on the <A href=\"ftp://ftp.sanger.ac.uk/pub2/wormbase/WS$WS_current\">current WS$WS_current release</A> of the <I>C. elegans</I> AceDB database.\n";
  close (PARAGRAPH);

  # update the 'current_release.shtml' file
  open(RELEASE,">$www/$WS_name/current_release.shtml") || croak "Coudn't write to current_release.shtml\n";
  print RELEASE "The current release is wormpep${WS_current}\n";
  close(RELEASE);

  # update the history of wormpep releases, i.e. wormpep_release.txt
  # need to use file from previous release 
  open (TABLE, "<$www/$WS_previous_name/wormpep_release.txt") || carp "Can't read from the file: $www/$WS_name/wormpep_release.txt\n\n";
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


