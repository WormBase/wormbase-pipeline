#!/usr/local/bin/perl5.6.0
#
# unpdate_website.pl
# 
# by Keith Bradnam aged 12 and half
#
# Usage : update_website.pl
#
# A script to finish the last part of the weekly build by updating all of the
# relevant WormBase and Wormpep web pages.

#################################################################################
# load modules etc.                                                             #
#################################################################################

$|=1;
use IO::Handle;
use Getopt::Std;
use Cwd;
use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Symbol 'gensym';
use Ace;
use lib "/nfs/WWW/perl";
use SangerWeb;

##############################
# Script variables           #
##############################

my $maintainers      = "dl1\@sanger.ac.uk krb\@sanger.ac.uk kj2\@sanger.ac.uk";
my $rundate          = `date +%y%m%d`; chomp $rundate;
my $runtime          = `date +%H:%M:%S`; chomp $runtime;
my $cvs_version      = &get_cvs_version("$0");
my @chrom            = qw ( I II III IV V X ); 
my $WS_current       = &get_wormbase_version;
my $release_date     = &get_wormbase_release_date("long");
my $release_date2    = &get_wormbase_release_date("short");
my $WS_previous      = $WS_current - 1;
my $WS_name          = &get_wormbase_version_name;
my $WS_previous_name = "WS".$WS_previous;
my $wormpep_version  = &get_wormpep_version;

# file path info
my $www_root         = "/nfs/WWW/htdocs/Projects/C_elegans";
my $www              = "/nfs/WWW/htdocs/Projects/C_elegans/WORMBASE";
my $gff              = "/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS";
my $dbpath           = "/wormsrv2/autoace";
our $log             = "/wormsrv2/logs/update_website.$rundate";


###########################################################################################################
# Main subroutine calls    
###########################################################################################################

&create_log_file;

&create_web_directories;

&create_top_level_web_pages; # e.g. release_header.shtml, finished_seqs.shtml, unfinished_seqs.shtml etc.

&create_wormpep_page; # the wormpep page accessible from the main wormbase page

&copy_overlapcheck_files; # copies files created in /wormsrv2/autoace/CHECKS and converts to HTML etc.

&copy_EST_files; # copies EST_analysis.html and then edits this file then copies other EST files

&copy_GFF_files; # copies some of the GFF files in /wormsrv2/autoace/GFF_SPLITS

&create_GFF_intron_files; # copies the GFF_introns_confirmed_CDS* files from previous WS release, updates totals

&update_wormpep_pages; # update the main wormpep web pages

&create_loci_designations; # previously done by the list_loci_designations script




#########################
# final bit of tidying up
#########################

# update 'current' symlink
print LOG "\nChanging 'current symbolic link to point to new release\n";
system("rm -f $www/current") && die "Couldn't remove 'current' symlink\n";
system("ln -s $www/$WS_name/ current") && die "Couldn't create new symlink\n";

print LOG "\n\nC'est finis\n\n";

# mail log
&mail_maintainer("WormBase Report: unpdate_website.pl",$maintainers,$log);

close(LOG);
exit(0);






###########################################################################################################
# Create log file          
###########################################################################################################

sub create_log_file{
    
  open (LOG,">$log") || die "Cannot open logfile $!\n";
  LOG->autoflush();
  
  print LOG "# update_website.pl\n\n";     
  print LOG "# run details    : $rundate $runtime\n";
  print LOG "# WormBase version : $WS_name\n";
  print LOG "# cvs version      : $cvs_version\n";
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


  foreach my $file (@filenames) {
    my $line_total;
    my @line_stats;


    foreach my $chrom (@chrom) {
      print LOG "Copying file CHROMOSOME_$chrom.$file to $www/$WS_name/Checks\n";
      system ("cp -f /wormsrv2/autoace/CHECKS/CHROMOSOME_$chrom.$file $www/$WS_name/Checks/")
	&& die "Could not copy /wormsrv2/autoace/CHECKS/CHROMOSOME_$chrom.$file $!\n";
      
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
    
    open ($fh, "$www/$WS_previous_name/Checks/$file.html") || die "Cannot open old html file $!\n";
    open ($newfh, ">$www/$WS_name/Checks/$file.html") || die "Cannot open new html file $!\n";
    
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
    close($fh)    || die "Cannot close old html file $!\n";
    close($newfh) || die "Cannot close new html file $!\n";
  }
}


###########################################################################################################
# Creates: release_header.shtml, finished_seqs.shtml, unfinished_seqs.shtml, DNA_table.shtml,
# dbcomp.shtml, and wormpep.shtml
###########################################################################################################

sub create_top_level_web_pages{

  # create release_header.shtml  
  print LOG "\ncreate_top_level_web_pages\n";
  print LOG "--------------------------\n";
  

  print LOG "Creating $www/$WS_name/release_header.shtml\n";
  open (RELEASE_HEADER, ">$www/$WS_name/release_header.shtml") || die "Couldn't create release_header.shtml\n";
  print RELEASE_HEADER  "<P>\n";
  print RELEASE_HEADER  "<FONT SIZE=+2 COLOR=purple><B>$WS_name release</B></FONT><BR>\n";
  print RELEASE_HEADER  "<FONT SIZE=-2 COLOR=red><B>Last updated $release_date.</B></FONT>\n";
  print RELEASE_HEADER  "</P>\n";
  close (RELEASE_HEADER);



  # create DNA.table.shtml
  print LOG "Creating $www/$WS_name/DNA.table.shtml\n";
  my ($DNA_tot, $DNA_d, $DNA_dp, $DNA_n, $DNA_np, $DNA_g, $DNA_gp, $DNA_a, $DNA_ap, $DNA_t, $DNA_tp, $DNA_c, $DNA_cp);
  open (DNA_table, ">$www/$WS_name/DNA_table.shtml") || die "Couldn't create DNA_table.shtml\n\n";
  print DNA_table "<P>\n";
  
  open (COMPOSITION, "</wormsrv2/autoace/CHROMOSOMES/composition.all") || die "Failed to open composition.all\n";
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
    
  print DNA_table "<SPACER TYPE=\"horizontal\" SIZE=\"50\">\n";
  print DNA_table "<TABLE WIDTH=\"40%\" CELLSPACING=\"0\" CELLWIDTH=\"0\" BORDER=\"1\"><TR BGCOLOR=\"yellow\"><TH>&nbsp</TH> <TH>No.</TH> <TH>%</TH></TR>\n";
  print DNA_table "<TR> <TD ALIGN=\"CENTER\">Total</TD> <TD ALIGN=\"RIGHT\">$DNA_tot</TD> <TD ALIGN=\"RIGHT\">100</TD>     </TR>\n";
  print DNA_table "<TR> <TD ALIGN=\"CENTER\">-</TD>     <TD ALIGN=\"RIGHT\">$DNA_d</TD>   <TD ALIGN=\"RIGHT\">$DNA_dp</TD> </TR>\n";
  print DNA_table "<TR> <TD ALIGN=\"CENTER\">G</TD>     <TD ALIGN=\"RIGHT\">$DNA_g</TD>   <TD ALIGN=\"RIGHT\">$DNA_gp</TD> </TR>\n";
  print DNA_table "<TR> <TD ALIGN=\"CENTER\">A</TD>     <TD ALIGN=\"RIGHT\">$DNA_a</TD>   <TD ALIGN=\"RIGHT\">$DNA_ap</TD> </TR>\n";
  print DNA_table "<TR> <TD ALIGN=\"CENTER\">T</TD>     <TD ALIGN=\"RIGHT\">$DNA_t</TD>   <TD ALIGN=\"RIGHT\">$DNA_tp</TD> </TR>\n";
  print DNA_table "<TR> <TD ALIGN=\"CENTER\">C</TD>     <TD ALIGN=\"RIGHT\">$DNA_c</TD>   <TD ALIGN=\"RIGHT\">$DNA_cp</TD> </TR>\n";
  print DNA_table "<TR> <TD ALIGN=\"CENTER\">N</TD>     <TD ALIGN=\"RIGHT\">$DNA_n</TD>   <TD ALIGN=\"RIGHT\">$DNA_np</TD> </TR>\n";
  print DNA_table "</TABLE>\n";
  print DNA_table "</P>\n";
  close (DNA_table);


  # create finished_seqs.shtml
  print LOG "Creating $www/$WS_name/finished_seqs.shtml\n";
  my $db = Ace->connect(-path=>$dbpath) || do { print "Connection failure: ",Ace->error; die();};
  my $count = $db->fetch(-query=> 'find genome_sequence where Finished');

  open (FINISHED_SEQS, ">$www/$WS_name/finished_seqs.shtml") || die "Couldn't create finished_seqs.shtml\n";
  print FINISHED_SEQS  "<P>\n";
  print FINISHED_SEQS "$count projects\n";
  print FINISHED_SEQS  "</P>\n";
  close (FINISHED_SEQS);



  # create unfinished_seqs.shtml
  print LOG "Creating $www/$WS_name/unfinished_seqs.shtml\n";
  my %who = ('HX', 'Hinxton', 'RW', 'St Louis');

  $count = $db->fetch(-query=> 'find genome_sequence where Remark="Unfinished fragment"');

  open (UNFINISHED_SEQS, ">$www/$WS_name/unfinished_seqs.shtml") || die "Couldn't creat unfinished_seqs.shtml\n";
  print UNFINISHED_SEQS  "<P>\n";
  print UNFINISHED_SEQS "$count projects\n";
  print UNFINISHED_SEQS  "</P>\n";
  print UNFINISHED_SEQS "<SPACER TYPE=\"horizontal\" SIZE=\"50\">\n";
  print UNFINISHED_SEQS "<TABLE WIDTH=\"40%\" CELLSPACING=\"0\" CELLWIDTH=\"0\" BORDER=\"1\"><TR BGCOLOR=\"yellow\" ><TH>Project</TH> <TH>Source</TH> <TH>EMBL_Acc</TH> <TH>Length (bp)</TH></TR>\n";
  my $i = $db->fetch_many(-query=> 'find genome_sequence where Remark="Unfinished fragment"');  
  while (my $obj = $i->next) {
    my $project = $obj;
        
    # Database information
    my $dbaccno = $obj->at('DB_info.Database[3]');
    if (!$dbaccno) {
      print "## Warning: No accession no\n";
      $dbaccno = "n/a";
    }
        
    # DNA length
    my $DNA_len = $obj->DNA(2);
    
    # Laboratory
    my $lab = $obj->From_Laboratory(1);
    
    print UNFINISHED_SEQS "<TR> <TD ALIGN=\"CENTER\">$project</TD> <TD ALIGN=\"CENTER\">$who{$lab}</TD> <TD ALIGN=\"RIGHT\">$dbaccno</TD> <TD ALIGN=\"RIGHT\">$DNA_len</TD> </TR>\n";
        
  }
  print UNFINISHED_SEQS "</TABLE>\n";
  close (UNFINISHED_SEQS);




  # create dbcomp.shtml
  print LOG "Creating $www/$WS_name/dbcomp.shtml\n";

  system ("cp /wormsrv2/autoace/COMPARE/current.out $www/$WS_name/WS.dbcomp_output") && die "Couldn't copy current.out\n";

  open (DB_comp, ">$www/$WS_name/dbcomp.shtml") || die "Couldn't create dbcomp.shtml\n\n";
  print DB_comp "<P>\n";
  print DB_comp "<SPACER TYPE=\"horizontal\" SIZE=\"50\">\n";
  print DB_comp "<TABLE WIDTH=\"40%\" CELLSPACING=\"0\" CELLWIDTH=\"0\" BORDER=\"1\"><TR BGCOLOR=\"yellow\" ><TH>Class</TH> <TH>$WS_previous_name</TH> <TH>$WS_name</TH> <TH>Diff</TH></TR>\n";

  open (COMP, "</wormsrv2/autoace/COMPARE/current.dbcomp") || die "Failed to open composition.all\n\n";
  while (<COMP>) {
    next if (/\+/);
    next if (/ACEDB/);
    next if (/change/);

    my ($nowt,$class,$count_1,$count_2,$diff) = split (/\|/,$_);
    $class =~ s/^\s+//g;
    print DB_comp "<TR> <TD ALIGN=\"RIGHT\"><A  href=\"/cgi-bin/Protozoa/wormbase_dbcomp.pl?class=$class\" TARGET=\"_blank\">$class</A></TD> <TD ALIGN=\"RIGHT\">$count_1</TD> <TD ALIGN=\"RIGHT\">$count_2</TD> <TD ALIGN=\"RIGHT\">$diff</TD> </TR>\n";

  }
  close (COMP);
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

  open (WORMPEP, ">$www/$WS_name/wormpep.shtml") || die "Failed to create wormpep.shtml\n\n";
  print WORMPEP "<P>\n";
  print WORMPEP "<SPACER TYPE=\"horizontal\" SIZE=\"50\">\n";
  print WORMPEP "<TABLE WIDTH=\"40%\" CELLSPACING=\"0\" CELLWIDTH=\"0\" BORDER=\"1\"><TR BGCOLOR=\"yellow\" ><TH>Release</TH> <TH>Sequences</TH> <TH>Letters</TH> <TH>Alt. Splices</TH></TR>\n";


  my ($wp_seq,$wp_let);
  open (WP_1, "</wormsrv2/WORMPEP/wormpep$wormpep_version/wormpep_current.log") || die "Failed to open wormpep.log\n";
  while (<WP_1>) {
    if (/==\> (\S+) sequences totalling (\S+) letters/) {
      ($wp_seq,$wp_let) = ($1,$2);
    }
  }
  close (WP_1);

  
  # get details from last wormpep log file
  my $wp_alt;
  my @possible_logs = `ls -t /wormsrv2/logs/make_wormpep.$wormpep_version*`;
  my $wormpeplog  = "$possible_logs[0]";
  open (WP_2, "<$wormpeplog") || die "Failed to open wormpep.log\n";
  while (<WP_2>) {
    if (/\((\d+)\) alternate splice_forms/) {
      ($wp_alt) = ($1);
    }
  }
  close (WP_2);
  print WORMPEP "<TR><TD>wormpep$wormpep_version</TD> <TD>$wp_seq</TD> <TD>$wp_let</TD> <TD>$wp_alt</TD> </TR>\n";
  print WORMPEP "</TABLE>\n";
  print WORMPEP "</P><PRE>\n";
  

  my (@changed, @lost, @new, @reappeared);
  open (WP_3, "</wormsrv2/WORMPEP/wormpep$wormpep_version/wormpep.diff$wormpep_version") || die "Failed to open wormpep.history\n";
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
    open(EST_HTML,">$www/$WS_name/Checks/$file.shtml") || die "Couldn't create $file.shtml\n";
    
    print EST_HTML &SangerWeb::header();
    
    print EST_HTML "<TABLE WIDTH=\"100%\" CELLSPACING=\"0\" CELLPADDING=\"0\"><TR VALIGN=\"top\" BGCOLOR=\"darkblue\" ><TD WIDTH=\"100%\"><BR><H2 align=\"center\"><FONT COLOR=\"white\">";
    print EST_HTML "EST analysis: \"$files{$file}\"";
    print EST_HTML "</FONT></H2></TD></TR></TABLE>\n";

    print EST_HTML "<P><PRE>";
    open (EST_TXT, "</wormsrv2/autoace/CHECKS/$file.txt") || die "Couldn't open EST_total.txt\n";

    my $line = `wc -l /wormsrv2/autoace/CHECKS/$file.txt`;
    # take line count and add to array
    my ($new) = ($line =~ /(\d+)/);
    push(@line_counts,$new);

    while (<EST_TXT>) {
     print EST_HTML "$_";
    }
    print EST_HTML "</PRE></P>";
    print EST_HTML &SangerWeb::footer();
    close(EST_TXT);
    close(EST_HTML);
  }

  # need to add 'n/a' into the line_counts array to avoid problems later on
  splice(@line_counts,2,0,"n/a");


  print LOG "updating EST_analysis.shtml from $www/$WS_previous/Checks to $www/$WS_name/Checks, adding new info\n"; 

  open(EST_OLD, "<$www/$WS_previous_name/Checks/EST_analysis.html") || die "Couldn't read old version of EST_analysis.html\n";
  open(EST_NEW, ">$www/$WS_name/Checks/EST_analysis.html") || die "Couldn't create new version of EST_analysis.html\n";

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
  my @gff_files = ("clone_ends", "clone_path", "exon","intron_confirmed_CDS", "clone_acc", "genes", "repeats", "intron", "rna");
  foreach my $chrom (@chrom) {
    foreach my $file (@gff_files){
      system("sort -u $gff/CHROMOSOME_$chrom.$file.gff | gff_sort > $www/$WS_name/GFF/CHROMOSOME_$chrom.$file.gff")
	&& die "Couldn't copy CHROMOSOME_$chrom.$file.gff\n";
    }
  }

  print LOG "Copying across GFF files from /wormsrv2/autoace/CHECKS/\n";
  system ("cp /wormsrv2/autoace/CHECKS/*.gff $www/$WS_name/GFF/") && die "Could not copy GFF files from autoace/CHECKS $!\n";
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
   open(OLDFILE, "<$www/$WS_previous_name/Checks/GFF_introns_confirmed_CDS_$lab.html") || die "Couldn't read old $lab GFF intron file\n";
   open(NEWFILE, ">$www/$WS_name/Checks/GFF_introns_confirmed_CDS_$lab.html") || die "Couldn't create new $lab GFF intron file\n";

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
  open(LOGFILE,"</wormsrv2/WORMPEP/wormpep$wormpep_version/wormpep_current.log") || die "Couldn't open wormpep log file\n";
  my $text = <LOGFILE>;
  close(LOGFILE);

  # grab details of wormpep release, number of sequences etc
  my $count = $text;
  my $letters = $text;
  $count =~ s/.*==> (.*) sequences.*/$1/;
  $letters =~ s/.*sequences totalling (.*) letters/$1/;
  # calculate number of splice variants by looking for proteins ending in 'A' in wormpep.table file
  my $alt_spliced = `cut -f 1 /wormsrv2/WORMPEP/wormpep$wormpep_version/wormpep.table$wormpep_version | sed 's/.*[0-9B-Z]\$//' | grep \".\"| wc -l`; 
  $alt_spliced =~ s/\s+//g;

  # create release_paragraph.shtml
  open (PARAGRAPH, ">$www_root/wormpep/release_paragraph.shtml") || die "Can't open the file: $www/release_paragraph.shtml\n\n";
  print PARAGRAPH "The current Wormpep database, Wormpep$wormpep_version (released $release_date), contains $letters residues in $count protein sequences (of which $alt_spliced have splice variants). Wormpep$wormpep_version is based on the <A href=\"ftp://ftp.sanger.ac.uk/pub/wormbase/WS$WS_current\">current WS$WS_current release</A> of the <I>C. elegans</I> AceDB database.\n";
  close (PARAGRAPH);

  # update the 'current_release.shtml' file
  system("rm -f $www_root/wormpep/current_release.shtml") && die "Couldn't remove current_release.shtml\n";
  open(RELEASE,">$www_root/wormpep/current_release.shtml") || die "Coudn't write to current_release.shtml\n";
  print RELEASE "The current release is Wormpep$wormpep_version\n";
  close(RELEASE);

  # update the history of wormpep releases, i.e. wormpep_release.txt

  open (TABLE, "<$www_root/wormpep/wormpep_release.txt") || die "Can't read from the file: $www_root/wormpep/wormpep_release.txt\n\n";
  open (TABLE2, ">$www_root/wormpep/tmp_table.txt") || die "Can't create the file: $www_root/wormpep/tmp_table.txt\n\n";
  print TABLE2 "$wormpep_version\t$release_date2\t$count";
  while(<TABLE>){
    print TABLE2 $_;
  }
  close(TABLE); 
  close(TABLE2);
  system("mv $www_root/wormpep/tmp_table.txt $www_root/wormpep/wormpep_release.txt") && die "Couldn't rename file\n";


  # write a new table for the wormpep_download.shtml page    
  print LOG "Updating wormpep_download.shtml page\n";
  my $rows = $wormpep_version + 5;

  open (LIST, ">$www_root/wormpep/releases.shtml") || die "Can't open the file: $www_root/wormpep/releases.shtml\n\n";
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

    if (($wp_rel < 8) || ($wp_rel == $wormpep_version)) {
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


###########################################################################################################
# Makes the pages of loci designations for each letter of the alphabet
###########################################################################################################

sub create_loci_designations{
  print LOG "\nRunning list_loci_desidnations\n";
  system ("/wormsrv2/scripts/list_loci_designations") && die "Cannot execute list_loci_desinations $!\n";

} 
