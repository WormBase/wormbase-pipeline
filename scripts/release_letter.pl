#!/usr/local/bin/perl5.6.0 -w                    # perl5.6.0 and -w flag
#
# release_letter.pl                            
# 
# by Keith Bradnam                              
#
# Last updated by: $Author: ar2 $                      # These lines will get filled in by cvs and helps us
# Last updated on: $Date: 2002-07-05 16:12:20 $         
#
#Generates a release letter at the end of build.
#Three subroutines are called during the build - 
#release_wormpep by make_wormpep
#release_composition by dump_chromosomes.pl
#release_databases by dbcomp
#These write to a file in autoace/RELEASE_LETTER and are incorperated in to the letter at the end. 
#This allows for overwriting during the build when errors are fixed and scripts rerun




#test the subs#############################
#&release_composition;
#&release_wormpep(19462,20639,1177);#$number_cds $number_total $number_alternate
#&release_databases;
###########################################

use strict;                                      # always use
use lib "/wormsrv2/scripts/";                    # Try to use Wormbase.pm where it useful to do so
use Wormbase;

# Try to keep different parts of code cleanly separated using comments...

##############
# variables  #                                                                   
##############

my $ver = &get_wormbase_version;
my $old_ver = $ver -1;

# Most checking scripts should produce a log file that is a) emailed to us all 
# and b) copied to /wormsrv2/logs

my $maintainer  = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $date        = `date`;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
our $log        = "/wormsrv2/logs/finish_build.$rundate";

#make the release letter
my $release_letter = "/wormsrv2/autoace/RELEASE_LETTERS/letter.WS$ver";
open (RL,">$release_letter");

print RL "New Release of acedb WS$ver, Wormpep$ver and Wormrna$ver $date\n======================================================================\n\n";
print RL "This directory includes:\n
i)\tdatabase.WS$ver.*.tar.gz     -   compressed data for new release\n
ii)\tmodels.wrm.WS$ver           -   the latest database schema (also in above database files)\n
iii)\tCHROMOSOMES/subdir         -   contains 3 files (DNA, GFF & AGP per chromosome)\n
iv)\tWS$ver-WS$old_ver.dbcomp    -   log file reporting difference from last release\n
v)\twormpep$ver.tar.gz           -   full Wormpep distribution corresponding to WS$ver\n
vi)\twormrna$ver.tar.gz          -   latest WormRNA release containing non-coding RNA's in the genome\n
vii)\tconfirmed_genes.WS$ver.gz  -   DNA sequences of all genes confirmed by EST &/or cDNA\n\n\n";

print RL "Release notes on the web:\n-------------------------\n";
print RL "http://www.sanger.ac.uk/Projects/C_elegans/WORMBASE\n\n\n\n";

my $release_dir = ("/wormsrv2/autoace/RELEASE_LETTERS");
my @release_files = ("$release_dir/dbases","$release_dir/composition","$release_dir/wormpep");

#include all the pre-generated reports
my $file = shift(@release_files);
while (defined($file))
       {
	 open (READIN, "<$file") || die "cant open $file\n";
	 while(<READIN>) {
	     print RL "$_";
	   }
	 close READIN;
	 print RL "\n\n\n";
	 $file = shift(@release_files);
       }


# Get the GeneModel corrections
$ver = 81; $old_ver = 80;
my %cam_introns;
$cam_introns{$ver} = `grep CHROMO /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$ver/GFF/CHROMOSOME_*.check_intron_cam.gff | wc -l`;
$cam_introns{$old_ver} = `grep CHROMO /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$old_ver/GFF/CHROMOSOME_*.check_intron_cam.gff | wc -l`;
$cam_introns{change} = $cam_introns {$ver} - $cam_introns {$old_ver};

my %stl_introns;
$stl_introns{$ver} = `grep CHROMO /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$ver/GFF/CHROMOSOME_*.check_intron_stl.gff | wc -l`;
$stl_introns{$old_ver} = `grep CHROMO /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$old_ver/GFF/CHROMOSOME_*.check_intron_stl.gff | wc -l`;
$stl_introns{change} = $stl_introns {$ver} - $stl_introns {$old_ver};

print RL "GeneModel correction progress WS$old_ver -\> WS$ver\n-----------------------------------------\n
Confirmed introns not is a CDS gene model;\n\n\t\t+---------+--------+\n\t\t| Introns | Change |\n\t\t+---------+--------+\n";
printf RL ("Cambridge\t|  %5d  |  %4d  |\n", $cam_introns{$ver},$cam_introns{change});
printf RL ("St Louis \t|  %5d  |  %4d  |\n", $stl_introns{$ver},$stl_introns{change});
print RL "\t\t+---------+--------+\n\n\n";


#Members of known repeat families that overlap predicited exons
my %cam_repeats;
$cam_repeats{$ver} = `grep match /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$ver/Checks/CHROMOSOME_*.repeat_in_exon_cam | wc -l`;
$cam_repeats{$old_ver} = `grep match /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$old_ver/Checks/CHROMOSOME_*.repeat_in_exon_cam | wc -l`;
$cam_repeats{change} = $cam_repeats{$ver} - $cam_repeats{$old_ver};

my %stl_repeats;
$stl_repeats{$ver} = `grep match /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$ver/Checks/CHROMOSOME_*.repeat_in_exon_stl | wc -l`;
$stl_repeats{$old_ver} = `grep match /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/WS$old_ver/Checks/CHROMOSOME_*.repeat_in_exon_stl | wc -l`;
$stl_repeats{change} = $stl_repeats{$ver} - $stl_repeats{$old_ver};

print RL "Members of known repeat families that overlap predicted exons;\n\n\t\t+---------+--------+\n\t\t| Introns | Change |\n\t\t+---------+--------+\n";
printf RL ("Cambridge\t|  %5d  |  %4d  |\n", $cam_repeats{$ver},$cam_repeats{change});
printf RL ("St Louis \t|  %5d  |  %4d  |\n", $stl_repeats{$ver},$stl_repeats{change});
print RL "\t\t+---------+--------+\n\n\n";



#Synchronisation with GenBank / EMBL
#grep ERROR /wormsrv2/autoace/yellow_brick_road/CHROMOSOME_X.agp_seq.log
my @chromosomes = ("I","II","III","IV","V","X");
my $csome = shift @chromosomes;
print RL "Synchronisation with GenBank / EMBL:\n------------------------------------\n\n";
my $check = 0;
while ($csome)
  {
    my $errors = `grep ERROR /wormsrv2/autoace/yellow_brick_road/CHROMOSOME_$csome.agp_seq.log`;
    while( $errors =~ m/for\s(\p{IsUpper}\w+)/g ) {
       print RL "CHROMOSOME_$csome\tsequence $1\n";
       $check = 1;
      }
    $csome = shift @chromosomes;
  }
if ($check == 0){
  print RL "No synchronisation issues\n\n";
}

print RL "----------------------------------------------------------------\n\n\n\n";
print RL "New Data:\n---------\n\n\n";
print RL "New Fixes:\n----------\n\n\n";
print RL "Known Problems:\n--------------\n\n\n";
print RL "Other Changes:\n--------------\n\n";
print RL "Proposed Changes / Forthcoming Data:\n------------------------------------\n\n\n";


print RL "Quick installation guide for UNIX/Linux systems\n-----------------------------------------------\n\n";
print RL "1. Create a new directory to contain your copy of WormBase,\n\te.g. /users/yourname/wormbase\n\n";
print RL "2. Unpack and untar all of the database.*.tar.gz files into\n\tthis directory. You will need approximately 2-3 Gb of disk space.\n\n";
print RL "3. Obtain and install a suitable acedb binary for your system\n\t(available from www.acedb.org).\n\n";
print RL "4. Use the acedb 'xace' program to open your database, e.g.\n\ttype 'xace /users/yourname/wormbase' at the command prompt.\n\n";
print RL "5. See the acedb website for more information about acedb and\n\tusing xace.\n\n";


print RL "____________  END _____________\n";

print "DONT FORGET TO FILL IN THE LAST FEW FIELDS IN THE LETTER\n found at /wormsrv2/autoace/RELEASE_LETTERS/letter.WS$ver\n";

my $name = "Release Letter for WS$ver";
&mail_maintainer($name,$maintainer,$release_letter);


##################################################################################
# Databases used in build                                                        #
##################################################################################
sub release_databases
  {
    #GET THE DATABASES USED IN THIS BUILD

    #get the non_cambridge DB info from file
    open (DBS,"/wormsrv2/autoace/logs/Primary_databases_used_in_build");
    my ($stlace,$citace, $brigdb, $cshace);
    my %dates;
    while(<DBS>)
      {
	chomp;
	my @info = split(/ : /,$_);
	$dates{$info[0]} = $info[1];
      }
    foreach my $key( keys %dates)
      {
	my $oldstyle = $dates{$key};
	my $newstyle = "20".substr($oldstyle,0,2)."-".substr($oldstyle,2,2)."-".substr($oldstyle,4);
	$dates{$key} = $newstyle;
      }
    
    #get the Cambridge dates directly from block file 1
    my @date = &find_file_last_modified("/wormsrv2/camace/database/block1.wrm");
    $dates{camace} = $date[0];
    @date = &find_file_last_modified("/wormsrv2/camace/database/block1.wrm");
    $dates{genace} = $date[0];
    
    foreach my $key( keys %dates)
      {
	print "$key $dates{$key}\n";
      }

    #PARSE THE RELEASE LETTER FOR LAST BUILD INFO
    my $old_ver = &get_wormbase_version -1;
    my $ver = $old_ver+1;
    my %old_dates;
    my $located = 0;
    my $count =0;  #determines how many lines to read = no databases
    open (OLD_LETTER,"/wormsrv2/autoace/RELEASE_LETTERS/letter.WS$old_ver");
    while(<OLD_LETTER>)
      {
	if( ($located == 1) && ($count <= 6))
	  {
	    chomp;
	    my @info = split(/ : | - /,$_);
	    #this will put some crap in the hash from the 1st two line of the section but it no matter
	    $old_dates{$info[0]} = $info[1];
	    $count++;
	  }
	elsif ($_ =~ m/Primary/){
	  $located = 1;
	}
      }  
    my $dbaseFile = "/wormsrv2/autoace/RELEASE_LETTERS/dbases";
    open (WRITE, ">$dbaseFile");
    
    print WRITE "Primary databases used in build WS$ver\n------------------------------------\n";
    foreach my $key(sort keys %dates)
      {
	print WRITE "$key : $dates{$key}";
	if ("$dates{$key}" gt "$old_dates{$key}") {
	  print WRITE " - updated\n";
	}
	elsif ("$dates{$key}" lt "$old_dates{$key}"){
	  print WRITE "you're using a older version of $key than for WS$old_ver ! ! \n";
	}
	else{
	  print "\n";
	}
      }
    close WRITE;

    my $name = "Database update report";
    my $maintainer = "All";
    &mail_maintainer($name,$maintainer,$dbaseFile);
 }

##################################################################################
# Returns the date yyyy-mm-dd and time hh:mm:ss file was last modified           #
##################################################################################
sub find_file_last_modified
  { 
    my $filename = shift;
    open (FILE,"<$filename") || die "cant open file $filename\n";
    my @fileinfo = stat FILE;
    my @date = localtime($fileinfo[9]);
    close(FILE);
    my $year = sprintf("%d-%02d-%02d\n",$date[5]+1900,$date[4]+1,$date[3]);
    my $time = "$date[2]:$date[1]:$date[0]";
    chomp $year;

    my @last_modified = ($year, $time);
    return @last_modified;
  }




##################################################################################
# DNA Sequence composition                                                       #
##################################################################################
sub release_composition
  {
    #get the old info from current_DB
    my $ver = &get_wormbase_version;
    my $old_ver = $ver -1;
    my %old_data;
    my $old_letter ="/wormsrv2/WS$old_ver/CHROMOSOMES/composition.all";
    open (OLD, "<$old_letter") || die "cant open data file - $old_letter";
    while (<OLD>) {
      chomp;
	if ($_ =~ m/(\d+)\s+total$/){
	  #my $tot = $1;$tot =~ s/,//g; # get rid of commas
	  $old_data{Total} = $1;
	}
	elsif ($_ =~ m/^\s+([\w-]{1})\s+(\d+)/){
	  print "adding to old_data $1  $2\n";
	  $old_data{"$1"} = $2;
	}
    }
    print "keys are ",(keys %old_data),"\n";
    close (OLD);
    

    print "\n\n$ver data\n";
    #now get the new stuff to compare
    my $new_letter ="/wormsrv2/WS$ver/CHROMOSOMES/composition.all"; #should be autoace ?
    my %new_data;
    open (NEW, "<$new_letter") || die "cant open data file - $new_letter";
    while (<NEW>) {
      chomp;
      if ($_ =~ m/(\d+)\s+total$/){
	$new_data{Total} = $1;
      }
      elsif ($_ =~ m/^\s+([\w-]{1})\s+(\d+)/){
	print "adding to new_data $1  $2\n";
	$new_data{"$1"} = $2;
      }1
    }
    close NEW;

    #now check the differences
    my %change_data;
    my $compositionFile = "/wormsrv2/autoace/RELEASE_LETTERS/composition";
    open (COMP_ANALYSIS, ">$compositionFile") || die "cant open $compositionFile";
    print COMP_ANALYSIS "Genome sequence composition:\n----------------------------\n\n";
    print COMP_ANALYSIS "       \tWS$ver       \tWS$old_ver      \tchange\n";
    print COMP_ANALYSIS "----------------------------------------------\n";
    foreach my $key(keys %old_data) {
	$change_data{$key} = $new_data{$key} - $old_data{$key};
      }
   
    my @order = ("a","c","g","t","n","-","Total");
    foreach(@order){
      if ("$_" eq "Total"){
	print COMP_ANALYSIS "\n";
      }
      printf COMP_ANALYSIS ("%-5s\t%-8d\t%-8d\t%+4d\n", $_, $new_data{$_}, $old_data{$_}, $change_data{$_});
    }

    if ($change_data{"-"} > 0){
      print COMP_ANALYSIS "Number of gaps has increased - please investigate ! \n";
    }
    
    if ($change_data{"total_length"} < 0) {
      print COMP_ANALYSIS "Total number of bases has decreased - please investigate ! \n";
    }
    close COMP_ANALYSIS;

    my $name = "Sequence composition report";
    my $maintainer = "ar2\@sanger.ac.uk";#"All";
    &mail_maintainer($name,$maintainer,$compositionFile);
  }

##################################################################################
#  Wormpep                                                                       #
##################################################################################
sub release_wormpep       #($number_cds $number_total $number_alternate )
  {  
    my ($number_cds, $number_total, $number_alternate) = @_;
    my $ver = &get_wormbase_version;
    my $old_ver = $ver -1;
    
    #extract data from new wormpep files
    my $lost = `more /wormsrv2/WORMPEP/wormpep$ver/wormpep.diff$ver | grep 'lost' | wc -l`;
    my $new = `more /wormsrv2/WORMPEP/wormpep$ver/wormpep.diff$ver | grep 'new' | wc -l`;
    my $changed = `more /wormsrv2/WORMPEP/wormpep$ver/wormpep.diff$ver | grep 'changed' | wc -l`;
    my $appeared = `more /wormsrv2/WORMPEP/wormpep$ver/wormpep.diff$ver | grep 'appear' | wc -l`;
    my $entries = `cat /wormsrv2/WORMPEP/wormpep$ver/wormpep.diff$ver | wc -l`;
    my $net = $new + $appeared - $lost;
    my $codingDNA;



    #get no of coding bases from log file
    open (THIS_LOG,"/wormsrv2/WORMPEP/wormpep$ver/wormpep_current.log");
    while(<THIS_LOG>){
	if( $_ =~ m/(\d+(,\d+)*)\s+letter/){
	 $codingDNA = $1;}
      }

    #write new letter
    my $wormpepFile = "/wormsrv2/autoace/RELEASE_LETTERS/wormpep";
    open (LETTER, ">$wormpepFile") || die "cant open $wormpepFile\n";


    print LETTER "\n\nWormpep data set:\n----------------------------\n";
    print LETTER"\nThere are $number_cds CDS in autoace, $number_total when counting $number_alternate alternate splice forms.\n
The $number_total sequences contain $codingDNA base pairs in total.\n\n";
   
    print LETTER "Modified entries      $changed";
    print LETTER "Deleted entries       $lost";
    print LETTER "New entries           $new";
    print LETTER "Reappeared entries    $appeared\n";
    printf LETTER "Net change  %+d",$net;

    #get the number of CDS's in the previous build
    open (OLD_LOG,"/wormsrv2/WORMPEP/wormpep$old_ver/wormpep_current.log");
    my $oldCDS;
    while(<OLD_LOG>){
      if( $_ =~ m/(\d+,?\d+)\s+seque/){
	$oldCDS = $1;
	$oldCDS =~ s/,//g;
      }
    }
    close OLD_LOG;

    #check
    my $mail;
    if ($lost + $new + $changed + $appeared != $entries) {
	print LETTER "cat of wormpep.diff$ver does not add up to the changes (from $0)";
      }
    if ($oldCDS + $net != $number_total){
      print LETTER"The differnce between the total CDS's of this ($number_total) and the last build ($oldCDS) does not equal the net change $net\nPlease investigate! ! \n";
    }
       
    close LETTER; 

    my $name = "Wormpep release stats";
    my $maintainer = "All";
    &mail_maintainer($name,$maintainer,$wormpepFile);
  }
#
#                   All these should add up to cat /wormsrv2/WORMPEP/wormpepxx/wormpepxx.diff | wc -l
#                    Calculate net changes (new + reappeared - deleted)
#                    Check whether net change corresponds to change in total CDSs from last week
#


# Always good to cleanly exit from your script
exit(0);
