#!/usr/local/bin/perl5.8.0 -w
#
# update_caltech.pl
#
# by Chao-Kung Chen [030113]
#
# Automatically update Geneace with Erich's functional annotation update
#
# Last updated on: $Date: 2003-12-08 12:58:28 $
# Last updated by: $Author: ck1 $

use strict;                    
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Cwd 'chdir';
use Getopt::Long;

###################################################
# check user is wormpub otherwise script won't run
###################################################

my $user = `whoami`; chomp $user;
if ($user ne "wormpub"){
  print "\nYou have to be wormpub to run this script!\n\n";
  exit(0);
}

# touch logfile for run details
$0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");


###################################################
# variables and command-line options with aliases 
###################################################

my ($debug, $help, $update, $merge, $recipients);
$recipients ="bastiani\@its.caltech.edu, ck1\@sanger.ac.uk, emsch\@its.caltech.edu, kimberly\@minerva.caltech.edu, krb\@sanger.ac.uk";
my $tace = &tace;   # tace executable path

GetOptions ("d|debug=s"  => \$debug,
	    "u|update"   => \$update,
	    "m|merge"    => \$merge,
            "h|help"     => \$help,
           );

if (!$merge && !$update){
  print "\nYou need to specify [-m | -merge] or [-u | -update] to proceed!\n\n";
  exit(0);
}

# touch logfile for run details
$0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");


my $script = $0;
my $rundate = `date +%y%m%d`; chomp $rundate;
my $log = "/wormsrv2/logs/update_caltech.$rundate.$$";
my $nlog = "/wormsrv2/logs/new_gene_name_merge.$rundate.$$";

if ($help){
  system("perldoc update_caltech.pl");
  exit (0);
}

if ($debug){
  print "DEBUG : $debug\n";
  $recipients = "$debug\@sanger.ac.uk";
}
	
my $caltech="/wormsrv1/geneace/ERICHS_DATA/FTP";
my $updatedir="/wormsrv1/geneace/ERICHS_DATA";
my (@dates, @ftp_date, @last_date);
my (@info, @exceptions, %exceptions, @other_names, %other_names, @other_main, %other_main, %locus_cds);

open(LOG, ">$log") || die "Can't write to log file $log\n";
open (NLOG, ">$nlog") || die "Can't write to log file $nlog\n";

system("chmod 777 $log $nlog");
my $found = 0;   # counter for found gene name merge

######################
#   start working
######################

# see POD for -update and -merge options
my $count = &download if $update;  

&process_main_other_name;                                 # get list of main names, their other name(s) and corresponding seq. names

$found = &check_gene_name_merge if $merge;                # email a notice of new gene name merge

&main_other_name_assignment if ($update && $count != 0);  # update geneace with latest func. annots; assign annots to main name, if not already
                                                          # email Caltech what needs to be changed 

# Mail to people, but only to person running script in debug mode
mail_maintainer("Functional annotation update feedback", $recipients, $log) if $update;
mail_maintainer("New gene name merge notice", $recipients, $nlog) if $found == 1;

# when no new gene name merge is found
if ($merge && $found == 0){
    $recipients ="ck1\@sanger.ac.uk"; # notify ck1 if no update
    $recipients ="$debug\@sanger.ac.uk" if $debug; 
    mail_maintainer("New gene name merge notice", $recipients, $nlog);
}

exit(0);


########################
#      subroutines
########################

sub download {  

  # check for new file to upate on FTP site


  print LOG "This file is generated automatically. If you have spotted any bug, please contact ck1\@sanger.ac.uk\n";
  print LOG "--------------------------------------------------------------------------------------------------\n\n";
  
  chdir $caltech;
  system("echo '\$ caltech' | ftp -i caltech.wormbase.org") && print LOG "Failed to download file\n";
  
  @ftp_date=dataset($caltech);    # [0] is the date of the file , [1] is the filename
  @last_date=dataset($updatedir); # [0] is the date of the file , [1] is the filename
   
  if ($last_date[0] != $ftp_date[0]){
    print LOG "UPDATE file $ftp_date[1] avilable on FTP\n\n";
    system("mv $ftp_date[1] ../");
    return $count = 1;
  }  
  if ($last_date[0] == $ftp_date[0]){
    print LOG "No new update on FTP site\n";
    $recipients ="ck1\@sanger.ac.uk, krb\@sanger.ac.uk"; # notify ck1 & krb if no update
    $recipients ="$debug\@sanger.ac.uk" if $debug;       # notify ck1 & krb if no update
    mail_maintainer("Functional annotation update feedback", $recipients, $log); 
    exit(0) if !$merge;
    return $count = 0;
  }
}

sub process_main_other_name {
  
  ########################################
  # process locus / other_name information
  ########################################
  
  my $ga_dir = "/wormsrv1/geneace";
  #my $ga_dir = "/nfs/disk100/wormpub/DATABASES/BACKUPS/geneace_backup.031107";

  my $locus_has_other_name_to_cds=<<EOF;
Table-maker -p "/wormsrv1/geneace/wquery/locus_has_other_name_to_cds.def" quit
EOF
  my $locus_to_CDS=<<EOF;
Table-maker -p "/wormsrv1/geneace/wquery/locus_to_CDS.def" quit
EOF
  
  @other_main=CGC_loci_and_other_name($locus_has_other_name_to_cds, $locus_to_CDS, $ga_dir);
  
  %other_main=%{$other_main[0]};
  %locus_cds=%{$other_main[1]};

  push(@info, "\n-----------------------------------------------------------------------------\n"); 
  push(@info, "FYI - CURRENT LIST ($rundate) OF LOCI THAT ARE BOTH MAIN NAME AND OTHER NAME:\n");
  push(@info, "-----------------------------------------------------------------------------\n");

  foreach (sort keys %other_main){
    my $locus = $other_main{$_}->[0];
    my @cds   = $other_main{$_}->[1];
    
    if (exists $locus_cds{$_}){
      push(@info, "$_ (@{$locus_cds{$_}}) is different from $_ which is an other_name of $locus (@cds)\n");
      push(@info, "Annotation needs to be assigned to $_ (@{$locus_cds{$_}})\n\n");
      push(@exceptions, $_);
    }
    else {
      push(@other_names, $_);
    }
  }

  # hash set for quick look up
  foreach (@exceptions){$exceptions{$_}++}
  foreach (@other_names){$other_names{$_}++}
}

sub check_gene_name_merge {  

  # write current list of exception loci (as CGC main name and non-CGC other name as well)
  my $exception_loci_list = "/wormsrv1/geneace/ERICHS_DATA/exception_loci_list.$rundate";
  open(EXCP, ">$exception_loci_list") || print NLOG "\nCannot write to $exception_loci_list!\n";
  foreach (sort keys %exceptions){
    print EXCP "$_\n";
  }
  close EXCP;

  # write current list of main loci to a file
  my $main_loci_list = "/wormsrv1/geneace/ERICHS_DATA/main_loci_list.$rundate";
  open(MAIN, ">$main_loci_list") || print NLOG "\nCannot write to $main_loci_list!\n";
  foreach (sort keys %locus_cds){
    print MAIN "$_ -> @{$locus_cds{$_}}\n";
    
  }
  close MAIN;
  
  # grep the two main loci files from two different dates
  my @file_date;
  my @loci_files = glob"/wormsrv1/geneace/ERICHS_DATA/main_loci_list*";
  foreach (@loci_files){
    if ($_ =~ /.+main_loci_list\.(\d+)/){
      push(@file_date, $1);
    }
  }

  @file_date = sort {$a <=> $b} @file_date;
  # diff the two main loci files to see what is changed
  my @name_diff = `diff $loci_files[0] $loci_files[1]`;
  
  # process diff result for the names found, check also what is its main name now
  my @merger;

  foreach (@name_diff){
    chomp;
    #print $_, "\n";
    if ($_ =~ /^< (.+) -> (.+)/){  # older file
      my $name = $1;
      my @seqs = $2;
      print $name, "\n";
      # ensures that the name (gone) is an other-name of a main name and is by itself not a main name
      if (exists $other_main{$name} && !exists $locus_cds{$name}){
        $found = 1;
        push(@merger, "$name became an other_name of the main name $other_main{$name}->[0]($other_main{$name}->[1])\n");
      }
    }
  }

  # grep the two exceiption loci files from two different dates
  @loci_files = glob"/wormsrv1/geneace/ERICHS_DATA/exception_loci_list*";
  foreach (@loci_files){
    if ($_ =~ /.+main_loci_list\.(\d+)/){
      push(@file_date, $1);
    }
  }

  @file_date = sort {$a <=> $b} @file_date;
  # diff the two exception loci files to see what is changed
  @name_diff = `diff $loci_files[0] $loci_files[1]`;

  foreach (@name_diff){
    chomp;   
    print $_, "\n";
    if ($_ =~ /^> (.+)/){
      my $name = $1;
      print $name, "\n";
      # ensures that the exception name is an other-name of a main name and is by itself a CGC main name
      $found = 1; 
      push(@merger, "$name became an other_name of the main name $other_main{$name}->[0] ($other_main{$name}->[1]), this is different from the CGC main name $name ($locus_cds{$name}->[0])\n");
      push(@merger, "$name became an other_name of the main name $other_main{$name}->[0] ($other_main{$name}->[1]), this is different from the CGC main name $name (not connected to sequence yet)\n") if $locus_cds{$name}->[0] eq "NA";
    }
  }

  if ($found == 1){
    print NLOG "This notice is generated by script, please contact ck1\@sanger.ac.uk if you have spotted any error.\n\n";
    print NLOG "-----------------------------------\n";
    print NLOG "New gene name merge(s) as of $rundate\n";
    print NLOG "-----------------------------------\n\n";
    print NLOG @merger; 
    print NLOG "\nPlease update other related annotations accordingly. Thanks.\n" if $found == 1;
  }
  else {
    $found = 0;
    print NLOG "This notice is generated by script, please contact ck1\@sanger.ac.uk if you have spotted any error.\n\n";
    print NLOG "No new gene name merge.\n";
  }
	
  #remove previous files: exception_ / main_loci_list.date  (there should always be one for each of these in this folder)
  system("rm -f /wormsrv1/geneace/ERICHS_DATA/*loci_list.$file_date[0]");

  return $found;
}

sub main_other_name_assignment {

  open (IN, "$updatedir/$ftp_date[1]");
  my $modify ="$updatedir/$ftp_date[1].modified";
  
  open(OUT, ">$modify");

  my @assignment;
  push (@assignment,"-----------------------------------------------------------------\n");
  push (@assignment,"B: Detailed FUNCTIONAL ANNOTATION ASSIGNMENT BASED ON YOUR UPDATE\n");
  push (@assignment,"-----------------------------------------------------------------\n");
  
  while(<IN>){
    if ($_ =~ /^Locus : \"(.+)\"/){
      
      my $flag = $1;
      ######################################################### 
      # check if locus is an other_name and also a locus object
      #########################################################
      
      if ($exceptions{$flag}){
	print $flag, "\n";
	my $locus = $other_main{$flag}->[0]; # main
	my @cds   = $other_main{$flag}->[1] if $other_main{$flag}->[1];
	
	# locus attached to cds
	if ($other_main{$flag}->[1]){   
	  push (@assignment, "$flag (@{$locus_cds{$1}}) is different from $flag which is an other_name of $locus (@cds)\n");
	  push (@assignment, "Functional annotation is assigned to $flag (@{$locus_cds{$flag}} in geneace)\n\n");
	  print OUT "Locus : \"$flag\"\n";
	}
	
	# locus not attached to cds
	else {                             
	  push (@assignment, "$flag (@{$locus_cds{$flag}}) is different from $flag which is an other_name of $locus (no CDS connection yet)\n");
	  push (@assignment, "Functional annotation is assigned to $flag (@{$locus_cds{$flag}} in geneace)\n\n");
	  print OUT "\nLocus : \"$flag\"\n";
	}
      }
      
      ##################################################################### 
      # check if locus should become an other_name of a CGC locus object
      #####################################################################
      
      elsif ($other_names{$flag}){
	my $locus = $other_main{$flag}->[0]; # main name
	my @cds   = $other_main{$flag}->[1]; # seq name if available
	
	# other name has main name attached to cds
	if ($other_main{$flag}->[1] ne "NA"){
	  push (@assignment, "$flag is an other_name of $locus (@cds)\n");
	  push (@assignment, "Functional annotation is assigned to $locus\n\n");
	  print OUT "-R Locus : \"$flag\" \"$locus\"\n";
	  print OUT "\nLocus : \"$locus\"\n";
	  print OUT "Other_name \"$flag\"\n";
	}
	
	# other name has main name not attached to cds
	else {
	  push (@assignment, "$flag is an other_name of $locus (no CDS connection yet)\n");
	  push (@assignment, "Functional annotation is assigned to $locus\n\n");
	  print OUT "-R Locus : \"$flag\" \"$locus\"\n";
	  print OUT "\nLocus : \"$locus\"\n";
	  print OUT "Other_name \"$flag\"\n";	
	}
      }
      else {
	print OUT "$_";
      }
    }
    else {
      print OUT "$_";
    }
  }

  my $diff_out = `diff $updatedir/$ftp_date[1] $modify`;
  
  if (!$diff_out){ 
    print LOG "------------------------------------------\n\n";
    print LOG "EVERYTHING WENT OK\n\n";
    print LOG "THANKS.\n\n";
   # print LOG @info;
  }
  else {
    print LOG "\n\n";
    print LOG "------------------------------------\n";
    print LOG "A: CHANGES MADE FOR YOUR UPDATE FILE\n";
    print LOG "------------------------------------\n";  
    print LOG $diff_out, "\n";
    print LOG "------------------------------------\n\n";
    print LOG "Please update accordingly\n";
    print LOG "THANKS.\n\n\n";
    print LOG @assignment;
    print LOG "\n";
    print LOG @info;
  }

  my $command=<<END;
find CDS * where concise_description OR detailed_description OR provisional_description
show -a -T -f /wormsrv1/geneace/ERICHS_DATA/CDS_TS_dump.ace
edit -D Concise_description
edit -D Detailed_description
edit -D Provisional_description

find Transcript * where concise_description OR detailed_description OR provisional_description
show -a -T -f /wormsrv1/geneace/ERICHS_DATA/TRANSCRIPT_TS_dump.ace
edit -D Concise_description
edit -D Detailed_description
edit -D Provisional_description

find Pseudogene * where concise_description OR detailed_description OR provisional_description
show -a -T -f /wormsrv1/geneace/ERICHS_DATA/PSEUDOGENE_TS_dump.ace
edit -D Concise_description
edit -D Detailed_description
edit -D Provisional_description

find locus * where concise_description OR detailed_description OR provisional_description
show -a -T -f /wormsrv1/geneace/ERICHS_DATA/loci_TS_dump.ace
edit -D Concise_description
edit -D Detailed_description
edit -D Provisional_description

pparse $modify

save
quit
END


  my $geneace_dir="/wormsrv1/geneace";
  
  if (!$debug){
    open (Load_GA,"| $tace -tsuser \"Functional_annotation\" $geneace_dir >> $log") || die "Failed to upload to Geneace";
    print Load_GA $command;
    close Load_GA;
    
    if($!){
      print "##########################################\n";
      print "\nPhenotype annotation is now updated!\n\n";       
      print "\nIf everthing is OK, REMEMBER to remove\n"; 
      print "loci_TS_dump.ace and seq_TS_dump.ace\n";
      print "in /wormsrv1/geneace/ERICHS_DATA\n\n";
      print "##########################################\n\n"; 
    }
    else{
      print "######################################\n";
      print "Mission not 100% successful, Mr. Bond!\n";
      print "######################################\n\n";
    }
    
    ########################################################
    # move modified file/last update file to ARCHIVE folder
    #######################################################
    
    chdir $updatedir;
    system("mv $last_date[1]* ARCHIVE/");
  }  
}

sub dataset {
  my $dir = shift;
  my ($date, $name);
  opendir(DIR, $dir) || die "Can't read directory";
  my @dir=readdir DIR;
  splice(@dir, 0,2);
  closedir (DIR);
  foreach (@dir){
    if ($_ =~ /^annots-(\d+)(\w{3,3})(\d+)_CDS\.ace$/){
      $name = $_;
      my $mon = $2;
      if ($mon eq "jan"){$mon = "01"}
      if ($mon eq "feb"){$mon = "02"}
      if ($mon eq "mar"){$mon = "03"}
      if ($mon eq "apr"){$mon = "04"}
      if ($mon eq "may"){$mon = "05"}
      if ($mon eq "jun"){$mon = "06"} 
      if ($mon eq "jul"){$mon = "07"}
      if ($mon eq "aug"){$mon = "08"}
      if ($mon eq "sep"){$mon = "09"}
      if ($mon eq "oct"){$mon = "10"}
      if ($mon eq "nov"){$mon = "11"}
      if ($mon eq "dec"){$mon = "12"}
      $date = $1.$mon.$3;
    }
  }
  return $date, $name; 
}

sub CGC_loci_and_other_name {

  my ($def1, $def2, $dir) = @_;
  my ($main, $other_name, %other_main, $locus, @seqs, $cds, %locus_cds);

  open (FH1, "echo '$def1' | tace $dir |") || die "Couldn't access geneace\n";

  # hash, key = other-name, value = locus (main), seqs (or NA, if not available)
  while (<FH1>){
    chomp $_;
    if ($_ =~ /\"(.+)\"\s+\"(.+)\"\s+\"(.+)\"/){
      $main = $1;
      $other_name = $2;
      $other_name =~ s/\\//;
      if ($3){
	$cds = $3;
	push(@{$other_main{$other_name}}, $main, $cds);
      }
      else {
	push(@{$other_main{$other_name}}, $main, "NA");
      }
    }
  }
  close FH1;

  open (FH2, "echo '$def2' | tace $dir |") || die "Couldn't access geneace\n";
  while (<FH2>){
    chomp $_;
    if ($_ =~ /^\"/){
      # a locus maybe linked to zero, or 1 or 2 or 3 seq. classes
      @seqs = split(/\s+/, $_) if $_ =~ /^\"/;
      $locus = $seqs[0]; $locus =~ s/\"//g;
      
      shift @seqs; # only seq.
      
      # hash, key = locus (main), value = (seq. linked to it)
      if (@seqs){
	foreach (@seqs){
	  $cds = $_; $cds =~ s/\"//g;
	  push(@{$locus_cds{$locus}}, $cds);
	}
      }
      else {
	push(@{$locus_cds{$locus}}, "NA");
      }
    }
  }
  close FH2;
  return \%other_main, \%locus_cds;
}

__END__

=head2 NAME - update_caltech.pl  


            This script automatically checks for available functional annotation update at Caltech FTP site.
            It emails to ck1 and krb no matter update is available or not.
            When update is available log messages is also included in the mail notice. 
            Erich at Caltech will also receive mail notice automatically.

            The script also checks for loci that are other_names in Caltech update file. 
            If an other_name is still a separate object, it will be listed in the mail message 
            and need to be looked at by hand.
            If a locus has been merged to a main name, that locus will be changed to its 
            main name in the update ace file. 

=head3 <USAGE> 
 

=head2 Options: [d or debug] [h or help]

 
B<-debug:>   
            Debug mode, follow by user name, eg. -d ck1 or -debug ck1
            Email will be sent to only the person running the script


B<-help:>      Read this POD

B<-merge:>     Check what is new for main, other-name assignments
            (gene name merge) and send email to Caltech curators

B<-update:>    Check Caltech FTP site for updates on functional annotations

=head3 <RUN update_caltech.pl>

            The script is now set to run 
                 on Sunday mornings at 007 with the -u (update) switch
                 from Mon-Fri at 8pm with the -m (merge) switch
            
           

