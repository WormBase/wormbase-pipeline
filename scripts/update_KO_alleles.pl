#!/usr/local/bin/perl5.6.1 -w

# update_KO_alleles.pl

# by Chao-Kung Chen [030113]

# Last updated on: $Date: 2003-03-03 12:13:39 $
# Last updated by: $Author: ck1 $


# Automatically update KO sequenced alleles

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Cwd 'chdir';
use Getopt::Long;

#######################
# check user is wormpub
#######################

my $user = `whoami`; chomp $user;
if ($user ne "wormpub"){
  print "\nYou have to be wormpub to run this script!\n\n";
  exit(0);
}

# touch logfile for run details
$0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

####################################
# variables and command line options
####################################

my ($debug, $help);
my $rundate = `date +%y%m%d`; chomp $rundate;
my $tace = &tace;
my $allele_dir = "/wormsrv1/geneace/ALLELE_DATA";
my $ftp_dir = "$allele_dir/KO_FTP";
my $allele_carers = "ck1\@sanger.ac.uk, krb\@sanger.ac.uk";

my $script = $0; $script =~ s/.pl//;

# debug mode
GetOptions ("d|debug=s"  => \$debug,
	    "h|help"     => \$help,
	   );

if ($debug){
  print "\nDEBUG = $debug\n\n";
  $allele_carers = "$debug\@sanger.ac.uk";
}

if ($help){
  system("perldoc $0")
}

##########
# log file
##########

my $log = "/wormsrv2/logs/KO_allele_update.$rundate.$$";
open(LOG, ">>$log") || die "Can't write to log file $log!";
print LOG "This KO_allele update notice has been sent automatically\n\n";
print LOG "$script.pl started to run on $rundate\n";
print LOG "=======================================================================================\n";

#######################################
# checking updates on KO FTP site
# if both update file(s) are available: 
#
# sequenced_KO_alleles.ace
# PCR_primers.ace
#
# --> download
#######################################

chdir $ftp_dir;
system("echo '\$ knockout' | ftp -i aceserver.biotech.ubc.ca") && die "print LOG $!\n";

my $updates=data_set($ftp_dir, "ftp");
#print "FTP SEQ: $$updates[0]\n"; # current update for sequenced_KO_alleles.ace.rundate on FTP
#print "FTP PCR: $$updates[1]\n"; # current update for PCR_primers.ace.rundate on FTP

my $last_updates=data_set($allele_dir, "last");
#print "GA SEQ: $$last_updates[0]\n"; # last update for sequenced_KO_alleles.ace.rundate in geneace
#print "GA PCR: $$last_updates[1]\n"; # last update for PCR_primers.ace.rundate in geneace

my ($ftp_seq_date, $ftp_pcr_date);

my $seq_update = 0;
my $pcr_update = 0;

if ($$updates[0] eq $$last_updates[0]){  # seqd allele file
  print LOG "No update for sequenced_KO_alleles.ace\n"; 
  
}
else {
  $$updates[0] =~ /^sequenced_KO_alleles.ace.(\d+)/;
  $ftp_seq_date = $1;
  $$last_updates[0] =~ /^sequenced_KO_alleles.ace.(\d+)/;
  my $last_update = $1;
  push(my @dates, $ftp_seq_date, $last_update);
  my @order = sort {$a <=> $b} @dates; 
  if ($order[0] eq $ftp_seq_date){print LOG "No new sequenced_KO_alleles.ace update\n"}
  if ($order[1] eq $ftp_seq_date){ $seq_update = 1; print LOG "New sequenced_KO_alleles.ace update available\n"}
}

if ($$updates[1] eq $$last_updates[1]){ # PCR file
  print LOG "No update for PCR_primers.ace\n"; 
}
else {
  $$updates[1] =~ /^PCR_primers.ace.(\d+)/;
  $ftp_pcr_date = $1; 
  $$last_updates[1] =~ /^PCR_primers.ace.(\d+)/;
  my $last_update = $1;
  push(my @dates, $ftp_pcr_date, $last_update);
  my @order = sort {$a <=> $b} @dates; 
  if ($order[0] eq $ftp_pcr_date){print LOG "No new PCR_primers.ace update\n"}
  if ($order[1] eq $ftp_pcr_date){$pcr_update = 1; print LOG "New PCR_primers.ace available\n"}
}  

#######################
# download update files
#######################
  
if ($seq_update = 1 && $pcr_update = 1){
  system ("mv $$updates[0]  $$updates[1] ../") && print LOG "Couldn't move $$updates[0] and $$updates[1] to $allele_dir\n";
  chdir $allele_dir;
  system ("mv $$last_updates[0] $$last_updates[1] ARCHIVE") && print LOG "Couldn't move $$updates[0] and $$updates[1] to 
$allele_dir/ARCHIVE\n"; 
}

#####################################################################################
# exit program if not both sequenced_KO_alleles.ace and PCR_primers.ace are available 
#####################################################################################

if ($seq_update == 0 || $pcr_update == 0){
  print LOG "Since not both sequenced_KO_alleles.ace and PCR_primers.ace updates are available,\n"; 
  print LOG "no knockout allele update is done\n";
  exit(0);
}

###########################################
# make del.ace files for downloaded file(s)
###########################################

my $seqfile = "$allele_dir/$$updates[0]";
my $pcrfile = "$allele_dir/$$updates[1]";
my $delseq  = "$allele_dir/sequenced_KO_alleles.del.ace.$ftp_seq_date";
my $delpcr  = "$allele_dir/PCR_primers.del.ace.$ftp_pcr_date";

open(SEQ, $seqfile) || die "Can't read in file!";
open(PCR, $pcrfile) || die "Can't read in file!";
open(OSEQ, ">$delseq") || die "Can't write to in file!";
open(OPCR, ">$delpcr") || die "Can't write to in file!";

while(<SEQ>){
  chomp;
  if ($_ =~ /^Allele : .+/){print OSEQ "\n$_\n"}
  else {
    if ($_ eq ""){}
    else {print OSEQ "-D $_\n"}
  }
}  
while(<PCR>){
  chomp;
  if ($_ =~ /^PCR_product : .+/ || $_ =~ /^Oligo : .+/){print OPCR "\n$_\n"}
  else {
    if ($_ eq ""){}
    else {print OPCR "-D $_\n"}
  }
}

######################################
# upload ace files to Geneace
#
# make backup allele class w/ TS
# first load del.ace from last update
# then load current update file
# delete last del.ace file(s)
######################################

########################
# get last del.ace files
########################

my @dels=data_set($allele_dir, "delACE");

#print LOG "Last Seq'd allele del.ace: $dels[0]\n"; 
#print LOG "Last PCR_primer del.ace: $dels[1]\n";

################
# upload geneace 
################

my $command=<<END;
find allele *
show -a -T -f /wormsrv1/geneace/ALLELE_DATA/allele_TS_dump.ace.$rundate

pparse $dels[0]
pparse $dels[1]
pparse $seqfile
pparse $pcrfile 

save
quit
END

my $ga_dir="/wormsrv1/geneace";

open (Load_GA,"| $tace -tsuser \"KO_consortium_allele_update\" $ga_dir >> $log") || die "Failed to upload to Geneace";
print Load_GA $command;
close Load_GA;

system ("\\rm $dels[0] $dels[1]"); 

#############
# mail notice
#############

mail_maintainer($0, $allele_carers, $log);

#############
# subroutines
#############

sub data_set {
  my ($dir, $file) = @_;
  my (@dir, @update, @Supdate, @Pupdate, @P_last, @S_last, @Sdel, @Pdel, @Pd_last, @Sd_last);	

  opendir(DIR, $dir) || die "Can't read directory";
  @dir=readdir DIR;
  splice(@dir, 0,2);
  closedir (DIR);
  if ($file eq "ftp" || $file eq "last"){
    foreach (@dir){
      if ($_ =~ /^sequenced_KO_alleles.ace.(\d+)/){push(@Supdate, $1)}
      if ($_ =~ /^PCR_primers.ace.(\d+)/){push(@Pupdate, $1)}
    }
    if (scalar @Supdate>1){
      @S_last = sort {$a <=> $b} @Supdate;
      push(@update, "sequenced_KO_alleles.ace.$S_last[-1]"); # ensure to get the latest file of seqd allele in this folder
    }
    if (scalar @Supdate==1){
      push(@update, "sequenced_KO_alleles.ace.@Supdate");
    }
    if (scalar @Pupdate >1){
      @P_last = sort {$a <=> $b} @Pupdate;
      push(@update, "PCR_primers.ace.$P_last[-1]");   # ensure to get the latest file of PCR_primers in this folder
    }
    if (scalar @Pupdate==1){
      push(@update, "PCR_primers.ace.@Pupdate");
    }

    return \@update;
  }
  if ($file eq "delACE"){
    foreach (@dir){
      if ($_ =~ /^sequenced_KO_alleles.del.ace.(\d+)/){push(@Sdel, $1)}
      if ($_ =~ /^^PCR_primers.del.ace.(\d+)/){push(@Pdel, $1)}
    }
    @Pd_last = sort {$a <=> $b} @Pdel;
    @Sd_last = sort {$a <=> $b} @Sdel;
    if (scalar @Pd_last == 1 || scalar @Sd_last == 1){
      if (scalar @Pd_last == 1){
	print  LOG "\nNo last del.ace for PCR_primers.ace available.\n";
	print  LOG "Script stops.....NO PCR_primers UPDATE for now.\n";
      }
      if (scalar @Sd_last == 1){
	 print  LOG "\nNo last del.ace for sequenced_KO_alleles.ace available\n";
	 print  LOG "Script stops.....NO Sequenced alleles UPDATE for now.\n";
      }
      exit(0);
    }	
    else {return "sequenced_KO_alleles.del.ace.$Sd_last[-2]", "PCR_primers.del.ace.$Pd_last[-2]"}
  }
}


__END__

=head2 NAME - update_KO_alleles.pl  


            This script automatically checks for available Knockout consortium alleles update on their FTP site.

            When updates are available, it downloads two files: 
              PCR_primers.ace.new_date                    
              sequenced_KO_alleles.ace.new_date
            to /wormsrv1/geneace/ALLELE_DATA/KO_FTP/
            
            where new_date = latest version of these files on FTP site which is represented by yymmdd format            

            In the /wormsrv1/geneace/ALLELE_DATA/ foleder, there should always be these 4 files (so that the script will work): 
              last update files:
                1. PCR_primers.ace.old_date  
                2. sequenced_KO_alleles.ace.old_date   
                3. PCR_primers.del.ace.old_date   
                4. sequenced_KO_alleles.del.ace.old_date

            where old_date = last version of these files on FTP site which is represented by yymmdd format


            When updates are available, the first two files will be moved to /wormsrv1/geneace/ALLELE_DATA/ARCHIVE
            and replaced by 
              PCR_primers.ace.new_date                    
              sequenced_KO_alleles.ace.new_date
            The last two del.ace files will be deleted and two new del.ace files of the latest update will be created.


=head2 Options: [d or debug] [h or help]

 
B<-debug:>   
            Debug mode, follow by user name, eg. -d ck1 or -debug ck1

            No email will be sent to anyone including the person running the script,
            because there are too MANY emails already. You just need to look for 
            "/wormsrv2/logs/update_caltech.rundate.pID" file. 


B<-help:>      Quess what!


=head3 <RUN update_caltech.pl>

            Cron job is set to run on Mondays 007 hour.









