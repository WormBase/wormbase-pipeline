#!/usr/local/bin/perl5.6.1 -w

# update_KO_alleles.pl

# by Chao-Kung Chen [030113]

# Last updated on: $Date: 2003-02-06 18:12:14 $
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
#$0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

####################################
# variables and command line options
####################################

my $debug;
my $rundate = `date +%y%m%d`; chomp $rundate;
my $allele_dir = "/wormsrv1/geneace/ALLELE_DATA";
my $ftp_dir = "$allele_dir/KO_FTP";
my $log;
my $allele_carers = "ck1\@sanger.ac.uk, krb\@sanger.ac.uk";
my $script = $0; $script =~ s/.pl//;

# debug mode
GetOptions ("d|debug=s"  => \$debug);

if ($debug){
  print "\nDEBUG = $debug\n\n";
  $allele_carers = "$debug\@sanger.ac.uk";
}

##########
# log file
##########

open(LOG, ">/wormsrv2/logs/KO_allele_update.$rundate.$$") || die "Can't write to log file $log!";;
print LOG "This KO_allele update notice has been sent automatically\n\n";
print LOG "$script.pl started to run on $rundate\n";
print LOG "=======================================================================================\n";

####################################
# checking updates on KO FTP site
# if update file(s) available: 
#        (1) download
#        (2) rename downloaded files
####################################

chdir $ftp_dir;
#system("ftp -i -v aceserver.biotech.ubc.ca") && die "print LOG $!\n";

my $updates=data_set($ftp_dir, "ftp");

print "SEQ: $$updates[0]\n"; # current FTP update for KO_sequenced_alleles.rundate.ace
print "PCR: $$updates[1]\n"; # current FTP update for KO_PCR_primers.rundate.ace

my ($last, $current, @dates);

open (DATE, "$allele_dir/last_KO_update") || die "Can't read in file!";

#=start
while(<DATE>){
  if ($_ =~ /^sequenced_KO_alleles.ace.(\d+)/){
    $last = $1;
    if ($$updates[0] =~ /^sequenced_KO_alleles.ace.(\d+)/){$current = $1}
    push(@dates, $last, $current);
    my @order = sort {$a <=> $b} @dates; print "SEQ: @order#########\n";
    if ($current eq $last || $current eq $order[0]){
      print LOG "\nLatest sequenced_KO_allele data on FTP from: $current  <--->  Last Sanger update from $last\n";
      print LOG "----> No new data to update\n";
    }
    else{
      print LOG "\nLatest sequenced_KO_allele data on FTP from: $current  <--->  Last Sanger update from $last\n";	
      print LOG "New data on KO FTP site available..started to download...\n";

      ###################################
      # download Sequenced alleles update
      ###################################

      chdir $allele_dir;
      system("ftp -i -v aceserver.biotech.ubc.ca");
      $$updates[0] =~ /^sequenced_KO_alleles.ace.(\d+)/;
      system("mv $$updates[0] KO_sequenced_alleles.$1.ace#"); 	
   }
  }
  @dates =();
  if ($_ =~ /^PCR_primers.ace.(\d+)/){
    $last = $1;
    if ($$updates[1] =~ /^PCR_primers.ace.(\d+)/){$current = $1}
    push(@dates, $last, $current);
    my @order = sort {$a <=> $b} @dates; print "PCR: @order#########\n";
    if ($current eq $last || $current eq $order[0]){
      print LOG "\nLatest PCR_primers data on FTP from: $current  <--->  Last Sanger update from: $last\n";
      print LOG "No new data to update\n";
    }
    else{
      print LOG "\nLatest PCR_primers date on FTP from: $current  <--->  Last Sanger update from $last";
      print LOG "\nNew data on KO FTP site available..started to download...\n";

      #############################
      # download PCR_primers update
      #############################

      chdir $allele_dir;
      system("ftp -i -v aceserver.biotech.ubc.ca");
      $$updates[1] =~ /^PCR_primers.ace.(\d+)/;
      system("mv $$updates[1] KO_PCR_primers.$1.ace#");
   }
  }	
}		 

###########################################
# make del.ace files for downloaded file(s)
###########################################

open(SEQ, "$allele_dir/KO_sequenced_alleles.$rundate.ace") || die "Can't read in file!";
open(PCR, "$allele_dir/KO_PCR_primers.$rundate.ace") || die "Can't read in file!";
open(OSEQ, ">$allele_dir/KO_sequenced_alleles.$rundate.del.ace") || die "Can't read in file!";
open(OPCR, ">$allele_dir/KO_PCR_primers.$rundate.del.ace") || die "Can't read in file!";

while(<SEQ>){
  chomp;
  if ($_ =~ /^Allele :.+/){print OSEQ "\n$_\n"}
  elsif ($_ eq ""){last}
  else {print OSEQ "-D $_\n"}
}
while(<PCR>){
  chomp;
  if ($_ =~ /^PCR_product :.+/ || $_ =~ /^Oligo :.+/){print OPCR "\n$_\n"}
  elsif ($_ eq ""){last}
  else {print OPCR "-D $_\n"}
}

######################################
# upload ace files to Geneace
#
# make backup allele class w/ TS
# first load del.ace from last update
# then load current update file
# delete last del.ace file(s)
######################################

my @dels=($allele_dir, "delACE");

print "SEQd: $dels[0]\n";  # last del.ace for KO_sequenced_alleles.rundate.del.ace
print "PCRd: $dels[1]\n";  # last del.ace for KO_PCR_primers.rundate.del.ace


mail_maintainer($0, $allele_carers, $log);
#=end
#=cut

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
  if ($file eq "ftp"){
    foreach (@dir){
      print $_, "###\n";
      if ($_ =~ /^sequenced_KO_alleles.ace.(\d+)/){push(@Supdate, $1)}
      if ($_ =~ /^PCR_primers.ace.(\d+)/){push(@Pupdate, $1)}
    }
    if (scalar @Supdate>1){
      @S_last = sort {$a <=> $b} @Supdate;
      #print $S_last[-2], "##\n";
      push(@update, "sequenced_KO_alleles.ace.$S_last[-1]");
    }
    if (scalar @Supdate==1){
      push(@update, "sequenced_KO_alleles.ace.@Supdate");
    }
    if (scalar @Pupdate >1){
      @P_last = sort {$a <=> $b} @Pupdate;
      print @P_last, "##\n";
      push(@update, "PCR_primers.ace.$P_last[-1]");
    }
    if (scalar @Pupdate==1){
      push(@update, "PCR_primers.ace.@Pupdate");
    }

    return \@update;
  }
  if ($file eq "delACE"){
    foreach (@dir){
      if ($_ =~ /^sequenced_KO_alleles.\d+.del.ace/){push(@Sdel, $_)}
      if ($_ =~ /^^PCR_primers.\d+.del.ace/){push(@Pdel, $_)}
    }
    @Pd_last = sort {$a <=> $b} @Pdel;
    @Sd_last = sort {$a <=> $b} @Sdel;

    if (scalar @Pd_last == 1 || scalar @Sd_last == 1){
      if (scalar @Pd_last == 1){
	print LOG "\nNo last del.ace for KO_PCR_primers.rundate.del.ace available.\n";
	print LOG "Script stops.....NO PCR_primers UPDATE for now.\n";
      }
      if (scalar @Sd_last == 1){
	 print LOG "\nNo last del.ace for KO_sequenced_alleles.rundate.del.ace available\n";
	 print LOG "Script stops.....NO Sequenced alleles UPDATE for now.\n";
      }
      exit(0);
    }	
    else {return \$Pd_last[-2],  $Sd_last[-2]}
  }
}










