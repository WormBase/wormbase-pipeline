#!/usr/local/bin/perl5.6.1 -w

# Author: Chao-Kung Chen
# 2001-11-28
# What it does: (1) Backup MySQL BlastP and BlastX data to ~wormpub/MYSQL_DUMPS
#               (2) Delete old blast data   

use strict;


my $date = `date +%y%m%d`;
chomp $date;

my $start=`date +%H:%M:%S`; chomp $start;

my $worm01="worm01_$date";
my $wormprot="wormprot_$date";

# password is not included for security reasons

print "Backup Blast_X homols...\n\n";
system("mysqldump -h ecs1f -u wormadmin  --opt -p worm01 > /nfs/disk100/wormpub/MYSQL_DUMPS/$worm01");

print "Backup Blast_P homols...\n\n";
system("mysqldump -h ecs1f -u wormadmin --opt -p wormprot > /nfs/disk100/wormpub/MYSQL_DUMPS/$wormprot");

my @files =`ls /nfs/disk100/wormpub/MYSQL_DUMPS`;
foreach (@files){
  chomp;
  if( ($_ ne "worm01_$date") || ($_ ne "wormprot_$date") ){ 
   system ("rm -rf $_");
  }
}

my $end=`date +%H:%M:%S`; chomp $end;

print "\nBackup started at $start\n";
print "Backup finished at $end\n";


