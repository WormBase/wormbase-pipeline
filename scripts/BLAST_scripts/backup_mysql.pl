#!/usr/local/bin/perl5.6.1 -w

# Author: Chao-Kung Chen
# 2001-11-28
# What it does: (1) Backup MySQL BlastP and BlastX data to ~wormpub/MYSQL_DUMPS
#               (2) Delete old blast data   

use strict;
use Getopt::Long;

my ($dbpass, $dump_worm_dna, $dump_worm_pep, $dump_brigpep);

GetOptions ( 
	    'dbpass=s'     => \$dbpass,
	    'worm_dna'     => \$dump_worm_dna,
	    'worm_pep'     => \$dump_worm_pep,
	    'worm_brigpep' => \$dump_brigpep,
	   );

die "I need the password, please - try something like - \nbackup_mysql.pl -dbpass  xxxxxxxx\n" unless $dbpass;

unless ( defined $dump_worm_dna or defined $dump_worm_pep or defined $dump_brigpep ) {
  $dump_worm_dna = 1; 
  $dump_worm_pep = 1;
  $dump_brigpep = 1;
}

my $date = `date +%y%m%d`;
chomp $date;

my $start=`date +%H:%M:%S`; chomp $start;
print "Backup started at $start\n\n";

my $worm_dna="worm_dna_$date";
my $worm_pep="wormpep_$date";
my $worm_brigpep="worm_brigpep_$date";

# password is not included for security reasons

&dump_it("worm_dna",$worm_dna) if $dump_worm_dna;
&dump_it("worm_pep",$worm_pep) if $dump_worm_pep;
&dump_it("worm_brigpep",$worm_brigpep) if $dump_brigpep;

my $end=`date +%H:%M:%S`; chomp $end;

print "Backup finished at $end\n";


sub dump_it
  {
    my $database = shift;
    my $dump_file = shift;

    print "Backup $database to $dump_file...\n\n";
    system("mysqldump -h ecs1f -u wormadmin --opt -p$dbpass $database > /nfs/disk100/wormpub/MYSQL_DUMPS/$dump_file");
}
