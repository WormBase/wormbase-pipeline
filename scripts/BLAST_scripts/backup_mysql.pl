#!/usr/local/bin/perl5.6.1 -w

# Author: Chao-Kung Chen
# 2001-11-28
# What it does: (1) Backup MySQL BlastP and BlastX data to ~wormpub/MYSQL_DUMPS
#               (2) Delete old blast data   

use strict;
use Getopt::Long;

my ($dbpass, $dump_worm01, $dump_wormprot, $dump_brigprot);

GetOptions ( 
	    'dbpass=s' => \$dbpass,
	    'worm01'   => \$dump_worm01,
	    'wormprot' => \$dump_wormprot,
	    'worm_brigprot' => \$dump_brigprot,
	   );

die "I need the password, please - try something like - \nbackup_mysql.pl -pass  xxxxxxxx\n" unless $dbpass;

unless ( defined $dump_worm01 or defined $dump_wormprot or defined $dump_brigprot ) {
  $dump_worm01 = 1; 
  $dump_wormprot = 1;
  $dump_brigprot = 1;
}

my $date = `date +%y%m%d`;
chomp $date;

my $start=`date +%H:%M:%S`; chomp $start;
print "Backup started at $start\n\n";

my $worm01="worm01_$date";
my $wormprot="wormprot_$date";
my $worm_brigprot="worm_brigprot_$date";

# password is not included for security reasons

&dump_it("worm01",$worm01) if $dump_worm01;
&dump_it("wormprot",$wormprot) if $dump_wormprot;
&dump_it("worm_brigprot",$worm_brigprot) if $dump_brigprot;

my $end=`date +%H:%M:%S`; chomp $end;

print "Backup finished at $end\n";


sub dump_it
  {
    my $database = shift;
    my $dump_file = shift;

    print "Backup $database to $dump_file...\n\n";
    system("mysqldump -h ecs1f -u wormadmin --opt -p$dbpass $database > /nfs/disk100/wormpub/MYSQL_DUMPS/$dump_file");
}
