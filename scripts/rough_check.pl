#!/usr/local/bin/perl
#
# rough_checks.pl
#
#
#
#

$|=1;
use lib "/wormsrv2/scripts/";   
#use strict;
use Wormbase;
use IPC::Open2;
use IO::Handle;
use Getopt::Std;
use vars qw( $opt_n $opt_s $opt_e $opt_p $opt_d $opt_h);


 ##############################
 # Script variables (run)     #
 ##############################

my $maintainer = "All";
my $rundate = `date +%y%m%d`;   chomp $rundate;
my $runtime = `date +%H:%M:%S`; chomp $runtime;

$tace = &tace;

 ##############################
 # command-line options       #
 ##############################

$opt_n = "";       # write acefiles from autoace.config file
$opt_e = "";       # write non-elegans EST dataset
$opt_d = "";       # debug mode (maintainer = dl1)
$opt_h = "";       # Help/Usage page

getopts ('d');

# cheat for quiet debugging
# set dl1 as only mail recipient if running debug mode
$maintainer = "dl1\@sanger.ac.uk" if ($opt_d);

$dbpath = "/wormsrv2/autoace";

&predicted_genes;

exit(0);

 ########################################
 # Open logfile                         #
 ########################################

my $logfile = "/wormsrv2/logs/rough_check.$rundate.$$";
system ("/bin/touch $logfile");
open (LOG,">>$logfile") or die ("Could not create logfile\n");
LOG->autoflush();
open (STDOUT,">>$logfile");
STDOUT->autoflush();
open (STDERR,">>$logfile"); 
STDERR->autoflush();

print LOG "# rough_check\n\n";     
print LOG "# run details    : $rundate $runtime\n";
print LOG "\n";
print LOG "WormBase version : ${WS_version}\n";
print LOG "\n";
print LOG "======================================================================\n";
print LOG " -n : Write .ace files\n"                            if ($opt_n);
print LOG "    : from database $opt_s\n"                        if ($opt_s);
print LOG " -e : Write non-elegans nematode EST data sets\n"    if ($opt_e);
print LOG " -p : Write peptide data sets\n"                     if ($opt_p);
print LOG "======================================================================\n";
print LOG "\n";
print LOG "Starting make_acefiles .. \n\n";


&predicted_genes;


my $endtime = `date +%H:%M:%S`; chomp $endtime;
print LOG "Ended make_acefiles @ $endtime\n";
close STDERR;
close STDOUT;
close LOG;



 ##############################
 # mail $maintainer report    #
 ##############################

&mail_maintainer("WormBase Report: rough_check ",$maintainer,$logfile);

# hasta luego
exit (0);

#################################################################################
### Subroutines                                                               ###
#################################################################################

sub predicted_genes {

format LIST =
@<<<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<< [@<] @>> @>> @<<<<<<<<< @<<<<<<<<<<<<<<<<<<<<<< [@<<<<<<<<] 
$f[0],             $f[1],         $f[3],$f[2],$noest,$f[7],     $f[4], $confirmed
.

$~ = "LIST";
$ENV{'ACEDB'} = $dbpath;
$command=<<EOF;
Table-maker -p "/wormsrv2/autoace/wquery/roughcheck.def"
quit
EOF
    
open (TACE, "echo '$command' | $tace | ");
while (<TACE>) {
    next if (/\/\//);
    s/acedb\> //g;
    chomp;
    next if  ($_ eq "");
    s/\"//g;

    undef ($confirmed);
    $noest = 0;

    @f = split (/\t/);

#    print "'$f[0]' '$f[1]' '$f[2]' '$f[3]' '$f[4]' '$f[5]' '$f[6]' '$f[7]'\n"; 

    if ($f[6] eq "") {
	$noest = 0;
    } 
    else {
	$noest = $f[6];
    }

    if ($f[5] eq "Confirmed_by") {
	$confirmed = "Confirmed";
    }
    elsif ($noest >= 1) {
	$confirmed = "Supported";
    }
    
    write;


}


}
