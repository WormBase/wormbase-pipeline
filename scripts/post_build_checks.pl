#!/usr/local/bin/perl5.6.0 -w
#
# wrapper to call introncheck, estcheck and overlapcheck.pl on gff files, 
# recommended after database rebuilt :o)
#
# by Kerstin Jekosch
# 13/07/01
#
# N.B. Previously called gffcheck
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2002-11-22 10:08:16 $


use Getopt::Std;
use IO::Handle;
$|=1;
use lib "/wormsrv2/scripts/";
use Wormbase;
use strict;

###########
# options #
###########
#our($opt_a,$opt_o,$opt_e,$opt_i,$opt_l,$opt_h);

use vars qw / $opt_a $opt_o $opt_e $opt_i $opt_h /;

$opt_a="";   # does all the following
$opt_o="";   # performs overlapcheck.pl
$opt_e="";   # performs estcheck 
$opt_i="";   # performs introncheck
$opt_h="";   # Help/Usage page

getopts ('aeioh');
&usage if ($opt_h);
if ($opt_a) { $opt_e = 1; $opt_i = 1; $opt_o = 1;}

#################
# set variables #
#################

my $maintainers = "All";
my $rundate = `date +%y%m%d`; chomp $rundate;
my $runtime = `date +%H:%M:%S`; chomp $runtime;


my $WS_version = &get_wormbase_version;

#################
# open LOG file #
#################

my $logfile = "/wormsrv2/logs/post_build_checks.WS${WS_version}.${rundate}.$$";
open (LOG,">$logfile") || die "Cannot open logfile $!\n";
LOG->autoflush();

print LOG "# post_build_checks.pl\n\n";     
print LOG "# run details    : $rundate $runtime\n";
print LOG "WormBase version : WS${WS_version}\n";

print LOG "======================================================================\n";
print LOG "  -a : executes all of the following -e -i -o\n" if ($opt_a);
print LOG "  -e : executes estcheck\n"                      if ($opt_e);
print LOG "  -i : executes introncheck\n"                   if ($opt_i);
print LOG "  -o : executes overlapcheck.pl\n"                  if ($opt_o);
print LOG "======================================================================\n";
print LOG "\n";

#########################
# tiny little main loop #
#########################

&runestcheck     if ($opt_e);
&runintroncheck  if ($opt_i);
&runoverlapcheck if ($opt_o);

##############################
# mail $maintainer report    #
##############################

close LOG;

open (mailLOG, "|/usr/bin/mailx -s \"Report: post_build_checks.pl\" $maintainers ") || die "Cannot open mailLOG $!\n";
open (readLOG, "<$logfile") || die "Cannot open readLOG $!\n";
while (<readLOG>) {
    print mailLOG $_;
}
close readLOG || die "Cannot close readLOG $!\n";
close mailLOG || die "Cannot close mailLOG $!\n";

exit(0);

###############
# subroutines #
###############

sub runestcheck {
    print LOG "Starting estcheck\n";
    system ("/wormsrv2/scripts/estcheck");
    print LOG "Finished running estcheck\n";
}

sub runintroncheck {
    print LOG "Starting introncheck\n";
    system ("/wormsrv2/scripts/introncheck") or die "cant run /wormsrv2/scripts/introncheck\n";
    print LOG "Finished running introncheck\n";
}

sub runoverlapcheck {
    print LOG "Starting overlapcheck.pl\n";
    system ("/wormsrv2/scripts/overlapcheck.pl");
    print LOG "Finished running overlapcheck.pl\n";
}
   
sub usage {
    system("perldoc /wormsrv2/scripts/post_build_checks.pl") && die "Cannot help you, sorry $!\n";
    exit (0);
}
__END__

=pod

=head2   NAME - post_build_checks.pl

=head1 USAGE

=over 4

=item post_build_checks.pl [-options]

=back

post_build_checks.pl is a wrapper to drive scripts to check the gff files for confirmed introns (introncheck), 
inconsistencies in EST assignments (estcheck), overlapping genes, ESTs matching introns and repeats 
within genes (overlapcheck.pl).  

post_build_checks.pl mandatory arguments:

=over 4

=item none, (but it won't do anything)

=back

post_build_checks.pl OPTIONAL arguments:

=over 4

=item -a, executes all of the following -eiowl

=item -e, executes estcheck

=item -i, executes introncheck 

=item -o, execute overlapcheck.pl

=back

=cut
