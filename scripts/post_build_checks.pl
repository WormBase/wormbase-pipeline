#!/usr/local/bin/perl5.6.0 -w
#
# wrapper to call introncheck, estcheck and overlapcheck on gff files, 
# recommended after database rebuilt :o)
#
# by Kerstin Jekosch
# 13/07/01
#
# N.B. Previously called gffcheck


use Getopt::Std;
use IO::Handle;
$|=1;
use lib "/wormsrv2/scripts/";
use Wormbase;
use strict;

###########
# options #
###########
#our($opt_a,$opt_o,$opt_e,$opt_i,$opt_w,$opt_l,$opt_h);

use vars qw / $opt_a $opt_o $opt_e $opt_i $opt_w $opt_l $opt_h /;

$opt_a="";   # does all the following
$opt_o="";   # performs overlapcheck
$opt_e="";   # performs estcheck 
$opt_i="";   # performs introncheck
$opt_w="";   # performs copy2web.pl
$opt_l="";   # performs list_loci_designations
$opt_h="";   # Help/Usage page

getopts ('aeiowlh');
&usage if ($opt_h);
if ($opt_a) { $opt_e = 1; $opt_i = 1; $opt_o = 1; $opt_w = 1; $opt_l = 1;}

#################
# set variables #
#################

my $maintainers = "dl1\@sanger.ac.uk kj2\@sanger.ac.uk krb\@sanger.ac.uk";
my $rundate = `date +%y%m%d`; chomp $rundate;
my $runtime = `date +%H:%M:%S`; chomp $runtime;

# grab version number from cvs
my $cvs_version = &get_cvs_version($0);

my $WS_release = get_wormbase_version; 
my $WS_version = substr ($WS_release,2);

#################
# open LOG file #
#################

my $logfile = "/wormsrv2/logs/post_build_checks.WS${WS_version}.${rundate}.$$";
open (LOG,">$logfile") || die "Cannot open logfile $!\n";
LOG->autoflush();

print LOG "# post_build_checks.pl\n\n";     
print LOG "# run details    : $rundate $runtime\n";
print LOG "# cvs version: $cvs_version\n\n";
print LOG "WormBase version : WS${WS_version}\n";

print LOG "======================================================================\n";
print LOG "  -a : executes all of the following -eio\n" if ($opt_a);
print LOG "  -e : executes estcheck\n"                  if ($opt_e);
print LOG "  -i : executes introncheck\n"               if ($opt_i);
print LOG "  -o : executes overlapcheck\n"              if ($opt_o);
print LOG "  -l : executes list_loci_designations\n"    if ($opt_l);
print LOG "  -w : executes copy2web.pl\n"               if ($opt_w);
print LOG "======================================================================\n";
print LOG "\n";

#########################
# tiny little main loop #
#########################

&runestcheck     if ($opt_e);
&runoverlapcheck if ($opt_o);
&runintroncheck  if ($opt_i);
&run_list_loci_designations if ($opt_l);
&runcopy2web     if ($opt_w);


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
    system ("/wormsrv2/scripts/estcheck \&") && die "Cannot execute estcheck $!\n";
    print LOG "Run estcheck\n";
}

sub runintroncheck {
    system ("/wormsrv2/scripts/introncheck \&") && die "Cannot execute introncheck $!\n";
    print LOG "Run introncheck\n";
}

sub runoverlapcheck {
    system ("/wormsrv2/scripts/overlapcheck") && die "Cannot execute overlapcheck $!\n";
    print LOG "Run overlapcheck\n";
}
 
sub runcopy2web {
    system ("/wormsrv2/scripts/copy2web.pl \&") && die "Cannot execute copy2web.pl $!\n";
    print LOG "Run copy2web.pl\n";
}

sub run_list_loci_designations {
    system ("/wormsrv2/scripts/list_loci_designations \&") && die "Cannot execute list_loci_desinations $!\n";
    print LOG "Run list_loci_designations\n";
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
within genes (overlapcheck).  It can also run the list_loci_designations script to update the website
with the latest gene->sequence assignments

post_build_checks.pl mandatory arguments:

=over 4

=item none, (but it won't do anything)

=back

post_build_checks.pl OPTIONAL arguments:

=over 4

=item -a, executes all of the following -eiowl

=item -e, executes estcheck

=item -i, executes introncheck 

=item -o, execute overlapcheck

=item -w, execute copy2web.pl

=item -l, execute list_loci_designations

=back

=cut
