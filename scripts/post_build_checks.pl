#!/usr/local/bin/perl5.6.1 -w
#
# wrapper to call introncheck, estcheck and overlapcheck.pl on gff files, 
# recommended after database rebuilt :o)
#
# by Kerstin Jekosch
# 13/07/01
#
# N.B. Previously called gffcheck
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2003-01-24 10:18:20 $


use Getopt::Std;
use IO::Handle;
$|=1;
use lib "/wormsrv2/scripts/";
use Wormbase;
use strict;

###########
# options #
###########

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

our $log;
&create_log_files;

#########################
# tiny little main loop #
#########################

&runintroncheck  if ($opt_i);

&runestcheck     if ($opt_e);

&runoverlapcheck if ($opt_o);


##############################
# Tidy up                    #
##############################

close LOG;
&mail_maintainer("post_build_checks.pl",$maintainers,$log);

exit(0);




##############################################################
#
# Subroutines
#
##############################################################

sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate     = `date +%y%m%d`; chomp $rundate;
  $log        = "/wormsrv2/logs/$script_name.$rundate.$$";

  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name\n";
  print LOG "started at ",`date`,"\n";
  print LOG "=============================================\n";
  print LOG "\n";

}

##########################################


sub runestcheck {
  my $runtime = `date +%H:%M:%S`; chomp $runtime;
  print LOG "Starting estcheck at $runtime\n";
  system ("/wormsrv2/scripts/estcheck") && die "Can't run estcheck\n";
  $runtime = `date +%H:%M:%S`; chomp $runtime;
  print LOG "Finished running estcheck at $runtime\n\n";
}

sub runintroncheck {
  my $runtime = `date +%H:%M:%S`; chomp $runtime;
  print LOG "Starting introncheck at $runtime\n";
  system ("/wormsrv2/scripts/introncheck") && die "cant run /wormsrv2/scripts/introncheck\n";
  $runtime = `date +%H:%M:%S`; chomp $runtime;
  print LOG "Finished running introncheck at $runtime\n\n";
}

sub runoverlapcheck {
  my $runtime = `date +%H:%M:%S`; chomp $runtime;
  print LOG "Starting overlapcheck.pl at $runtime\n";
  system ("/wormsrv2/scripts/overlapcheck.pl") && die "Can't run overlapcheck.pl\n";
  $runtime = `date +%H:%M:%S`; chomp $runtime;
  print LOG "Finished running overlapcheck.pl at $runtime\n\n";
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
