#!/usr/local/bin/perl5.6.1 -w

# EMBLDump.pl :  makes EMBL dumps from camace.
#  by  ag3
#  Last updated on: $Date: 2002-12-09 13:50:03 $
#  Last updated by: $Author: ck1 $

# [dl1 000618] : Altered script to accomodate change of lab 'CB' -> 'HX'

$0 =~ s/^\/.+\///;
system ("touch /wormsrv2/logs/history/$0.`date +%y%m%d`");

use strict;
use IO::Handle;
use IPC::Open2;
use POSIX qw(:signal_h :errno_h :sys_wait_h);
use lib "/wormsrv2/scripts";
use Wormbase;

# Avoid filling process table with zombies

$SIG{CHLD} = \&REAPER;
sub REAPER {
  my $pid;
  $pid=waitpid(-1,&WNOHANG);
  $SIG{CHLD}=\&REAPER;
}


my $exec        = &giface;
my $dbdir       = "/wormsrv2/camace";
my $giface      = "$exec $dbdir";
my $outfilename = "/nfs/disk100/wormpub/tmp/EMBLdump.$$";

# Dump the EMBL file from camace

open (READ,$giface) or die ("Could not open $giface\n"); 
open (WRITE,$giface) or die ("Could not open $giface\n"); 

my $query=<<END;
Query Find Genome_Sequence From_Laboratory = HX AND Finished AND DNA
gif EMBL $outfilename
END
print "Query = $query\n";
print "Exec = $giface\n";
print WRITE $query;
close WRITE;
while (<READ>) {
  if ($_ =~ /\/\//) {next};
  if ($_ =~ /acedb/) {next};
}		    
close READ;

print "Outfile is $outfilename";
exit;








