#!/usr/local/bin/perl5.6.1 -w
#
# Dump_Subseq.pl
#
# by Dan Lawson
#
# Last updated on: $Date: 2002-12-09 14:54:10 $
# Last updated by: $Author: krb $
#
# Ace dumps a subsequence from a given target database
# And retrieves the coordinates from parent object

#####################################################################################################

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use IPC::Open2;
use Getopt::Std;

our ($opt_d, $opt_s, $opt_n, $keyset);
getopts ('n:d:s:');

if ($opt_n) {
    $keyset = $opt_n;
    $opt_d = "/wormsrv2/autoace";
    $opt_s = "dumped.subsequence";
}

if ((!$opt_d)||(!$opt_s)) {
 &PrintHelp;
 exit 0;
}

my $dbpath=$opt_d;
$opt_d =~ /camace/  && do {$dbpath = "/wormsrv2/camace";};
$opt_d =~ /autoace/ && do {$dbpath = "/wormsrv2/autoace";};
$opt_d =~ /stlace/  && do {$dbpath = "/wormsrv2/stlace";};
my $cwd = `/bin/pwd`;
chomp $cwd;
print "DBPATH: $dbpath\n";
my $tace= &tace;
my $database = "$tace $dbpath";

my $outfile = "$cwd/" ."$opt_s".".ace";

print "** Output file : $outfile **\n";
if ($opt_s!~/(\w+)\.\w+/) {
  print "Sorry - $opt_s is not a subsequence\n";
  exit 0;
}

if ($keyset) {
    open (KEYSET, $keyset) || die "Could not open the keyset file\n";
    while (<KEYSET>) {
	if (/^Sequence : \"(\S+)\"/) {
	    $opt_s = $1;
	    print "Retrieve Subsequence $opt_s\n";
	    &retrieve_subsequence;
	}
    }
    close (KEYSET);
} else {
    &retrieve_subsequence;
}
exit 0;

################################################################################################

sub retrieve_subsequence {

my @object;
print "Retrieving subsequence $opt_s ..\n";
open2(\*READ,\*WRITE,$database) or die ("Could not open the connection with $opt_d\n");
my $QUERY=<<END;
find Sequence $opt_s
show -a
END
print WRITE $QUERY;
close WRITE;

my $SOURCE;
while (<READ>) {  
  /\/\// && next;
  /acedb/ && next;
  /Source\s+/ && do {$SOURCE=$_;chomp $SOURCE; print "Source: $SOURCE\n"};
  push (@object, $_);
}
close READ;
$SOURCE =~ s/Source\s+//;
$SOURCE =~ s/\"//mg;

print "Retrieving parent $SOURCE ..\n";
open2(\*READ,\*WRITE,$database) or die ("Could not open the connection with $opt_d\n");
my $QUERY2=<<END;
find Sequence $SOURCE
show -a
END
print WRITE $QUERY2;
close WRITE;

open (OUTFILE,">>$outfile");
print OUTFILE "\nSequence : \"$SOURCE\" \n";
while (<READ>) {
  /\/\// && next;
  /acedb/ && next;
  /Subsequence\s+\"$opt_s\"/ && print OUTFILE $_;
}
close READ;


print OUTFILE "\n";
foreach (@object) {
    chomp;
    if ($_ eq "") {next;}
    print OUTFILE "$_\n";
}
close OUTFILE;

print " .. DONE. Results are in $outfile\n";

}

#---------------------
# Prints documentation
#
sub PrintHelp {
   exec ('perldoc',$0);
}


__END__

=pod

=head1   NAME - Dump_Subseq.pl

=head2 USAGE

Dump_Subseq.pl will extract a Subsequence object (with coordinates from the parent)
from an ACEDB database.

Dump_Subseq.pl mandatory arguments [n excludes d and s]

=over 4

=item -d Database: one of autoace, camace, stlace

=item -s Subsequence_object_name: requested subsequence

=item -n keyset_file_name: dumps a keyset of Subsequences from camace

=back

=cut






