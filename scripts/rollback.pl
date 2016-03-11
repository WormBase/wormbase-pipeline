#!/software/bin/perl -w

# Usage: rollback.pl infile outfile
# Function: change a .ace file into a .ace file that removes the data that the first .ace file loaded in

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Carp;

my $infile = shift @ARGV;
my $outfile = shift @ARGV;

open (IN, "< $infile") || die "Can't open $infile\n";
open (OUT, "> $outfile") || die "Can't open $outfile\n";


my $object = 1; # next line is a object start line

while (my $line = <IN>) {
  if ($line =~ /^\s*$/) {
    $object = 1;
    next;
  }
  if ($object) {
    print OUT "\n$line";
    $object = 0;
  } else {
    if ($line =~ /^\/\//) { # // comment
      print OUT "$line";
    } else {
      print OUT "-D $line";
    }
  }
}
close(IN);
close (OUT);



