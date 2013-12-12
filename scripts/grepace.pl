#!/software/bin/perl -w

# A utility to grep for patetrns that should (-pattern) or should not (-v) occur in a block of ace
# -pattern and -v can be used at the same time

use strict;                              
use Getopt::Long;

my ($test, $pattern, $v, $file);

GetOptions (
	    "pattern:s" => \$pattern,
	    "v:s"       => \$v,
            "file:s"    => \$file,
            )
                    or die("invalid commandline option\n");
            ;



# Change input separator to paragraph mode, but store old mode in $oldlinesep
my $oldlinesep = $/;
$/ = "";

open (FILE, "<$file") or die("cant open $file : $!\n");
while (my $record = <FILE>) {
  while ($record =~ /^\s*\/\//) {$record =~ s/^\s*\/\/.*?\n//} # strip out any comments at the start of the record
  if (defined $pattern && defined $v) {
    if ($record =~ /$pattern/ && $record !~ /$v/) {print $record}
  } else {
    if ((defined $pattern && $record =~ /$pattern/) || (defined $v && $record !~ /$v/) ) {print $record}
  }
}
close FILE;
  
# reset input line separator
$/= $oldlinesep;
