#!/usr/bin/env perl
#
# cleanse_gff.pl
# 
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-08-15 09:04:49 $

# 'Cleanes' GFF by:
# (a) removing lines with source = '.'
# (b) ensures all lines have 9 columns ("." is added to col 9 if not the case)

use strict;
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

my ($debug, $test, $store, $wormbase, $species);

my ($gff3,$infile,$outfile, $lines_modified, $lines_removed);

GetOptions (
  "debug=s"    => \$debug,
  "test"       => \$test,
  "store:s"    => \$store,
  "species:s"  => \$species,
  "infile:s"   => \$infile,
  "outfile:s"  => \$outfile,
  "gff3"       => \$gff3,
    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
			     );
}


# establish log file.
my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

while (<$gff_in_fh>) {
  if (/^\#/) {
    print $gff_out_fh $_;
    next;
  }
  chomp;
    
  my @l = split(/\t+/, $_);
  if (scalar(@l) == 8) {
    # no group field, so explicit add an empty one
    push @l, ".";
    $lines_modified++;
  } elsif (scalar(@l) != 9) {
    $log->log_and_die(sprintf("Found line in GFF with bad number of columns:\n$_\n"));
  }
  
  if ($l[1] eq '.') {
    if ($l[2] =~ 'Clone_left_end' or $l[2] =~ 'Clone_right_end') {
      # do this just so that everything has a non '.' source
      $l[1] = 'Clone';
    } else {
      $lines_removed++;
      next;
    }
  }

  print $gff_out_fh join("\t", @l), "\n";
}
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");

$log->write_to("Finished cleansing : $lines_modified lines modified, $lines_removed lines removed\n");

$log->mail();
exit(0);

