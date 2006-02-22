#!/usr/local/bin/perl5.8.0 -w 
#
#   remap_gff_between_releases.pl                 
# 
# by Gary Williams                         
#
# This remaps the clone positions and genomic
# sequence of two relreases
#
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2006-02-22 17:04:42 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
#use Ace;
#use Sequence_extract;
#use Coords_converter;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($release1, $release2, $version, $gff, $output);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store"      => \$store,
	    "gff=s"      => \$gff,
	    "output=s"   => \$output,
	    "release1=i"  => \$release1,
	    "release2=i"  => \$release2,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


if (! defined $release1 || ! defined $release2) {
  die "Specify the release numbers to use\n";
}
if (! defined $gff || ! defined $output) {
  die "Specify the input and output files\n";
}

##########################
# read in the mapping data
##########################

my @mapping_data = &read_mapping_data($release1, $release2);


##########################
# MAIN BODY OF SCRIPT
##########################


open (OUT, "> $output") || die "Can't open $output";
open (GFF, "< $gff") || die "Can't open GFF file $gff\n";

while (my $line = <GFF>) {
  chomp $line;
  if ($line =~ /^\s*$/) {next;}
  if ($line =~ /^#/) {next;}
  my @f = split /\t/, $line;

      my ($chromosome, $start, $end, $sense) = ($f[0], $f[3], $f[4], $f[6]);

      ($f[3], $f[4], $f[6]) = remap($chromosome, $start, $end, $sense, @mapping_data);

      $line = join "\t", @f;
      print OUT $line,"\n";
}

close (GFF);
close (OUT);

# Close log files and exit
$log->write_to("Finished.\n");

$log->mail();
exit(0);






##############################################################
#
# Subroutines
#
##############################################################


# read_mapping_data
# Usage: @mapping_data = &read_mapping_data($release1, $release2);
# Function: reads the data used to remap across releases
# Args: the first and last wormbase release numbers to use

# the mismatch_start value is the start of the mismatch, it is the first position which doesn't match
# the mismatch_end value is the base past the end of the mismatch region, the first base which matches again
# ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped) 

#
# Input files look like:
# Chromosome: I
# 4765780 4765794 14      4765780 4765794 14      0


sub read_mapping_data {

  # array (one for each release) of hashes (keyed by chromosome number) of list (one for each difference) of list (one for each field)
  # access like: $fields_hashref = @{ $mapping_data[$release]{$chrom}[$next_difference] }
  my @mapping_data;

  foreach my $release (($release1+1) .. $release2) {
    my %chroms;
    my $infile = "/nfs/disk100/wormpub/CHROMOSOME_DIFFERENCES/sequence_differences.WS$release";
    open (IN, "< $infile") || die "Can't open $infile\n";
    my $chrom;
    while (my $line = <IN>) {
      chomp $line;
      if ($line =~ /Chromosome:\s+(\S+)/) {
	$chrom = $1;
      } else {
	my @fields = split /\t/, $line;
	#print "fields=@fields\n";

	push @{ $mapping_data[$release]{$chrom} }, [@fields];

	# debug
	#my $a = $mapping_data[$release]{$chrom} ;
	#foreach my $b (@$a) {
	#  print "just pushed array @$b\n";
	#}
      }
    }
    close(IN);
  }

  return @mapping_data;
}

##########################################

# remap
# Usage: ($f[3], $f[4]) = remap($chromosome, $start, $end, $sense, @mapping_data);
# Function: does the remapping of a pair of values
# Args: $chromosome, $start, $end, $sense, @mapping_data

sub remap {
  my ($chromosome, $start, $end, $sense, @mapping_data) = @_;
  
  if ($chromosome =~ /CHROMOSOME_(\S+)/) {$chromosome = $1;}

  foreach my $release ($release1 .. $release2) {

    if (exists $mapping_data[$release]{$chromosome}) {
      foreach  my $fields (@{$mapping_data[$release]{$chromosome}}) {
	#print "$release $chromosome fields= @$fields \n";

# the mismatch_start value is the start of the mismatch, it is the first position which doesn't match
# the mismatch_end value is the base past the end of the mismatch region, the first base which matches again
# ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped) 
	my ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped) = @$fields;

	if ($flipped) {	   # is the feature inside a flipped region?  

	  if ($start >= $mismatch_start1 && $end < $mismatch_end1) { 

	    if ($sense eq '+') {$sense = '-';} else {$sense = '+';} # flip the sense
	    $start = $mismatch_start1 + $mismatch_end1 - $start; # flip the start and end positions
	    $end = $mismatch_start1 + $mismatch_end1 - $end;
	    if ($start > $end) {
	      my $tmp = $start;
	      $start = $end;
	      $end = $tmp;
	    }
	  }

	} else {
	  # if the start or end are beyond the start of the change region, apply any shift

	  if ($start >= $mismatch_end1) { # if past the end of the change region, shift it
	    $start += $len2 - $len1;
	  } elsif ($start >= $mismatch_start1 && $start - $mismatch_start1 > $len2 ) { # if was in the change region and now out, set it to the end
	    $start = $mismatch_start1 + $len2;
	  }

	  if ($end >= $mismatch_end1) { # if past the end of the change region, shift it
	    $end += $len2 - $len1;
	  } elsif ($end >= $mismatch_start1 && $end - $mismatch_start1 > $len2) { # if was in the change region and now out, set it to the end
	    $end = $mismatch_start1 + $len2;
	  }


	}

      }
    } else {
      #print "no change: doesn't exist: $release $chromosome\n";
    }
  }


  return ($start, $end, $sense);
}

##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################

# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - remap_gff_between_releases.pl

=head1 USAGE

=over 4

=item remap_gff_between_releases.pl [options]

=back

This script reads in a GFF file and the numbers of two releases of
wormbase and maps the chromosomal locations of the first release to
the second release.

script_template.pl MANDATORY arguments:

=over 4

=item -release1 The first (earlier) database to convert from e.g. 140

=back

=item -release2 The second (later) database to convert to e.g. 155

=back

=item -gff the name of the GFF file to read in and convert

=back

=item -outfile the name of the converted GFF file to write out

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Gary Williams

=back

=cut
