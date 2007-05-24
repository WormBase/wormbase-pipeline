#!/usr/bin/env perl
#===============================================================================
#
#         FILE: Remap_Sequence_Change.pm
#
#  DESCRIPTION: Functions for remapping locations based on changes in
#               chromosome sequences between database releases.
#
#        FILES:  ---
#         BUGS:  ---
#        NOTES: 
#      $Author: gw3 $
#      COMPANY:
#     $Version:  $
#      CREATED: 2006-02-27
#        $Date: 2007-05-24 09:05:03 $
#===============================================================================
package Remap_Sequence_Change;

use strict;
use warnings;
#use Modules::Coords_converter;
use Coords_converter;
use Wormbase;

##########################################################
# 
# Name:      read_mapping_data
# Usage:     @mapping_data = &read_mapping_data($release1, $release2);
# Function:  reads the data used to remap across release
# Args:      $release1, $release2, the first and last wormbase release
#                  numbers to use e.g. 140, 150 to convert data made using wormbase
#                  release WS140 to the coordinates of release WS150
# Returns:   the mapping data, for use in remap_gff()
#

# the mismatch_start value is the start of the mismatch, it is the first position which doesn't match
# the mismatch_end value is the base past the end of the mismatch region, the first base which matches again
# ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped)
                                                                                                                                                            
#
# Input files look like:
# Chromosome: I
# 4765780 4765794 14      4765780 4765794 14      0
                                                                                                                                                            
                                                                                                                                                            
sub read_mapping_data {

  my ($release1, $release2) = @_;
                                                                                                                                           
  # array (one for each release) of hashes (keyed by chromosome number) of list (one for each difference) of list (one for each field)
  # access like: $fields_hashref = @{ $mapping_data[$release]{$chrom}[$next_difference] }
  my @mapping_data;
                                                                                                                                                            
  foreach my $release (($release1+1) .. $release2) {
    my %chroms;
    my $infile = "/lustre/cbi4/work1/wormpub/CHROMOSOME_DIFFERENCES/sequence_differences.WS$release";
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


##########################################################
# 
# Name:      remap_test
# Usage:     $changed = remap_test($release1, $release2, @mapping_data);
# Function: test to see if there have been indel or reverse of
#           orientation changes between the given versions.
#           If so then we will need to run the remapping programs.
# Args:      $release1, $release2, the first and last wormbase release
#                  numbers to use e.g. 140, 150 to convert data made using wormbase
#                  release WS140 to the coordinates of release WS150
#            @mapping_data - data as returned by read_mapping_data
# Returns:   1 if there have been changes, 0 if not


sub remap_test {
  my ($release1, $release2, @mapping_data) = @_;

  my $changed = 0;

  foreach my $release (($release1+1) .. $release2) {
              
    foreach my $chromosome (keys %{$mapping_data[$release]}) {
      if (exists $mapping_data[$release]{$chromosome}) {
	foreach  my $fields (@{$mapping_data[$release]{$chromosome}}) {

# The mismatch_start value is the start of the mismatch, it is the first position which doesn't match.
# The mismatch_end value is the base past the end of the mismatch region, the first base which matches again
# ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped)
	  my ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped) = @$fields;

	  # if a region has been flipped in orientation then we will need to remap anythng in that region
	  if ($flipped) {$changed = 1;}

	  # if a region has had an indel then we will need to remap anything downstream of it
	  if ($len1 != $len2)  {$changed = 1;}

	}
      }
    }
  }

  return $changed;
}


##########################################################
# 
# Name:      remap_ace
# Usage:     ($new_start, $new_end, $indel, $change) = remap_ace($chromosome, $start, $end, $release1, $release2, @mapping_data);
# Function:  does the remapping of a pair of location values for an ACE file
# Args:      $chromosome, the chromosome number, e.g. 'III'
#            $start, the start value of the chromosomal location coordinate
#            $end, the end value of the chromosomal location coordinate ($end < $start for reverse sense)
#            $release1, $release2, the first and last wormbase release
#                  numbers to use e.g. 140, 150 to convert data made using wormbase
#                  release WS140 to the coordinates of release WS150
#            @mapping_data - data as returned by read_mapping_data
# Returns:   $new_start, $new_end, $new_sense - the updated location coordinates
#            $indel - true if indels affect this location
#            $change - true if any sort of changes affect this location
#

sub remap_ace {
  my ($chromosome, $start, $end, $release1, $release2, @mapping_data) = @_;
                      
  my $sense = "+";
  my ($indel, $change);

  if ($start > $end) {
    $sense = "-";
    my $tmp = $start;
    $start = $end;
    $end = $tmp;
  }

  ($start, $end, $sense, $indel, $change) = remap_gff($chromosome, $start, $end, $sense, $release1, $release2, @mapping_data);

  if ($sense eq '-') {
    my $tmp = $start;
    $start = $end;
    $end = $tmp;
  }
                                                   
  return ($start, $end, $indel, $change);

}

##########################################################
# 
# Name:      remap_gff
# Usage:     ($new_start, $new_end, $new_sense, $indel, $change) = remap_gff($chromosome, $start, $end, $sense, $release1, $release2, @mapping_data);
# Function:  does the remapping of a pair of location values for a GFF file
# Args:      $chromosome, the chromosome number, e.g. 'III'
#            $start, the start value of the chromosomal location coordinate
#            $end, the end value of the chromosomal location coordinate (always >= $start)
#            $sense, the sense "+" or "-" of the coordinate
#            $release1, $release2, the first and last wormbase release
#                  numbers to use e.g. 140, 150 to convert data made using wormbase
#                  release WS140 to the coordinates of release WS150
#            @mapping_data - data as returned by read_mapping_data
# Returns:   $new_start, $new_end, $new_sense - the updated chromosomal location coordinates
#            $indel - true if indels affect this location
#            $change - true if any sort of changes affect this location
#

sub remap_gff {
  my ($chromosome, $start, $end, $sense, $release1, $release2, @mapping_data) = @_;

  my $indel = 0;		# true if indels affect this location
  my $change = 0;		# true if non-indel base changes affect this location
                                                                                                                                                            
  if ($chromosome =~ /CHROMOSOME_(\S+)/) {$chromosome = $1;}
                                                                                                                                                            
  foreach my $release (($release1+1) .. $release2) {
                                                                                                                                                            
    if (exists $mapping_data[$release]{$chromosome}) {
      foreach  my $fields (@{$mapping_data[$release]{$chromosome}}) {
        #print "$release $chromosome fields= @$fields \n";

# The mismatch_start value is the start of the mismatch, it is the first position which doesn't match.
# The mismatch_end value is the base past the end of the mismatch region, the first base which matches again
# ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped)
        my ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped) = @$fields;

# N.B. mismatch values are in the normal perl coordinate system starting at position 0.
# Convert them to the GFF coordinates system starting at position 1.
	$mismatch_start1++;
	$mismatch_end1++;
                                                                                                                                                            
        if ($flipped) {    # is the feature inside a flipped region?
                                                                                                                                                            
          if ($start >= $mismatch_start1 && $end < $mismatch_end1) {
                                                                                                                                                            
            if ($sense eq '+') {$sense = '-';} else {$sense = '+';} # flip the sense
            $start = $mismatch_start1 + $mismatch_end1 - $start; # flip the start and end positions
            $end = $mismatch_start1 + $mismatch_end1 - $end;
            if ($start > $end) {
              my $tmp = $start;
              $start = $end;
              $end = $tmp;
            }

          # does the edge of the flipped region overlap our location?
          } elsif ($start >= $mismatch_start1 && $start < $mismatch_end1 ||		
		   $end >= $mismatch_start1 && $end < $mismatch_end1) {
            # don't change the location, but note that we have changes
	    $change = 1;                                        
	  }
                                                                                                                                                            
        } else {
	  
	  # if there is a change inside our location, note it
	  if ($start < $mismatch_end1 && $end >= $mismatch_start1) {
	    $change = 1;
	  }

	  # note the length of our location so we can see any indels occurring
	  my $location_length = $end - $start;

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

	  # see if we have any indels in our location
	  if ($location_length != $end - $start) {
	    $indel = 1;
	  }
        }
                                                                                                                                                            
      }
    } else {
      #print "no change: doesn't exist: $release $chromosome\n";
    }
  }
                                                                                                                                                            
                                                                                                                                                            
  return ($start, $end, $sense, $indel, $change);
}
##########################################################
# 
# Name:      unmap_gff
# Usage:     ($new_start, $new_end, $new_sense, $indel, $change) = remap_gff($chromosome, $start, $end, $sense, $release2, $release1, @mapping_data);
# Function:  does the unmapping (remapping backwards to past versions) of a pair of location values for a GFF file
# Args:      $chromosome, the chromosome number, e.g. 'III'
#            $start, the start value of the chromosomal location coordinate
#            $end, the end value of the chromosomal location coordinate (always >= $start)
#            $sense, the sense "+" or "-" of the coordinate
#            $release2, $release1, the last and first wormbase release
#                  numbers to use e.g. 150, 140 to convert data made using wormbase
#                  release WS150 to the coordinates of release WS140
#            @mapping_data - data as returned by read_mapping_data
# Returns:   $new_start, $new_end, $new_sense - the updated chromosomal location coordinates
#            $indel - true if indels affect this location
#            $change - true if any sort of changes affect this location


sub unmap_gff {
  my ($chromosome, $start, $end, $sense, $release2, $release1, @mapping_data) = @_;

  my $indel = 0;		# true if indels affect this location
  my $change = 0;		# true if non-indel base changes affect this location
                                                                                                                                                            
  if ($chromosome =~ /CHROMOSOME_(\S+)/) {$chromosome = $1;}
                                                                                                                                                            
  foreach my $release ($release2 .. ($release1+1)) {
                                                                                                                                                            
    if (exists $mapping_data[$release]{$chromosome}) {
      foreach  my $fields (@{$mapping_data[$release]{$chromosome}}) {
        #print "$release $chromosome fields= @$fields \n";

# The mismatch_start value is the start of the mismatch, it is the first position which doesn't match.
# The mismatch_end value is the base past the end of the mismatch region, the first base which matches again
# ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped)
        my ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped) = @$fields;

# N.B. mismatch values are in the normal perl coordinate system starting at position 0.
# Convert them to the GFF coordinates system starting at position 1.
	$mismatch_start2++;
	$mismatch_end2++;
                                                                                                                                                            
        if ($flipped) {    # is the feature inside a flipped region?
                                                                                                                                                            
          if ($start >= $mismatch_start2 && $end < $mismatch_end2) {
                                                                                                                                                            
            if ($sense eq '+') {$sense = '-';} else {$sense = '+';} # flip the sense
            $start = $mismatch_start2 + $mismatch_end2 - $start; # flip the start and end positions
            $end = $mismatch_start2 + $mismatch_end2 - $end;
            if ($start > $end) {
              my $tmp = $start;
              $start = $end;
              $end = $tmp;
            }

          # does the edge of the flipped region overlap our location?
          } elsif ($start >= $mismatch_start2 && $start < $mismatch_end2 ||		
		   $end >= $mismatch_start2 && $end < $mismatch_end2) {
            # don't change the location, but note that we have changes
	    $change = 1;                                        
	  }
                                                                                                                                                            
        } else {
	  
	  # if there is a change inside our location, note it
	  if ($start < $mismatch_end2 && $end >= $mismatch_start2) {
	    $change = 1;
	  }

	  # note the length of our location so we can see any indels occurring
	  my $location_length = $end - $start;

          # if the start or end are beyond the start of the change region, apply any shift
          if ($start >= $mismatch_end2) { # if past the end of the change region, shift it
            $start += $len1 - $len2;
          } elsif ($start >= $mismatch_start2 && $start - $mismatch_start2 > $len1 ) { # if was in the change region and now out, set it to the end
            $start = $mismatch_start2 + $len1;
          }
                                                                                                                                                            
          if ($end >= $mismatch_end2) { # if past the end of the change region, shift it
            $end += $len1 - $len2;
          } elsif ($end >= $mismatch_start2 && $end - $mismatch_start2 > $len1) { # if was in the change region and now out, set it to the end
            $end = $mismatch_start2 + $len1;
          }

	  # see if we have any indels in our location
	  if ($location_length != $end - $start) {
	    $indel = 1;
	  }
        }
                                                                                                                                                            
      }
    } else {
      #print "no change: doesn't exist: $release $chromosome\n";
    }
  }
                                                                                                                                                            
                                                                                                                                                            
  return ($start, $end, $sense, $indel, $change);
}

##########################################################
# 
# Name:      remap_clone
# Usage:     ($clone_id, $start, $end, $indel, $change) = 
#	Remap_Sequence_Change::remap_clone($wormbase, $clone_id, $start, $end, $version, $current_converter, $autoace_converter, @mapping_data);
# Function:  does the remapping of a pair of ACE clone location values between the autoace version and the currentDB version
# Args:      $clone_id, the name of the clone
#            $start, the start value of the clone location coordinate
#            $end, the end value of the clone location coordinate (always >= $start)
#            $version, the lastest wormbase release (the autoace version number)
#            $current_converter, $autoace_converter - Coords converter objects
#              for currentDB and for autoace
#
#            @mapping_data - data as returned by read_mapping_data
# Returns:   $new_start, $new_end - the updated clone location coordinates
#            $indel - true if indels affect this location
#            $change - true if any sort of changes affect this location
#

sub remap_clone {
  my ($wormbase, $clone_id, $start, $end, $version, $current_converter, $autoace_converter, @mapping_data) = @_;

  my $indel = 0;		# true if indels affect this location
  my $change = 0;		# true if non-indel base changes affect this location

  my $prev_version = $version - 1;
  my ($new_chrom_start, $new_chrom_end);

# convert the clone location to the chromosomal location for version1
  my ($chromosome, $chrom_start, $chrom_end) = clone_to_chromosome($wormbase, $clone_id, $start, $end, $current_converter);
  #print "clone to chromosome: $chromosome, $chrom_start, $chrom_end\n";

# remap the chromosomal location
  ($new_chrom_start, $new_chrom_end, $indel, $change) = remap_ace($chromosome, $chrom_start, $chrom_end, $prev_version, $version, @mapping_data);
  #print "remapped chromosome: $chrom_start, $chrom_end to $new_chrom_start, $new_chrom_end, $indel, $change\n";# if ($new_chrom_start != $chrom_start);

# convert the chromosomal location for version2 to the clone location
  my ($new_clone, $new_start, $new_end) = chromosome_to_clone($wormbase, $chromosome, $new_chrom_start, $new_chrom_end, $autoace_converter);
  #print "chromosome to clone: $new_clone, $new_start, $new_end\n";# if ($new_chrom_start != $chrom_start);

  return ($new_clone, $new_start, $new_end, $indel, $change);
                                                     
}

##########################################################
# 
# Name:      clone_to_chromosome
# Usage:     ($chrom, $chrom_start, $chrom_end) = clone_to_chromosome($wormbase, $clone_id, $start, $end, $current_converter);
# Function:  converts clone coords to chromosomal coords
# Args:      $clone_id, the name of the clone
#            $start, the start value of the clone location coordinate
#            $end, the end value of the clone location coordinate (always >= $start)
#            $current_converter - coords converter object for the current database
# Returns:   $chrom, chromosome
#            $new_start, $new_end, $new_sense - the updated clone location coordinates
#            $indel - true if indels affect this location
#            $change - true if any sort of changes affect this location
#


sub clone_to_chromosome {
  my ($wormbase, $clone_id, $start, $end, $current_converter) = @_;
  my ($chrom, $chrom_start, $chrom_end);


  ($chrom, $chrom_start) = $current_converter->Coords_2chrom_coords($clone_id, $start);
  ($chrom, $chrom_end) = $current_converter->Coords_2chrom_coords($clone_id, $end);

  return ($chrom, $chrom_start, $chrom_end);
}

##########################################################
# 
# Name:      chromosome_to_clone
# Usage:     ($clone, $start, $end, $sense) = chromosome_to_clone($wormbase, $chromosome, $chrom_start, $chrom_end, $autoace_converter);
# Function:  converts  chromosomal coords to clone coords
# Args:      $chromosome, the chromosome number
#            $start, the start value of the chromosome location coordinate
#            $end, the end value of the chromosome location coordinate (always >= $start)
#            $autoace_converter - coords converter object for the autoace database
# Returns:   $chrom, chromosome
#            $new_start, $new_end
#

sub chromosome_to_clone {
  my ($wormbase, $chromosome, $start, $end, $autoace_converter) = @_;
  my ($clone_id, $clone_start, $clone_end);

  return $autoace_converter->LocateSpan($chromosome, $start, $end);
  
}



##########################################################
# 
# Name:      write_changes
# Usage:     $changes = write_changes($release, @mapping_data);
# Function:  returns the chromosomal changes to put in the release notes
# Args:      $release - the release to display changes for
#            @mapping_data - data as returned by read_mapping_data
# Returns:   text of the changes in a readable form
#

sub write_changes {
  my ($release, @mapping_data) = @_;
 
  my $text;
  my $title = "Chromosomal Changes:\n--------------------\n";
  my $no_changes = "There are no changes to the chromosome sequences in this release.\n";
  my $any_changes = 0;
  my @chromosomes = qw(I II III IV V X);
           
  foreach my $chromosome (@chromosomes) {
                                                                                                                                                            
    if (exists $mapping_data[$release]{$chromosome}) {

      $any_changes = 1;
      $text .= "\nChromosome: $chromosome\n";

      foreach  my $fields (@{$mapping_data[$release]{$chromosome}}) {
# The mismatch_start value is the start of the mismatch, it is the first position which doesn't match.
# The mismatch_end value is the base past the end of the mismatch region, the first base which matches again
# ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped)
        my ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped) = @$fields;

# N.B. mismatch values are in the normal perl coordinate system starting at position 0.
# Convert them to the human-readable coordinates system starting at position 1.
	$mismatch_start1++;
#	$mismatch_end1++; # don't increment this so it becomes the end of the changed region
	$mismatch_start2++;
#	$mismatch_end2++; # don't increment this so it becomes the end of the changed region
             
	$text .= "$mismatch_start1 $mismatch_end1 $len1   ->   $mismatch_start2 $mismatch_end2 $len2";
                                                                                                                                               
        if ($flipped) {
	  $text .= " REVERSED";
	}
	
	$text .= "\n";

      }
    }
  }

  if (!$any_changes) {
    $text = $no_changes;
  }

  return $title . $text;
}


1;

__END__

=pod

=head2 NAME - Remap_Sequence_Change.pm

=head2 USAGE

 use Remap_Sequence_Change;
 Remap_Sequence_Change::Function()

=over 4

=head2 DESCRIPTION 

Module for remapping chromosomal sequence locations across releases.

=back

=head3 FUNCTIONS

=over 4

=item @mapping_data = &read_mapping_data($release1, $release2);

Reads in the mapping information from the files created during the
build of each release.

=back

=over 4

=item ($new_start, $new_end, $indel, $change) = remap_ace($chromosome, $start, $end, @mapping_data);

Maps a start and end location pair of ACE file chromosomal
coordinates, to the new coordinates.

=back

=item ($new_start, $new_end, $new_sense, $indel, $change) = remap_ace($chromosome, $start, $end, $sense, @mapping_data);

Maps a start and end location pair of GFF file chromosomal
coordinates, with their sense ("+" or "-") to the new coordinates.

=back

=head3 AUTHOR
$Author: gw3 $

=cut
