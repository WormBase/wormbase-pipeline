#!/usr/bin/env perl
#
#   remap_variation_gff_between_releases.pl                 
# 
# This remaps variation GFF between releases
#
# There are a few special cases for variations (mainly insertions), hence custom script
#
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2013-08-07 15:25:31 $      

use strict;                                      
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

use Modules::Remap_Sequence_Change;

######################################
# variables and command-line options # 
######################################

my ($debug, $test, $verbose, $store, $wormbase,$species);
my ($release1, $release2, $genome_diffs_dir);
my ($ignore_list, $only_list);

GetOptions (
  "debug=s"    => \$debug,
  "test"       => \$test,
  "verbose"    => \$verbose,
  "store:s"    => \$store,
  "release1=i" => \$release1,
  "release2=i" => \$release2,
  'species=s'  => \$species,
  'genomediff=s' => \$genome_diffs_dir,
  'ignorelist=s' => \$ignore_list,
  'only=s'       => \$only_list,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species,
			     );
}

$genome_diffs_dir = $wormbase->genome_diffs if not defined $genome_diffs_dir;

my (%only, %ignore);

if ($only_list) {
  open(my $fh, $only_list);
  while(<$fh>) {
    /^(\S+)$/ and $only{$1} = 1;
  }
}
if ($ignore_list) {
  open(my $fh, $ignore_list);
  while(<$fh>) {
    /^(\S+)$/ and $ignore{$1} = 1;
  }
}


if (! defined $release1 || ! defined $release2) {
  die "Specify the release numbers to use\n";
}


##########################
# read in the mapping data
##########################

my $assembly_mapper = Remap_Sequence_Change->new($release1, $release2, $wormbase->species, $genome_diffs_dir);

##########################
# MAIN BODY OF SCRIPT
##########################


while (my $line = <>) {
  chomp $line;
  
  next if $line !~ /\S/;
  next if $line =~ /^\#/;

  my @f = split /\t+/, $line;

  my ($var) = $f[8] =~ /Variation\s+\"(\S+)\"/;

  next if $ignore_list and exists $ignore{$var};
  next if $only_list and not exists $only{$var};

  my ($chromosome, $start, $end, $sense) = ($f[0], $f[3], $f[4], $f[6]);

  $chromosome =~ s/^CHROMOSOME_//;
  if ($wormbase->species eq 'elegans') {
    $chromosome = "CHROMOSOME_${chromosome}";
  }


  
  # some checks for malformed GFF files
  if (not defined $chromosome or
      not defined $start or 
      not defined $end or
      not defined $sense) {
    die("Malformed line (missing values? spaces instead of TABs?): $line\n");
  }

  if ($start !~ /^\d+$/) {
    die("Malformed line (non-numeric start?): $line\n");
  }
  if ($end !~ /^\d+$/) {
    die("Malformed line (non-numeric end?): $line\n");
  }
  if ($sense !~ /^[\+|\-]$/) {
    die("Malformed line (invalid sense?): $line\n");
  }

  my ($new_chr, $new_st, $new_en, $new_strand, $indel, $change, $start_del, $end_del, $dummy1, $dummy2, $dummy3); 

  # map start and end separately so we can see what happened
  ($new_st, $dummy1, $new_strand, $indel, $change, $start_del, $dummy2) 
      = $assembly_mapper->remap_gff($chromosome, $start, $start, $sense);

  ($new_en, $dummy1, $new_strand, $indel, $change, $end_del, $dummy2) 
      = $assembly_mapper->remap_gff($chromosome, $end, $end, $sense);

  if (not $start_del and not $end_del) {
    # start and and map cleanly, but warn if there is a change in length of the feature
    if ($new_en - $new_st != $f[4] - $f[3]) {
      print "# $var - extremities mapped cleanly, but feature ended up different length than before\n";
    } else {
      print "# $var -mapped completely cleanly\n";
    }
  } else {
    # either the start or has been deleted in the version
    if ($start == $end) {
      # single bp feature that's been deleted; we'll have to give up on this one
      print "# $var - single bp feature deleted in new assembly; should be suppressed\n";
    } else {
      my ($new_start_del, $new_end_del);
      if ($start_del) {
        ($new_st, $dummy1, $new_strand, $indel, $change, $new_start_del, $dummy2) 
            = $assembly_mapper->remap_gff($chromosome, $start - 1, $start - 1, $sense);
      }
      if ($end_del) {
        ($new_en, $dummy1, $new_strand, $indel, $change, $new_end_del, $dummy2) 
            = $assembly_mapper->remap_gff($chromosome, $end + 1, $end + 1, $sense);
      }
      if ($new_start_del or $new_end_del) {
        print "# $var - at least one of extremities deleted and could not be rescued by shifting\n";
      } else {
        if ($end - $start != $new_en - $new_st) {
          print "# $var - at least one of extremities deleted AND shifting rescues but end up with weird feature size\n";
        } else {
          print "# $var - at least one of extremities deleted BUT rescued by shifting resulting in clean mapping\n";
        }
      }
    }
  }
  
  $f[3] = $new_st;
  $f[4] = $new_en;
  
  $line = join "\t", @f;
  print $line,"\n";
}

exit(0);



