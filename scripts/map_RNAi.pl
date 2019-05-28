#!/usr/bin/env perl
#
# map_RNAi.pl
#
# Associate RNA experiments with Transcripts/Genes based
# on genomic location
#
# Version: $Version: $
# Last updated by: $Author: klh $
# Last updated on: $Date: 2014-10-09 16:03:15 $

use strict;
use warnings;
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Modules::Map_Helper;
use Ace;

###############
# variables   #
###############

my (
  $store,
  $test,
  $debug,
  $species,
  $noload,
  $acefile,
  $ace_fh,
  $gffdir,
  $wormbase,
  $database,
  $min_overlap,
  @chromosomes,
  %results,
    );

$min_overlap = 50;

GetOptions(
  "debug=s"       => \$debug,
  'store=s'       => \$store,
  'species=s'     => \$species,
  "test"          => \$test,
  "noload"        => \$noload,
  "acefile=s"     => \$acefile,
  'chrom=s@'      => \@chromosomes,
  'gffdir=s'      => \$gffdir,
  'database=s'    => \$database,
  'minoverlap=s'  => \$min_overlap,
    );

#################################################
# config - what are we going to compare against?
#################################################
#                        Class          GFF_file_prefix              GFF_source            GFF_type    GFF_name_type
my $to_search_against = {
  
			 Predicted_gene => [ ['curated', 'curated', 'exon', 'CDS'] ],
			 
			 Transcript => [ 
					['7kncRNA',                  '7kncRNA',            'exon'],
					['asRNA',                    'asRNA',              'exon'],
					['Coding_transcript',        'Coding_transcript',  'exon'],
					['lincRNA',                  'lincRNA',            'exon'],
					['miRNA_primary_transcript', 'miRNA_primary_transcript',      'exon'],
					['ncRNA',                    'ncRNA',              'exon'],
					['Non_coding_transcript',    'non_coding_transcript',         'exon'],
					['piRNA',                    'piRNA',              'exon'],
					['rRNA',                     'rRNA',               'exon'],
					['scRNA',                    'scRNA',              'exon'],
					['snlRNA',                   'snRNA',              'exon'], # in GFF files as 'snRNA' - we do not distingush between snRNA and snlRNA in GFF
					['snoRNA',                   'snoRNA',             'exon'],
					['snRNA',                    'snRNA',              'exon'], 
					['stRNA',                    'stRNA',              'exon'], 
					['tRNA',                     'tRNA',               'exon'],
				       ],
			 
};



############################
# recreate configuration   #
############################
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("cant restore wormbase from $store\n"); 
}
else { 
  $wormbase = Wormbase->new(-debug => $debug, 
                            -test => $test,
                            -organism => $species,
      );
}

my $log = Log_files->make_build_log($wormbase);

##############
# Paths etc. #
##############

$database  = $wormbase->autoace    if not $database;
$gffdir    = $wormbase->gff_splits if not $gffdir;
$acefile = $wormbase->acefiles . "/RNAi_mappings.ace" if not defined $acefile;

if (not @chromosomes) {
  @chromosomes = $wormbase->get_chromosome_names(-mito => 1, -prefix => 1);
}


foreach my $class (keys %$to_search_against) {
  $log->write_to("Mapping to $class\n");

  my $fm = Map_Helper->new();

  foreach my $stanza (@{$to_search_against->{$class}}) {
    my ($file_prefix, $gff_source, $gff_type, $name_tag) = @$stanza;
    $name_tag = $class if not defined $name_tag;

    my @files;
    if ($wormbase->assembly_type eq 'contig') {
      push @files, "${gffdir}/${file_prefix}.gff";
    } else {
      foreach my $chr (@chromosomes) {
        push @files, "${gffdir}/${chr}_${file_prefix}.gff";
      }
    }

    foreach my $file (@files) {
      if (not -e $file) {
        warn("Could not find $file - skipping\n");
        next;
      }
      $log->write_to(" Reading $file\n");
      $fm->populate_from_GFF($file, $gff_source, $gff_type, sprintf('%s \"(\S+)\"', $name_tag));
    }

    $log->write_to(" Building index...\n");
    $fm->build_index();

    $log->write_to(" Searching features...\n");

    &get_RNAi_from_gff( \%results, $fm, $class);
  }
}


########################
# produce output files #
########################
my %tran2gene = $wormbase->FetchData('worm_gene2geneID_name');
open($ace_fh, ">$acefile") or $log->log_and_die("Could not open $acefile for reading\n");

my (%gene2rnai, %assoc_counts);

foreach my $rnai ( keys %results ) {
  my %genes;

  print $ace_fh "\nRNAi : $rnai\n";
  foreach my $rnai_type (keys %{$results{$rnai}}) {
    foreach my $class (keys %{$results{$rnai}->{$rnai_type}}) {
      foreach my $obj (keys %{$results{$rnai}->{$rnai_type}->{$class}}) {
        print $ace_fh "$class \"$obj\" Inferred_automatically \"RNAi_${rnai_type}\"\n";
        if (not exists $tran2gene{$obj}) {
          $log->log_and_die("Could not find parent gene for $obj\n");
        }
        my $gene = $tran2gene{$obj};
        $genes{$gene}->{$rnai_type} = 1;
        $gene2rnai{$gene}->{$rnai}->{$rnai_type} = 1;
      }
    }
  }
  foreach my $gene (keys %genes) {
    foreach my $tp (keys %{$genes{$gene}}) {
      print $ace_fh "Gene \"$gene\" Inferred_automatically \"RNAi_${tp}\"\n"; 
    }
  }
}



foreach my $gene ( keys %gene2rnai) {
  print $ace_fh "\nGene : $gene\n";
  foreach my $rnai (sort keys %{$gene2rnai{$gene}}) {
    my @types = keys %{$gene2rnai{$gene}->{$rnai}};
    foreach my $type (@types) {
      $assoc_counts{$type}++;

      print $ace_fh "Experimental_info RNAi_result \"$rnai\" Inferred_automatically \"RNAi_${type}\"\n";
    }
  }
}

$wormbase->load_to_database( $database, $acefile, "RNAi_mappings", $log )
    unless $noload;
    

####################################
# print some statistics to the log #
####################################

$log->write_to("\n\nStatistics\n");
$log->write_to("----------\n\n");
foreach my $type (sort keys %assoc_counts) {
  my $count = $assoc_counts{$type};
  $log->write_to("No. of $type RNAi links to genes written to database: $count\n");
}

$log->mail();
exit(0);

##############################################################
# Subroutines
##############################################################


#######################################
sub get_RNAi_from_gff {
  my ($results, $map, $class) = @_;

  foreach my $chr ($wormbase->get_chromosome_names(-prefix => 1, -mito => 1)) {
    my $gff_file = "${gffdir}/${chr}_RNAi_primary.gff";

    open(my $rnaifh, $gff_file) or $log->log_and_die("Could not open $gff_file for reading\n");

    while (<$rnaifh>) {
      /^\#/ and next;
      chomp;
      my @l = split(/\t+/, $_);

      next unless $l[1] eq "RNAi_primary" and $l[2] eq 'RNAi_reagent';

      my ($name) = $l[8] =~ /\"RNAi:(\S+.+)\"\s+\d+\s+\d+$/;

      my @hits = @{$map->search_feature_segments($l[0], $l[3], $l[4], undef, $min_overlap)};
      map { $results->{$name}->{'primary'}->{$class}->{$_} = 1 } @hits;
    
    }
  }
}

__END__
