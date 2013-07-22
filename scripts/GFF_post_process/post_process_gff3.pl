#!/usr/bin/env perl
#
# post_process_gff3.pl
#
# Transforms the raw GFF3 Ace dumps into a purer form of GFF3, including:
# - giving the WormBase canonical genes/transcripts/CDSs/ncRNAs source 'WormBase'
# - Correcting coordinates of between-base features for GFF3 conventions
# - Strip the class name prefix from all of the Name attrbutes 
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-07-22 16:02:35 $

use strict;
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

my %source_map = (
  Coding_transcript        => 'WormBase',
  curated                  => 'WormBase',
  gene                     => 'WormBase',
  Pseudogene               => 'WormBase',
  miRNA_mature_transcript  => 'WormBase',
  curated_miRNA            => 'WormBase',
  rRNA                     => 'WormBase',
  snoRNA_mature_transcript => 'WormBase',
  snRNA_mature_transcript  => 'WormBase',
  ncRNA                    => 'WormBase',
  snlRNA                   => 'WormBase',
  scRNA                    => 'WormBase',
    );

my %between_base_feature_types = (
  SL1                                 => 1,
  SL2                                 => 1,
  polyA_site                          => 1,
  insertion_site                      => 1,
  transposable_element_insertion_site => 1,
    );


######################################
# variables and command-line options #
######################################

my (  $debug, $test, $store, $wormbase, $species );
my ( $infile, $outfile, $changed_lines );

GetOptions(
    'debug=s'   => \$debug,
    'test'      => \$test,
    'store:s'   => \$store,
    'species:s' => \$species,
    'infile:s'  => \$infile,
    'outfile:s' => \$outfile,
);

if ($store) {
    $wormbase = retrieve($store) or croak("Can't restore wormbase from $store\n");
} 
else {
    $wormbase = Wormbase->new(
        -debug => $debug,
        -test  => $test,
        -organism => $species,
    );
}

my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

while(<$gff_in_fh>) {
  /^\#/ and do {
    print $gff_out_fh $_;
    next;
  };
  chomp;
  my $line = $_;
  my @l  = split(/\t+/, $line);

  #
  # Change source
  #
  if (exists $source_map{$l[1]}) {
    $l[1] = $source_map{$l[1]};
  }

  #
  # Correct location of between-base features
  #
  if (exists $between_base_feature_types{$l[2]}) {
    if ($l[7] eq '+') {
      $l[4]--;
    } else {
      $l[3]++;
    }
  }

  #
  # Remove class name prefix from Name attribute
  # 
  my @attr = split(/;/, $l[8]);

  @attr = map { 
    s/Name=\S+:(\S+)/Name=\1/;  
    s/^note=/Note=/;
    $_
  } @attr;
  $l[8] = join(";", @attr);

  my $new_line = join("\t", @l);
  $changed_lines++ if $new_line ne $line;

  print $gff_out_fh "$new_line\n";
} 
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
$log->write_to("Finished processing : $changed_lines lines modified\n");
$log->mail();
exit(0);
