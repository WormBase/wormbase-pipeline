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
# Last updated on: $Date: 2015-02-23 17:20:49 $

use strict;
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

my %source_map = (
  gene                      => 'WormBase',
  Coding_transcript         => 'WormBase',
  curated                   => 'WormBase',
  non_coding_transcript     => 'WormBase',
  Pseudogene                => 'WormBase',
  miRNA_mature              => 'WormBase',
  miRNA_primary_transcript  => 'WormBase',
  pre_miRNA                 => 'WormBase',
  tRNA                      => 'WormBase',
  tRNA_Pseudogene           => 'WormBase',
  rRNA                      => 'WormBase',
  rRNA_Pseudogene           => 'WormBase',
  ncRNA                     => 'WormBase',
  snRNA                     => 'WormBase',
  scRNA                     => 'WormBase',
  stRNA                     => 'WormBase',
  snoRNA                    => 'WormBase',
  lincRNA                   => 'WormBase',
  asRNA                     => 'WormBase',
  piRNA                     => 'WormBase',
  circRNA                   => 'WormBase',
  transposable_element_gene => 'WormBase_transposon',
  Transposon_CDS            => 'WormBase_transposon',
  Transposon_Pseudogene     => 'WormBase_transposon',
  Transposon                => 'WormBase_transposon',
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

my (  $debug, $test, $store, $wormbase, $species, $tecred );
my ( $infile, $outfile, $changed_lines );

GetOptions(
    'debug=s'   => \$debug,
    'test'      => \$test,
    'tecred'      => \$tecred,
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

$species = $wormbase->species;
my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");


#### TECRED modification START ####

# Read in all TECRED IDs and store in a hash
my %tecredid;
if ($tecred) {
	print "Doing TECRED processing\n";
	while(<$gff_in_fh>) {
		chomp;
  		my $line = $_;
		my @l = split(/\t+/, $line);
		if ($l[1]=~/TEC_RED/) {
			my @trid =split(/\s+/,$l[8]);
			$tecredid{$trid[0]}+=1;
			#print "$trid[0]\t$tecredid{$trid[0]}\n";
		}
	}
}
close($gff_in_fh) or $log->log_and_die("Could not close $outfile after writing\n");

#### TECRED modification END ####


open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

while(<$gff_in_fh>) {
  /^\#/ and do {
    if (/sequence-region/) {
      if ($species eq 'elegans') {
        s/CHROMOSOME_//;
      } elsif ($species eq 'briggsae') {
        s/chr//;
      }      
    }
    print $gff_out_fh $_;
    next;
  };
  chomp;
  my $line = $_;
  my @l  = split(/\t+/, $line);

  #
  # Strip prefix
  #
  if ($species eq 'elegans') {
    $l[0] =~ s/^CHROMOSOME_//;
  } elsif ($species eq 'briggsae') {
    $l[0] =~ s/^chr//;
  }

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
    # zero-length features are immediately to the right of start (=end), regardless of strand
    $l[4] = $l[3];
  }

  #
  # Fix up the tecreds if needed https://github.com/WormBase/website/issues/8039
  # Remove orientation, and mark up multimapped TEC-REDs
  #
  if ($tecred) {
	if ($l[1]=~/TEC_RED/) {
  		$l[6]='.';
		#print "$line\n";
		
		# Mark multimapping
		my @trid =split(/\s+/,$l[8]);
		if ($tecredid{$trid[0]} >1) {
			#print "$trid[0]\t multimap\n";
			$l[8]=$l[8] . "; multimapping=TRUE";
		}
		else {
			$l[8]=$l[8] . "; multimapping=FALSE";
			#print "$trid[0]\t$tecredid{$trid[0]}\n";
		}	
	}
  }

  #
  # Clean up attributes
  # 
  my @attr = split(/;/, $l[8]);
  @attr = grep { $_ !~ /Name=Method/ } @attr;

  @attr = map { 
    s/Name=[^:]+:(\S+)/Name=$1/;  
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
