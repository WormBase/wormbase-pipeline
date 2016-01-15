#!/usr/local/bin/perl5.8.0 -w
#
# make_isoforms.pl
# 
# by Gary Williams                         
#
# This makes isoforms in the region of genes
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2015-04-27 11:08:08 $      

# Things isoformer gets confused by or misses:
# - non-canonical spliced introns where the RNASeq intron is placed on the positive strand and so is missing from reverse-strand genes
# - retained introns
# - isoforms that start in the second (or more) exon
# - isoforms that terminate prematurely
# - two or more separate sites of TSL in the same intron cause multiple identical structuers to be created
# - existing structures where the Sequence span is not the same as the span of the exons so a new structure look unique when compared to it


use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
use Sequence_extract;
use Coords_converter;
use Modules::Isoformer;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($species, $database, $gff, $notsl, $outfile, $dogene, $chromosome);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "database:s" => \$database, # database being curated
	    "gff:s"      => \$gff, # optional location of the GFF file if it is not in the normal place
	    "notsl"      => \$notsl, # don't try to make TSL isoforms - for debugging purposes
	    "outfile:s"  => \$outfile, # output ACE file of isoform structures
	    "dogene:s"   => \$dogene, # for testing, specify one gene ID to process
	    "chromosome:s" => \$chromosome, # only do this chromsome
	    );

# always in debug mode
$debug = $ENV{USER};


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

my $currentdb = $wormbase->database('current');
$species = $wormbase->species;
if (!defined $database && $species eq 'elegans') {$database = $currentdb}
if (!defined $database) {$database = $wormbase->autoace()} # BUILD database for the species

# establish log file.
my $log = Log_files->make_build_log($wormbase);

my $Iso = Isoformer->new($wormbase, $log, $database, $gff, $notsl);

$Iso->interactive(0); # don't want any interactive stuff done at all


##########################
# MAIN BODY OF SCRIPT
##########################

# splices are held as array of splice structures

# splice node hash structure
# int - flag to ignore this splice
# char - seq - chromosome or contig
# int - start
# int - end
# char - sense
# int - score
# int array - list of possible child introns in the structure
# int - pos (current position in list of child nodes)

if (!defined $outfile) {
  $outfile = "/nfs/wormpub/CURATION_DATA/${species}_isoformer.ace";
}

open (ISOFORM, "> $outfile")  or die "cant open $Iso->outfile()\n";

# print out the Methods
$Iso->load_isoformer_method;
print ISOFORM $Iso->aceout();
$Iso->aceclear();

if (defined $chromosome) {
  process_genes_in_sequence($chromosome);
} else {
  my @chromosomes = $wormbase->get_chromosome_names(-mito =>1, -prefix => 1);
  foreach my $chromosome (@chromosomes) {
    process_genes_in_sequence($chromosome);
  }
}

print ISOFORM "\n\n// Finished.\n";

close(ISOFORM);

# close the ACE connection
$Iso->{db}->close;


$log->mail();
print "Finished.\n" if ($verbose);
exit(0);



##############################################################
#
# Subroutines
#
##############################################################

# process genes in the given sequence

sub process_genes_in_sequence {
  my ($sequence) = @_;
  
  print "\nLooking at $sequence\n";
  my $chrom_obj = $Iso->{db}->fetch(Sequence => "$sequence");
  my @genes = $chrom_obj->Gene_child;

  my $biotype;
  my $sense;
  
  foreach my $gene ( @genes ) {
    my $gene_name = $gene->name;
    if (defined $dogene && $dogene ne $gene->name) {next}
    my $gene_start = $gene->right->name;
    my $gene_end = $gene->right->right->name;
    my $gene_cds = $gene->Corresponding_CDS;
    my $gene_pseudogene = $gene->Corresponding_pseudogene;
    my $gene_transposon = $gene->Corresponding_transposon;
    my $gene_transcript = $gene->Corresponding_transcript;
    print "Doing $gene_name...\n";
    my $sense = '+';
    if ($gene_start > $gene_end) {
      ($gene_start, $gene_end) = ($gene_end, $gene_start);
      $sense = '-';
    }
    if (defined $gene_cds) {
      $biotype = "CDS";
    } elsif (defined $gene_pseudogene) {
      $biotype = 'Pseudogene';
      my $gene_class = $gene->Gene_class;
      if ($gene_class eq 'rrn') {next} # Ribosomal RNA causes problems because there are thousands of possible isoforms found in these.
    } elsif (defined $gene_transposon) {
      $biotype = 'Transposon';
    } else {
      $biotype = 'Transcript'; # ncRNA
      my $gene_class = $gene->Gene_class;
      if ($gene_class eq 'rrn') {next} # Ribosomal RNA causes problems because there are thousands of possible isoforms found in these.
    }
    #print "sense = $sense biotype = $biotype\n";

    my ($confirmed, $not_confirmed, $created, $warnings) = $Iso->make_isoforms_in_region($sequence, $gene_start, $gene_end, $sense, $biotype, $gene_name);

    print ISOFORM $Iso->aceout();
    $Iso->aceclear();

    if (@{$confirmed}) {print "\n*** The following structures were confirmed: @{$confirmed}\n";}
    if (@{$not_confirmed}) {print "\n*** THE FOLLOWING STRUCTURES WERE NOT CONFIRMED: @{$not_confirmed}\n";}
    if (@{$created}) {print "\n*** The following novel structures were created: @{$created}\n";}
    if (@{$warnings}) {print "\n@{$warnings}\n\n";}
    
  }

  # if this is a two-level reference sequence, get each clone in the chromosome
  my @subsequences = $chrom_obj->Subsequence;
  foreach my $subsequence (@subsequences) {
    my $subsequence_name = $subsequence->name;
    process_genes_in_sequence($subsequence_name);
  }


}
