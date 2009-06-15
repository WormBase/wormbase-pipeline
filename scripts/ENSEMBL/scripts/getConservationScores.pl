#!/bin/env perl

# based on the EnsEMBL compara ConservationScore example
# it dumps all toplevel sequences (chromosomes/contigs) into a bz2 compressed WIG file

# EXAMPLE_USAGE: get_ConservationScores.pl -species remanei -outfile /tmp/remanei.wig.bz2
# TAGS: perl510_clean

BEGIN {
	$ENV{ENSEMBL_REGISTRY}='/nfs/team71/worm/mh6/wormbase/scripts/ENSEMBL/etc/E_registry.pl';
}

use strict;
use lib '/software/worm/ensembl/ensembl/modules';
use lib '/software/worm/ensembl/ensembl-compara/modules';
use lib '/software/worm/ensembl/bioperl-live';
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);

use lib $ENV{'CVS_DIR'};
use Wormbase;

use Getopt::Long;

Bio::EnsEMBL::Registry->load_all;

my ($ofile,$wbspecies,$test,$store,$debug);
GetOptions(
    "outfile:s" => \$ofile,
    "species:s" => \$wbspecies,
    'test'      => \$test,
    'store:s'   => \$store,
    'debug:s'   => \$debug,
)||die($@);

my $wb;
if ($store) {
    $wb = retrieve($store) or croak("Can't restore wormbase from $store\n");
}
else {
    $wb = Wormbase->new(
        -debug => $debug,
        -test  => $test,
        -organism => $wbspecies
    );
}

my $species = $wb->full_name;
$wbspecies  = $wb->organism unless $wbspecies;

$ofile = $wb->chromosomes."/$wbspecies.wig.bz2" unless $ofile;
open OFH ,"|bzip2 -9> $ofile" || die($@);

# get method_link_species_set adaptor
my $mlss_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Compara', 'compara', 'MethodLinkSpeciesSet');

# get the method_link_species_set object for GERP_CONSERVATION_SCORE for 6 species
my $mlss = $mlss_adaptor->fetch_by_method_link_type_registry_aliases("GERP_CONSERVATION_SCORE", 
    ['Caenorhabditis elegans', 'Caenorhabditis briggsae',
     'Caenorhabditis remanei', 'Brugia malayi', 
     'Pristionchus pacificus', 'Caenorhabditis brenneri',
     'Caenorhabditis japonica']);


foreach my $seq_region ($wb->get_chromosome_names(-prefix => 1)){

    #get slice adaptor for $species
    my $slice_adaptor         = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'Slice')
       || throw("Registry configuration file has no data for connecting to <$species>");

    #create slice 
    my @_slice                = $slice_adaptor->fetch_by_region('toplevel', $seq_region);
    throw("No Slice can be created with coordinates $seq_region") if (!@_slice);

    print OFH 'track type=wiggle_0 name="GERP scores" description="single nucleotide ',
              'conservation scores from GERP/PECAN of 7 Nematodes" ',
              'visibility=full color=255,0,0 windowingFunction=mean smoothingWindow=5',"\n",
              "variableStep chrom=$seq_region span=1\n";

    my $window                = 1e6;
    my $slices                = split_Slices(\@_slice, $window, 0);

    #get conservation score adaptor
    my $cs_adaptor            = Bio::EnsEMBL::Registry->get_adaptor('Compara', 'compara', 'ConservationScore');			

    my $count                 = 1;
    while (my $slice          = shift @$slices){

        # To get one score per base in the slice, must set display_size to the size of
        # the slice.
        my $display_size         = $slice->end - $slice->start + 1; 
        my $scores               = $cs_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $display_size);

        print STDERR "processing $seq_region chunk ",$count++,"\n";

        my $start                = $slice->start-1;
        foreach my $score (@$scores) {
            if (defined $score->diff_score) {
                printf OFH ("%i %f\n",$score->position+$start,$score->expected_score - $score->observed_score);
            }
        }
    }
}
