#!/usr/bin/env perl

use strict;
use warnings;

use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::MethodLinkSpeciesSet;
use Bio::EnsEMBL::Compara::Utils::MasterDatabase;


my $compara_url = "$ENV{'HAL_MASTER_DB_URL'}";
my $source = 'wormbase';
my $url = "$ENV{'HAL_MLSS_URL'}";
my $curr_release = $ENV{"ENSEMBL_VERSION"};
my $collection = $ENV{"HAL_COLLECTION_NAME"};

my $compara_dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->go_figure_compara_dba($compara_url);

my $mlss_adaptor = $compara_dba->get_MethodLinkSpeciesSetAdaptor();
my $species_set_adaptor = $compara_dba->get_SpeciesSetAdaptor();
my $method_adaptor = $compara_dba->get_MethodAdaptor();

my $species_set = $species_set_adaptor->fetch_collection_by_name("collection-$collection");
my $pw_method = $method_adaptor->fetch_by_type('CACTUS_HAL_PW');


my @input_genome_dbs = @{$species_set->genome_dbs};
while (my $gdb1 = shift @input_genome_dbs) {
    foreach my $gdb2 (@input_genome_dbs) {

        my $species_set = $species_set_adaptor->fetch_by_GenomeDBs([$gdb1, $gdb2]);
        if (!defined $species_set) {
            $species_set = Bio::EnsEMBL::Compara::Utils::MasterDatabase::create_species_set([$gdb1, $gdb2]);
            $species_set->first_release($curr_release);
            $species_set_adaptor->store($species_set);
        }

        my $mlss_name = sprintf('%s %s', $species_set->name, $pw_method->display_name);
        my $mlss = Bio::EnsEMBL::Compara::MethodLinkSpeciesSet->new(
            -ADAPTOR => $mlss_adaptor,
            -METHOD => $pw_method,
            -NAME => $mlss_name,
            -SOURCE => $source,
            -SPECIES_SET => $species_set
        );

        # The url currently needs to be set directly, in order to circumvent a bug
        # whereby the url value is resolved to the full path before being stored.
        $mlss->{'original_url'} = $url;
        $mlss->{'url'} = $url;

        $mlss->first_release($curr_release);
        $mlss_adaptor->store($mlss);
    }
}
