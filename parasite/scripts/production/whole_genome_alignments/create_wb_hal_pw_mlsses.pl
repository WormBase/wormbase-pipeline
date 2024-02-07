#!/usr/bin/env perl
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 NAME

create_wb_hal_pw_mlsses.pl

=head1 SYNOPSIS

    perl create_wb_hal_pw_mlsses.pl \
        --url <compara_db_url> --release 111 --mlss_id <mlss_id>

=head1 EXAMPLE

    perl create_wb_hal_pw_mlsses.pl \
        --url mysql://ensro@mysql-ps-prod-1.ebi.ac.uk:4450/ensembl_compara_master_parasite --release 111 --mlss_id 7

=head1 OPTIONS

=over

=item B<[--help]>

Prints help message and exits.

=item B<[--url gene_tree_db_url]>

Compara database URL.

=item B<[--release ensembl_version]>

Ensembl release version.

=item B<[--mlss_id mlss_id]>

MLSS ID of the Cactus multiple alignment.

=back

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::Utils::MasterDatabase;


my ( $help, $url, $release, $mlss_id );
GetOptions(
    "help|?"     => \$help,
    "url=s"      => \$url,
    "release=i"  => \$release,
    "mlss_id=i"  => \$mlss_id,
) or pod2usage(-verbose => 2);

pod2usage(-exitvalue => 0, -verbose => 1) if $help;
pod2usage(-verbose => 1) if !$url or !$release or !$mlss_id;


sub get_short_dot_species_set_name {
    my ($species_set) = @_;
    my @individual_names;
    foreach my $gdb (@{$species_set->genome_dbs}) {
        my $species_name = $gdb->name;
        $species_name =~ s/\b(\w)/\U$1/g;
        $species_name =~ s/(\S)\S+\_/$1\./;
        $species_name = substr($species_name, 0, 5);
        push @individual_names, $species_name;
    }
    return join('-', @individual_names);
}


my $compara_dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->go_figure_compara_dba($url);

my $mlss_adaptor = $compara_dba->get_MethodLinkSpeciesSetAdaptor();
my $species_set_adaptor = $compara_dba->get_SpeciesSetAdaptor();
my $method_adaptor = $compara_dba->get_MethodAdaptor();

my $main_mlss = $mlss_adaptor->fetch_by_dbID($mlss_id);
my $hal_url = $main_mlss->get_original_url();

my $pw_method = $method_adaptor->fetch_by_type('CACTUS_HAL_PW');
my $main_species_set = $main_mlss->species_set;
my $source = 'wormbase';


my @input_genome_dbs = @{$main_species_set->genome_dbs};
while (my $gdb1 = shift @input_genome_dbs) {
    foreach my $gdb2 (@input_genome_dbs) {

        my $pw_species_set = $species_set_adaptor->fetch_by_GenomeDBs([$gdb1, $gdb2]);
        if (!defined $pw_species_set) {
            $pw_species_set = Bio::EnsEMBL::Compara::Utils::MasterDatabase::create_species_set([$gdb1, $gdb2]);
            $pw_species_set->first_release($release);
            $pw_species_set->name(get_short_dot_species_set_name($pw_species_set));
            $species_set_adaptor->store($pw_species_set);
        }

        my $pw_mlss = Bio::EnsEMBL::Compara::Utils::MasterDatabase::create_mlss($pw_method, $pw_species_set, $source, $hal_url);
        $pw_mlss->first_release($release);
        $mlss_adaptor->store($pw_mlss);
        $pw_mlss->store_tag('alt_hal_mlss', $mlss_id);
    }
}
