use strict;
use warnings;
use Test::More;
use ProductionMysql;
use List::Util qw/pairs unpairs/;
use Bio::EnsEMBL::DataCheck::Checks::RepeatFeatures;
use Bio::EnsEMBL::DataCheck::Checks::AnalysisDescription;
use Bio::EnsEMBL::DataCheck::Checks::PepstatsAttributes;
use Bio::EnsEMBL::DataCheck::Checks::GenomeStatistics;
use Bio::EnsEMBL::DataCheck::Checks::ProteinFeatures;
our @SPECIES = ProductionMysql->staging->species(@ARGV ? @ARGV :  "_core_$ENV{PARASITE_VERSION}_");
sub run_datacheck_for_all_species {
  my ($test_name, $class) = @_;
  subtest $test_name => sub {
    for my $species (@SPECIES){
      subtest $species => sub { 
        my $o = bless { dba => ProductionMysql->staging->adaptor($species, "MetaContainer")->{db}}, $class;
        $o->tests;
      };
    }
  };
}
run_datacheck_for_all_species("DNAFeatures", 'Bio::EnsEMBL::DataCheck::Checks::RepeatFeatures');
# TODO RNAFeatures
run_datacheck_for_all_species("Metadata and analysis descriptions", 'Bio::EnsEMBL::DataCheck::Checks::AnalysisDescription');
subtest "Core statistics" => sub {
  run_datacheck_for_all_species("Pepstats", 'Bio::EnsEMBL::DataCheck::Checks::PepstatsAttributes');
  run_datacheck_for_all_species("GenomeStatistics", 'Bio::EnsEMBL::DataCheck::Checks::GenomeStatistics');
};
run_datacheck_for_all_species("Protein features", 'Bio::EnsEMBL::DataCheck::Checks::ProteinFeatures');

my $s = <<EOF;
>load species
species.taxonomy_id
species.url
genebuild.start_date
genebuild.version
provider.name
provider.url
>BUSCO
assembly.busco_complete
assembly.busco_duplicated
assembly.busco_fragmented
assembly.busco_missing
assembly.busco_number
>CEGMA
assembly.cegma_complete
assembly.cegma_partial
>Xrefs
xref.timestamp
>Example gene
sample.gene_param
sample.gene_text
sample.location_param
sample.location_text
sample.search_text
sample.transcript_param
sample.transcript_text
>Example VEP
sample.vep_ensembl
sample.vep_hgvs
sample.vep_pileup
sample.vep_vcf
EOF
my @meta_keys_by_section = map {
  my ($h, @meta_keys) = split "\n";
  $h ? [$h => \@meta_keys] : ()
} split ">", $s;
my %results;
for my $species (@SPECIES){
  my $meta_adaptor = ProductionMysql->staging->adaptor($species, "MetaContainer");
  for my $p ( @meta_keys_by_section){
    for my $meta_key (@{$p->[1]}){
       my $vs = $meta_adaptor->list_value_by_key($meta_key);
       $results{$meta_key}{$species} = scalar @$vs;
    }
  }
}
for my $p (@meta_keys_by_section){
  my %missing;
  for my $meta_key (@{$p->[1]}){
    my @missing_species_for_key = map {$_->[0]} grep {not $_->[1]} pairs %{$results{$meta_key}};
    $missing{$meta_key} = \@missing_species_for_key if @missing_species_for_key;
  }
  ok (not(%missing), $p->[0]) or diag explain \%missing;
}

done_testing;
