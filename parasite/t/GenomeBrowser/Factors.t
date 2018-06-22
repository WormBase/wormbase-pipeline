use Test::More;

use File::Temp qw/tempdir/;
use GenomeBrowser::ArrayExpressMetadata;
use GenomeBrowser::RnaseqerMetadata;
use GenomeBrowser::Factors;

sub studies_grouped_by_factor_as_expected {
  my ($rnaseqer_characteristics_per_study, $array_express_factors_per_study, $expected, $description) = @_;

  my $root_dir = tempdir( CLEANUP => 1 );
  my $species = "schistosoma_mansoni";
  my $assembly = "schisto_v7.2";
  my $mock_rm = bless {
     $assembly => $rnaseqer_characteristics_per_study
  }, 'GenomeBrowser::RnaseqerMetadata';
  my $mock_aem = bless {
     primary_accession_to_factor_type => $array_express_factors_per_study,
     secondary_to_primary_accession => []
  }, 'GenomeBrowser::ArrayExpressMetadata';
  
  my $subject = GenomeBrowser::Factors->new($root_dir,$species,$assembly,$mock_rm,$mock_aem);

  is_deeply(\{@$subject}, $expected, $description);
}

studies_grouped_by_factor_as_expected({},{},[], "Null case");

done_testing();
