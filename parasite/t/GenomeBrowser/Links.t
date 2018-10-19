use GenomeBrowser::Links;
use Test::More;
my $links = 'GenomeBrowser::Links';
my $result = $links->misc_links("run_id", "study_id", "<data_location_url>", [1234, 5678]);
is_deeply($result,
   {
      "ENA study" => '<a href="https://www.ebi.ac.uk/ena/data/view/run_id">Study page: run_id</a>',
      "Mapping results" => '<a href="<data_location_url>">RNASeq-er processing directory: study_id</a>',
      "PubMed" => '<a href="https://www.ncbi.nlm.nih.gov/pubmed/1234">1234</a>, <a href="https://www.ncbi.nlm.nih.gov/pubmed/5678">5678</a>',
  },
   "unit tests are good"
) or diag explain $result;
done_testing();
