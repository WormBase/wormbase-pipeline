package GenomeBrowser::Resources::PubMed;
use parent GenomeBrowser::Resources::LocallyCachedResource;
use GenomeBrowser::Resources::RnaseqerMetadata;
use GenomeBrowser::Resources::ArrayExpressMetadata;
use List::MoreUtils qw(uniq);

sub _fetch {
    my ( $class, $species, $metadata ) = @_;
    my %data;
    for my $assembly( @{$metadata->{rnaseqer}->access}){
       for my $study_id ( @{ $metadata->{rnaseqer}->access($assembly) } ) {
          my $ena_pubmed_ids = $metadata->{ena}{$assembly}{$study_id}{pubmed} // [];
          my $geo_pubmed_ids = $metadata->{geo}{$assembly}{$study_id}{pubmed} // [];
          my $ae_pubmed_ids = $metadata->{array_express}->pubmed($study_id) //[];
          for my $pubmed_id ( uniq(@$ena_pubmed_ids, @$geo_pubmed_ids, @$ae_pubmed_ids)){
              $data{$assembly}{$study_id}{$pubmed_id} = &_short_and_full_paper_description_from_payload($class->get_xml(
                   "https://www.ncbi.nlm.nih.gov/pubmed/$pubmed_id?report=xml&format=text"
              ));
          } 
       }
    }
    return \%data;
}

sub _short_and_full_paper_description_from_payload {
    my $payload = shift;
    my $title = $payload->{MedlineCitation}{Article}{ArticleTitle};
    my @authors = @{$payload->{MedlineCitation}{Article}{AuthorList}{Author} || [] };
    my $first_author = @authors[0]->{LastName};
    my $last_author = @authors[-1]->{LastName};
    my $authors = $first_author ? $last_author ne $first_author ? "$first_author & $last_author" : $first_author :  "";
    my $year = $payload->{MedlineCitation}{Article}{Journal}{JournalIssue}{PubDate}{Year};
    my $short_description = "$authors, $year";
    return [$short_description, "$title ($short_description)"];
}
1;
