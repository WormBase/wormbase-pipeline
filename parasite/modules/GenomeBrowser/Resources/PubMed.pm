package GenomeBrowser::Resources::PubMed;
use parent GenomeBrowser::Resources::LocallyCachedResource;
use GenomeBrowser::Resources::RnaseqerMetadata;
use GenomeBrowser::Resources::ArrayExpressMetadata;
use List::MoreUtils qw(uniq);
use XML::Simple;

sub _fetch {
    my ( $class, $species, $metadata ) = @_;
    my %data;
    for my $assembly( @{$metadata->{rnaseqer}->access}){
       for my $study_id ( @{ $metadata->{rnaseqer}->access($assembly) } ) {
          my $ena_study_pubmed_ids = $metadata->{ena}{$assembly}{$study_id}{study_pubmed} // [];
          my $ena_bioproject_pubmed_ids = $metadata->{ena}{$assembly}{$study_id}{bioproject_pubmed} // [];
          my $geo_pubmed_ids = $metadata->{geo}{$assembly}{$study_id}{pubmed} // [];
          my $ae_pubmed_ids = $metadata->{array_express}->pubmed($study_id) //[];
          for my $pubmed_id ( uniq(@$ena_study_pubmed_ids,@$ena_bioproject_pubmed_ids, @$geo_pubmed_ids, @$ae_pubmed_ids)){
              next if $pubmed_id eq '2971468'; # Check if PRJNA392315 still refers to this paper in error
              $data{$assembly}{$study_id}{$pubmed_id} = &_short_and_full_paper_description_from_payload($class->get_xml(
                   "https://www.ncbi.nlm.nih.gov/pubmed/$pubmed_id?report=xml&format=text"
              ));
          } 
       }
    }
    return \%data;
}

sub _short_and_full_paper_description_from_payload {
    # PubMed formats this as string to encourage people to use their API
    # We are not encouraged enough, so we're going to parse twice.
    my $payload_string = shift;
    my $payload = XMLin($payload_string);

    my @authors = @{$payload->{MedlineCitation}{Article}{AuthorList}{Author} || [] };
    # Use regex: XML::Simple is being too simple.
    # E.g. 30049782: <ArticleTitle>Transgenerational Effects of Extended Dauer Diapause on Starvation Survival and Gene Expression Plasticity in <i>Caenorhabditis elegans</i>.</ArticleTitle>
    my ($title) = $payload_string =~ /<ArticleTitle>(.*)<\/ArticleTitle>/;
    $title = $title->{content} if ref $title eq 'HASH';
    $title = join "", @$title if ref $title eq 'ARRAY';
    if(ref $title){
       use Data::Dumper;
       die Data::Dumper::Dumper($payload);
    }
    my $authors = $payload->{MedlineCitation}{Article}{AuthorList}{Author};
    my @authors = $authors ? ref $authors eq 'ARRAY' ? @$authors : ($authors) : ();
    my $first_author = @authors[0]->{LastName};
    my $last_author = @authors[-1]->{LastName};
    my $authors = $first_author ? $last_author ne $first_author ? "$first_author & $last_author" : $first_author :  "";
    my $year = $payload->{MedlineCitation}{Article}{Journal}{JournalIssue}{PubDate}{Year};
    my $short_description = "$authors, $year";
    return [$short_description, "$title ($short_description)"];
}
1;
