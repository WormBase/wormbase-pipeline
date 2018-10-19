package GenomeBrowser::Links;

sub _link {
   my ($url, $name) = @_;
   return "<a href=\"$url\">$name</a>";
}
sub _pubmed_link {
   my ($id) = @_;
   return _link("https://www.ncbi.nlm.nih.gov/pubmed/$id", $id);
}

sub misc_links {
   my ($self, $study_id, $run_id, $data_location, $pubmed) = @_;
   
   my $result = {
      "Mapping results" =>
          _link($data_location, "RNASeq-er processing directory: $run_id"),
      "ENA study" =>
          _link("https://www.ebi.ac.uk/ena/data/view/$study_id", "Study page: $study_id"),
   };
   $pubmed //=[];
   $result->{PubMed} = join ", ", map {_pubmed_link($_)} @$pubmed if @$pubmed;
   return $result;
}

1;
