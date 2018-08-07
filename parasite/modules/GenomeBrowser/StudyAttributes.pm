
package GenomeBrowser::StudyAttributes;
use List::Util qw(sum);
use List::MoreUtils qw(uniq);
use parent GenomeBrowser::LocallyCachedResource;

sub _fetch {
    my ( $class, $species, $assembly, $rnaseqer_metadata ) = @_;
    my %data;
    for my $study_id ( @{ $rnaseqer_metadata->access($assembly) } ) {
        $data{$study_id} = &_properties_for_study_from_ena_payload(
            $class->get_xml(
                "https://www.ebi.ac.uk/ena/data/view/$study_id&display=xml")
        );
    }
    return \%data;
}
sub _clean_messy_text {
   my $name = shift;
   $name =~s/_+/ /g if scalar(split " ", $name) == 1;
   $name = ucfirst(lc($name)) if $name eq uc($name);
   return $name;
}
sub _pubmed_link {
    my $id = shift;
return "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/$id\">$id</a>";
}
sub _url_link {
   my ($label, $url) = @_;

#Are there any more URL links?
   (my $property_name = $label) =~ s/.*ArrayExpress/ArrayExpress/;

   return ($property_name , "<a href=\"$url\">$label</a>");
}
sub _properties_for_study_from_ena_payload {
    my $payload = shift;
    return {} unless $payload;
    my $result = {
        "study" => join( ": ",
            $payload->{STUDY}{IDENTIFIERS}{PRIMARY_ID},
            _clean_messy_text($payload->{STUDY}{DESCRIPTOR}{STUDY_TITLE})
        ),
        "ENA first public" => join( " ",
            map { $_->{TAG} eq 'ENA-FIRST-PUBLIC' ? $_->{VALUE} : () }
              @{ $payload->{STUDY}{STUDY_ATTRIBUTES}{STUDY_ATTRIBUTE} } ),
        "ENA last update" => join( " ",
            map { $_->{TAG} eq 'ENA-LAST-UPDATE' ? $_->{VALUE} : () }
              @{ $payload->{STUDY}{STUDY_ATTRIBUTES}{STUDY_ATTRIBUTE} } ),
    };
    
    my @pubmed_refs;
    for my $study_link (@{ $payload->{STUDY}{STUDY_LINKS}{STUDY_LINK} }){
       if(uc($study_link->{XREF_LINK}{DB}) eq 'PUBMED') {
          push @pubmed_refs, $study_link->{XREF_LINK}{ID};
       } elsif ($study_link->{URL_LINK}{LABEL} and $study_link->{URL_LINK}{URL}){
          my ($k, $v) = _url_link($study_link->{URL_LINK}{LABEL}, $study_link->{URL_LINK}{URL});
          $result->{$k} = $v;
       } else {
          #probably a link to an ENA something - skip
       }
    }    
    $result->{"PubMed"} = join(", ", map {_pubmed_link($_)} @pubmed_refs) if @pubmed_refs;
    
    my $sd = $payload->{STUDY}{DESCRIPTOR}{STUDY_DESCRIPTION};
    $result->{"Study description"} = $sd if $sd and length($sd) < 500 and $sd !~ /This data is part of a pre-publication release/;
    return $result;
}
1;
