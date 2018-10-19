
package GenomeBrowser::Resources::EnaMetadata;
use List::Util qw(sum);
use List::MoreUtils qw(uniq);
use parent GenomeBrowser::Resources::LocallyCachedResource;

sub _fetch {
    my ( $class, $species, $rnaseqer_metadata ) = @_;
    my %data;
    for my $assembly (@{ $rnaseqer_metadata->access}) {
        for my $study_id ( @{ $rnaseqer_metadata->access($assembly) } ) {
            my $data_for_study = &_data_for_study_from_ena_study_xml(
                $class->get_xml(
                    "https://www.ebi.ac.uk/ena/data/view/$study_id&display=xml")
            );
            if ($data_for_study->{bioproject}){
                $data_for_bioproject = &_data_for_bioproject_from_ena_bioproject_xml(
                    $class->get_xml(sprintf("https://www.ebi.ac.uk/ena/data/view/%s&display=xml", $data_for_study->{bioproject}))
                );
                #TODO more useful stuff: PubMed (uniq sum of both), maybe the descriptions?
                $data_for_study->{attributes}{submitting_centre} //= $data_for_bioproject->{submitting_centre}; 
            }
            $data{$assembly}{$study_id} = $data_for_study;
        }
    }
    return \%data;
}
sub _url_link {
   my ($label, $url) = @_;
# TODO can you rely on BioProjects to make ArrayExpress links?
# Are there any more URL links?
   (my $property_name = $label) =~ s/.*ArrayExpress/ArrayExpress/;

   return ($property_name , "<a href=\"$url\">$label</a>");
}
sub _data_for_bioproject_from_ena_bioproject_xml {
    my $payload = shift;
    return {} unless $payload;
    return {
       submitting_centre => $payload->{PROJECT}{center_name},
    };
}
sub _data_for_study_from_ena_study_xml {
    my $payload = shift;
    return {} unless $payload;
    my $attributes = {
        "ENA first public" => join( " ",
            map { $_->{TAG} eq 'ENA-FIRST-PUBLIC' ? $_->{VALUE} : () }
              @{ $payload->{STUDY}{STUDY_ATTRIBUTES}{STUDY_ATTRIBUTE} } ),
        "ENA last update" => join( " ",
            map { $_->{TAG} eq 'ENA-LAST-UPDATE' ? $_->{VALUE} : () }
              @{ $payload->{STUDY}{STUDY_ATTRIBUTES}{STUDY_ATTRIBUTE} } ),
    };

    $attributes->{submitting_centre} = $payload->{STUDY}{center_name} 
       unless uc($payload->{STUDY}{broker_name}) eq 'NCBI' and length ($payload->{STUDY}{center_name}) < 10;
    my @bioprojects;
    my $ids = $payload->{STUDY}{IDENTIFIERS}{EXTERNAL_ID};
    # XML::Simple is being a bit too simple
    # many ids -> array here: ERP016356
    # one id -> hash here: SRP093920
    my @ids = ref $ids eq 'ARRAY' ? @$ids : ref $ids eq 'HASH' ? ($ids) : ();
    for my $identifier (@ids){
       push @bioprojects, $identifier->{content} if uc($identifier->{namespace}) eq "BIOPROJECT" and $identifier->{label} eq "primary"; 
    }
    my ($bioproject, @other_bioprojects) = @bioprojects;
    die %$attributes if @other_bioprojects;
    my @pubmed_refs;
    for my $study_link (@{ $payload->{STUDY}{STUDY_LINKS}{STUDY_LINK} }){
       if(uc($study_link->{XREF_LINK}{DB}) eq 'PUBMED') {
          push @pubmed_refs, $study_link->{XREF_LINK}{ID};
       } elsif ($study_link->{URL_LINK}{LABEL} and $study_link->{URL_LINK}{URL}){
          my ($k, $v) = _url_link($study_link->{URL_LINK}{LABEL}, $study_link->{URL_LINK}{URL});
          $attributes->{$k} = $v;
       } else {
          #probably a link to an ENA something - skip
       }
    }
    
    return {
      attributes => $attributes,
      bioproject => $bioproject,
      pubmed => \@pubmed_refs,
      study_title => $payload->{STUDY}{DESCRIPTOR}{STUDY_TITLE},
      study_description => $payload->{STUDY}{DESCRIPTOR}{STUDY_DESCRIPTION},
    };
}
1;
