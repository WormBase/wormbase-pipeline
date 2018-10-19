
package GenomeBrowser::Resources::StudyAttributes;
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
            my $properties = $data_for_study->{properties};
            my $bioproject_id = $properties->{BioProject};
            if($bioproject_id and not $properties->{submitting_centre}){
              my $bioproject_properties = &_properties_for_study_from_ena_bioproject_xml(
                    $class->get_xml("https://www.ebi.ac.uk/ena/data/view/$bioproject_id&display=xml")
              );
              $properties = {%$properties, %$bioproject_properties};
            }
            $data{$assembly}{$study_id} = {%$data_for_study, properties => $properties};
        }
    }
    return \%data;
}
sub _clean_messy_text {
   my $name = shift;
   $name =~s/_+/ /g if scalar(split " ", $name) == 1;
   $name = ucfirst(lc($name)) if $name eq uc($name);
   return $name ? $name : ();
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
sub _properties_for_study_from_ena_bioproject_xml {
    my $payload = shift;
    my $properties = {};
    return $properties unless $payload;

    my $submitting_centre = $payload->{PROJECT}{center_name};
    $properties->{submitting_centre} = $submitting_centre if $submitting_centre;
    return $properties;
}
sub _data_for_study_from_ena_study_xml {
    my $payload = shift;
    return {} unless $payload;
    my $properties = {
        "ENA first public" => join( " ",
            map { $_->{TAG} eq 'ENA-FIRST-PUBLIC' ? $_->{VALUE} : () }
              @{ $payload->{STUDY}{STUDY_ATTRIBUTES}{STUDY_ATTRIBUTE} } ),
        "ENA last update" => join( " ",
            map { $_->{TAG} eq 'ENA-LAST-UPDATE' ? $_->{VALUE} : () }
              @{ $payload->{STUDY}{STUDY_ATTRIBUTES}{STUDY_ATTRIBUTE} } ),
    };

    $properties->{submitting_centre} = $payload->{STUDY}{center_name} 
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
    my ($bp, @other_bioprojects) = @bioprojects;
    die %$properties if @other_bioprojects;
    $properties->{BioProject} = $bp if $bp;
    my @pubmed_refs;
    for my $study_link (@{ $payload->{STUDY}{STUDY_LINKS}{STUDY_LINK} }){
       if(uc($study_link->{XREF_LINK}{DB}) eq 'PUBMED') {
          push @pubmed_refs, $study_link->{XREF_LINK}{ID};
       } elsif ($study_link->{URL_LINK}{LABEL} and $study_link->{URL_LINK}{URL}){
          my ($k, $v) = _url_link($study_link->{URL_LINK}{LABEL}, $study_link->{URL_LINK}{URL});
          $properties->{$k} = $v;
       } else {
          #probably a link to an ENA something - skip
       }
    }    
    $properties->{"PubMed"} = join(", ", map {_pubmed_link($_)} @pubmed_refs) if @pubmed_refs;
    
    my $sd = $payload->{STUDY}{DESCRIPTOR}{STUDY_DESCRIPTION};
    $properties->{"Study description"} = $sd if $sd and length($sd) < 500 and $sd !~ /This data is part of a pre-publication release/;
    my $st = _clean_messy_text($payload->{STUDY}{DESCRIPTOR}{STUDY_TITLE});
    return {
      properties => $properties,
      ($st ? (study_title => $st) : ()),
    };
}
1;
