
package GenomeBrowser::Studies;
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
    my @pubmed_refs =
      map { $_->{XREF_LINK}{DB} eq 'PUBMED' ? $_->{XREF_LINK}{ID} : () }
      @{ $payload->{STUDY}{STUDY_LINKS}{STUDY_LINK} };
    $result->{"PubMed references"} = join(", ", @pubmed_refs) if @pubmed_refs;
    
    my $sd = $payload->{STUDY}{DESCRIPTOR}{STUDY_DESCRIPTION};
    $result->{"Study description"} = $sd if $sd and length($sd) < 500 and $sd !~ /This data is part of a pre-publication release/;
    return $result;
}
1;
