package GenomeBrowser::Links;
use List::Util qw(sum);
use List::MoreUtils qw(uniq);
use parent GenomeBrowser::LocallyCachedResource;

sub _fetch {
    my ( $class, $species, $assembly, $rnaseqer_metadata ) = @_;
    my %data;
    for my $study_id ( @{ $rnaseqer_metadata->access($assembly) } ) {
        for my $run_id (@{ $rnaseqer_metadata->access($assembly, $study_id)}) {
            $data{$run_id}{"Link: RNASeq-er analysis results"} =
               $rnaseqer_metadata->data_location($run_id);
            $data{$run_id}{"Link: ENA study page"} =
               "http://www.ebi.ac.uk/ena/data/view/$study_id";
        }
    }
    return \%data;
}
1;
