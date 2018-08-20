package GenomeBrowser::Links;
use List::Util qw(sum);
use List::MoreUtils qw(uniq);
use parent GenomeBrowser::LocallyCachedResource;

sub _link {
   my ($url, $name) = @_;
   return "<a href=\"$url\">$name</a>";
}
sub _fetch {
    my ( $class, $species, $rnaseqer_metadata ) = @_;
    my %data;
    for my $assembly (@{ $rnaseqer_metadata->access}){
        for my $study_id ( @{ $rnaseqer_metadata->access($assembly) } ) {
            for my $run_id (@{ $rnaseqer_metadata->access($assembly, $study_id)}) {
                $data{$run_id}{"Mapping results"} =
                   _link($rnaseqer_metadata->data_location($run_id), "RNASeq-er processing directory: $run_id");
                $data{$run_id}{"ENA study"} =
                   _link("http://www.ebi.ac.uk/ena/data/view/$study_id", "Study page: $study_id");
            }
        }
    }
    return \%data;
}
1;
