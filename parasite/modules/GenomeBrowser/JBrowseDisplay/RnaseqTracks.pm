use strict;
use warnings;
package GenomeBrowser::JBrowseDisplay::RnaseqTracks;

use ProductionMysql;
use List::MoreUtils qw(uniq);
use List::Util qw(sum pairmap);
use File::Slurp qw/read_file/;
use JSON;
use GenomeBrowser::Deployment;
use Log::Any '$log';

# This is the data folder for jbrowse consumption
#
# input parameters:
#  - where to construct the folder
#  - corresponding data production location
sub new {
    my ( $class, $wbps_expression_dir ) = @_;
    return bless {
       wbps_expression_dir => $wbps_expression_dir
    }, $class;
}

sub track_selector_and_tracks_for_species_and_assembly {
    my ($self, $species, $assembly, %opts) = @_;
    my ($spe, $cies, $bp ) = split "_", $species;
    my $path = join("/", $self->{wbps_expression_dir}, $species, "${spe}_${cies}.studies.json");
    return unless -f $path;
    $log->info(__PACKAGE__ . "::track_selector_and_tracks_for_species_and_assembly $path");
    my @studies = @{from_json(read_file($path, { binmode => ':utf8' }))};
    return unless @studies;
    my @rnaseq_track_configs;
    my $has_multiple_categories = 1 < uniq map {$_->{study_category}} @studies;
    for my $study (@studies) {
        my $category = $has_multiple_categories ? "RNASeq: $study->{study_category}" : "RNASeq";
        for my $run (@{$study->{runs}}) {
          my $run_id = $run->{run_id};
          my $url    = GenomeBrowser::Deployment::sync_ebi_externally(
              $species, $assembly, $run_id,
              $run->{bigwig}, %opts );
          my $track_name = join(": ", grep {$_} (join("/", grep {$_} $run_id, $run->{replicate})), $run->{condition});
          my $attributes = {
             %{$study->{attributes}},
             %{$run->{attributes} // {}},
             track => $track_name,
             study => sprintf("%s: %s", $study->{study_id}, $study->{study_title}),
          };
          push @rnaseq_track_configs,
            {
              storeClass    => "JBrowse/Store/SeqFeature/BigWig",
              type          => "JBrowse/View/Track/Wiggle/XYPlot",
              category      => $run->{attributes} ? $category: "$category, unannotated",
              autoscale     => "local",
              ScalePosition => "right",
              urlTemplate => $url,
              key         => $track_name,
              label       => "RNASeq/$run_id",
              metadata    => $attributes,
            };
        }
    }
    return track_selector(column_headers_for_studies(@studies)), @rnaseq_track_configs;
}
sub track_selector {
    my ( $as ) = @_;
    my @as = @{$as};
    my %pretty;
    for my $a (@as) {
        ( my $p = $a ) =~ s/[\W_-]+/ /g;
        $pretty{$a} = ucfirst($p);
    }
    return {
        type             => "Faceted",
        displayColumns   => [ "track", @as ],
        selectableFacets => [
            "category","pubmed", "study",
            "submitting_centre",
            @as
        ],
        renameFacets => {
            pubmed => "PubMed",
            study                  => "Study",
            submitting_centre      => "Submitting centre",
            track  => "Track",
            %pretty
        }
    };
}

sub column_headers_for_studies {
  my $runs_total;
  my %occurrences_per_type;
  for my $study(@_){
    my %rnaseqer_characteristics;
    for my $run (@{$study->{runs}}) {
      $runs_total++;
      $occurrences_per_type{$_}++ for keys %{$run->{attributes}};    
    }
  }
  my @result = sort keys %occurrences_per_type;
  if($runs_total > 100) {
      @result = grep {$occurrences_per_type{$_} * 20 > $runs_total} @result;
  }
  return \@result;
}

1;
