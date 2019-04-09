use strict;
use warnings;
package GenomeBrowser::JBrowseDisplay::RnaseqTracks;

use PublicResources::Rnaseq;
use ProductionMysql;
use List::MoreUtils qw(uniq);
use List::Util qw(sum pairmap);
use List::MoreUtils qw(uniq);

use PublicResources::Rnaseq;
use GenomeBrowser::Deployment;

# This is the data folder for jbrowse consumption
#
# input parameters:
#  - where to construct the folder
#  - corresponding data production location
sub new {
    my ( $class, $resources_dir ) = @_;
    return bless {
       resources => PublicResources::Rnaseq->new($resources_dir)
    }, $class;
}

sub track_selector_and_tracks_for_species_and_assembly {
    my ($self, $species, $assembly, %opts) = @_;
    my @studies = $self->{resources}->get( $species, $assembly );
    return unless @studies;
    my @rnaseq_track_configs;
    for my $study (@studies) {
        for my $run (@{$study->{runs}}) {
          my $run_id = $run->{run_id};
          my $url    = GenomeBrowser::Deployment::sync_ebi_to_sanger(
              $species, $assembly, $run_id,
              $run->{data_files}{bigwig}, %opts );
          my $attributes = {
             %{$study->{attributes}},
             %{$run->{attributes}},
             track => join(": ", grep {$_} $run_id, $run->{run_description_short}),
             study => sprintf("%s: %s", $study->{study_id}, $study->{study_description_short}),
          };
          $attributes->{pubmed} = join(", " , pairmap {sprintf('<a href="https://www.ncbi.nlm.nih.gov/pubmed/%s">%s</a>', $a, $b->[1])} %{$study->{pubmed}}) if $study->{pubmed};
          $attributes->{study_description} = $study->{study_description_full} if $study->{study_description_full} ne $study->{study_description_short};
# We don't want both exact and approximate values to show, but we need the approximate values for facets
# So, delete exact values ( I don't know how to stop JBrowse from displaying some values) 
          delete $attributes->{library_size_reads};
          delete $attributes->{fraction_of_reads_mapping_uniquely};
          push @rnaseq_track_configs,
            {
              storeClass    => "JBrowse/Store/SeqFeature/BigWig",
              type          => "JBrowse/View/Track/Wiggle/XYPlot",
              category      => "RNASeq",
              autoscale     => "local",
              ScalePosition => "right",
              urlTemplate => $url,
              key         => "$run_id: ".$run->{run_description_full},
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
            "library_size_reads_approximate",
            "fraction_of_reads_mapping_uniquely_approximate",
            @as
        ],
        renameFacets => {
            pubmed => "PubMed",
            study                  => "Study",
            submitting_centre      => "Submitting centre",
            track  => "Track",
            library_size_reads_approximate    => "Library size (reads)",
            fraction_of_reads_mapping_uniquely_approximate  => "Fraction of reads mapping uniquely",
            %pretty
        }
    };
}

my @column_header_blacklist = (
 "synonym",
 "bioproject_id",
  "species",
  "organism",
  "replicate",
  "sample_name",
  "batch",
  "barcode",
  "insdc_center_name",
  "insdc_first_public",
  "insdc_secondary_accession",
  "insdc_status",
  "insdc_last_update",
  "label",
  "model",
  "package",
  "ncbi_submission_model",
  "ncbi_submission_package",
  "sample_comment",
  "sample_title",
  "geo_accession",
  "biological_replicate",
  "block",
  "zone", #schmidtea mediterranea
  "repplicate",
  "in_house_sample_code",
  "collected_by",
  "biomaterial_provider",
  "description_title",
  "treatment_sources",
  "population",
  "sample_name",
  "agarosemigrationtemperature",
  "agarosemigrationttime",
  "baermanntemperature",
  "base_calling_software_version",
  "culturetemperature",
  "culturetime",
  "library_id",
  "library_preparation",
  "wash",
);
sub not_in_blacklist {
  my $arg = shift;
  for (@column_header_blacklist){
     return if $arg eq $_;
  }
  return 1;
}
sub column_headers_for_studies {
  my %data;
  my $runs_total;
  my %occurrences_per_type;
  for my $study(@_){
    my %rnaseqer_characteristics;
    for my $run (@{$study->{runs}}) {
      $runs_total++;
      my %h = %{$run->{characteristics}};
      while (my ($k, $v) = each %h){
         $rnaseqer_characteristics{$k}{$v}++;
         $occurrences_per_type{$k}++;
      }
    }
    my @factors;
    for ( keys %rnaseqer_characteristics ) {
      my %d      = %{ $rnaseqer_characteristics{$_} };
      my @values = keys %d;
      my @counts = values %d;
      push @factors, $_ if @values > 1 or sum(@counts) == 1;
    }
    $data{$study->{study_id}} = \@factors;
  }
  my @result = map { @{$_} } ( values %data );
  @result = grep { not_in_blacklist($_) } @result;
  @result = uniq @result;
  @result = sort @result;
  if($runs_total > 100) {
      @result = grep {$occurrences_per_type{$_} * 20 > $runs_total} @result;
  }
  return \@result;
}


1;
