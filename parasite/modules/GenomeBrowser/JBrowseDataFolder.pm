
package GenomeBrowser::JbrowseDataFolder;
use Carp;
use File::Path qw(make_path);
use File::Slurp qw(write_file);
use JSON;
use SpeciesFtp;

# This is the data folder for jbrowse consumption
#
# input parameters: 
#  - where to construct the folder
#  - corresponding data production location


sub new {
  my ($class, $root_dir, $core_db) = @_;
  croak "Not enough args: @_" unless $root_dir and $core_db;

  my ($spe, $cies, $bioproject) = split "_", $core_db;
  my $species = join "_", $spe, $cies, $bioproject;
  my $dir = "$root_dir/$species";

  make_path $dir;

  return bless {
    species => $species,
    core_db => $core_db,
    dir => $dir,
  }, $class;
}

our $CONTENT_NAMES = {
  SEQUENCES => "seq",
  INDEXES => "names",
  TRACK_FILES_LOCAL => "tracks",
  CONFIG => "trackList.json",
  INCLUDES => "functions.conf",
}

sub path_to {
  my ($self, $name, @others) = @_;
  confess "No such JBrowse item: $name" unless $CONTENT_NAMES->{$name};
  return join "/", $self->{dir}, $CONTENT_NAMES->{$name}, @others; 
}

sub make_all {
  my ($self, $jbrowse_install, $processing_dir) = @_;

  my $jbrowse_tools= GenomeBrowser::JBrowseTools(install_location => $jbrowse_install, tmp => $processing_dir);
  unless -d $self->path_to("SEQUENCE"){
    $jbrowse_tools->prepare_sequence(
      output_path => $self->path_to("SEQUENCE"),
      input_path => SpeciesFtp->current_staging->path_to($self->{core_db}, "genomic.fa")
    );
  }
  for my $local_track (@$local_tracks){
    my $f = $self->path_to("TRACK_FILES_LOCAL", $local_track->{track_label});
    unless ( -f $f ) {
       $jbrowse_tools->track_from_annotation(
          %$local_track, 
          output_path => $f,
          input_path => SpeciesFtp->current_staging->path_to($self->{core_db}, "annotations.gff3")
       );
    }
  }
  unless -d $self->path_to("INDEXES") {
    $jbrowse_tools->index_names(output_path => $self->{dir});
  }
  unless -d $self->path_to("INCLUDES") {
    # TODO this configured colours of things.
    # Work on this much later. 
    # $self->copy_includes;
  }

  my %config = %$CONFIG_STANZA;
  my @tracks;
  push @tracks, $self->sequence_track;
  push @tracks, $self->gene_models_track;
  push @tracks, $self->feature_tracks;

  my ($attribute_query_order, @rnaseq_tracks) = GenomeBrowser::RnaseqTracks("$processing_dir/rnaseq", $self->{core_db});
  
  for my $rnaseq_track (@rnaseq_tracks) {
#TODO :)    
     deploy_rnaseq_track($rnaseq_track);
     
     push @tracks, $self->format_rnaseq_track($_);
  }
  
  
   $config{trackSelector} = {
     type => "Faceted",
     displayColumns => ["type", "category", @$attribute_query_order]
   } if @rnaseq_tracks;

  $config{tracks}=\@tracks; 

  write_file($self->path_to("CONFIG"), \%config);
}

our $CONFIG_STANZA = {

};

our $TRACK_STANZA = {
  storeClass => "JBrowse/Store/SeqFeature/BigWig",
  type => "JBrowse/View/Track/Wiggle/XYPlot",
  category => "RNASeq",
  autoscale => "local",
  ScalePosition => "right",
};

sub format_rnaseq_track {
  my ($self, %args) = @_;


}

our $local_tracks = [
  {
    'feature' => 'WormBase,WormBase_imported',
    'trackLabel' => 'Gene Models',
    'trackType' => 'CanvasFeatures',
    'category' => 'Genome Annotation',
    'type' => [qw/gene mRNA exon CDS five_prime_UTR three_prime_UTR tRNA rRNA pseudogene tRNA_pseudogene antisense_RNA lincRNA miRNA miRNA_primary_transcript mRNA piRNA pre_miRNA pseudogenic_rRNA pseudogenic_transcript pseudogenic_tRNA scRNA snoRNA snRNA ncRNA/]
  },
  {
    'feature' => 'ncrnas_predicted',
    'trackLabel' => 'Predicted non-coding RNA (ncRNA)',
    'trackType' => 'FeatureTrack',
    'category' => 'Genome Annotation',
    'type' => ['nucleotide_match']
  },
  {
    'feature' => 'RepeatMasker',
    'trackLabel' => 'Repeat Region',
    'trackType' => 'FeatureTrack',
    'category' => 'Repeat Regions',
    'type' => ['repeat_region']
  },
  {
    'feature' => 'dust',
    'trackLabel' => 'Low Complexity Region (Dust)',
    'trackType' => 'FeatureTrack',
    'category' => 'Repeat Regions',
    'type' => ['low_complexity_region']
  },
  {
    'feature' => 'tandem',
    'trackLabel' => 'Tandem Repeat (TRFs)',
    'trackType' => 'FeatureTrack',
    'category' => 'Repeat Regions',
    'type' => ['tandem_repeat']
  }
];
0;
