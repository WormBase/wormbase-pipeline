
package GenomeBrowser::Tracks;
use strict;
use Carp;
use File::Path qw(make_path);
use File::Slurp qw(write_file);
use JSON;
use SpeciesFtp;
use GenomeBrowser::JBrowseTools;
use GenomeBrowser::RnaseqTracks;
use GenomeBrowser::Deployment;
use ProductionMysql;

# This is the data folder for jbrowse consumption
#
# input parameters: 
#  - where to construct the folder
#  - corresponding data production location


sub new {
  my ($class, %args) = @_;
  
  $args{jbrowse_install} //= "/nfs/production/panda/ensemblgenomes/wormbase/software/packages/jbrowse/JBrowse-1.12.5";
  $args{root_dir} //= "$ENV{PARASITE_SCRATCH}/jbrowse/WBPS$ENV{PARASITE_VERSION}";

  make_path "$args{root_dir}/out";
  make_path "$args{root_dir}/JBrowseTools";
  make_path "$args{root_dir}/RnaseqTracks";
  return bless {
    dir => "$args{root_dir}/out",
    jbrowse_tools => GenomeBrowser::JBrowseTools->new(
       install_location => $args{jbrowse_install},
       tmp_dir => "$args{root_dir}/JBrowseTools",
       out_dir => "$args{root_dir}/out",
       species_ftp =>  $args{ftp_path} ? SpeciesFtp->new($args{ftp_path}) : SpeciesFtp->current_staging, 
    ),
    rnaseq_tracks => GenomeBrowser::RnaseqTracks->new("$args{root_dir}/RnaseqTracks"),
  }, $class;
}

# bugs: this is gonna be merged with itself on repeated runs
my $CONFIG_STANZA = {
   "names" => {
      "type" => "Hash",
      "url" => "names/"
   },
};

my $TRACK_STANZA = {
  storeClass => "JBrowse/Store/SeqFeature/BigWig",
  type => "JBrowse/View/Track/Wiggle/XYPlot",
  category => "RNASeq",
  autoscale => "local",
  ScalePosition => "right",
};


my $local_tracks = [
  {
    'feature' => [qw/WormBase WormBase_imported/],
    'trackLabel' => 'Gene Models',
    'trackType' => 'CanvasFeatures',
    'category' => 'Genome Annotation',
    'type' => [qw/gene mRNA exon CDS five_prime_UTR three_prime_UTR tRNA rRNA pseudogene tRNA_pseudogene antisense_RNA lincRNA miRNA miRNA_primary_transcript mRNA piRNA pre_miRNA pseudogenic_rRNA pseudogenic_transcript pseudogenic_tRNA scRNA snoRNA snRNA ncRNA/]
  },
  {
    'feature' => ['ncrnas_predicted'],
    'trackLabel' => 'Predicted non-coding RNA (ncRNA)',
    'trackType' => 'FeatureTrack',
    'category' => 'Genome Annotation',
    'type' => ['nucleotide_match']
  },
  {
    'feature' => ['RepeatMasker'],
    'trackLabel' => 'Repeat Region',
    'trackType' => 'FeatureTrack',
    'category' => 'Repeat Regions',
    'type' => ['repeat_region']
  },
  {
    'feature' => ['dust'],
    'trackLabel' => 'Low Complexity Region (Dust)',
    'trackType' => 'FeatureTrack',
    'category' => 'Repeat Regions',
    'type' => ['low_complexity_region']
  },
  {
    'feature' => ['tandem'],
    'trackLabel' => 'Tandem Repeat (TRFs)',
    'trackType' => 'FeatureTrack',
    'category' => 'Repeat Regions',
    'type' => ['tandem_repeat']
  }
];
sub make_all {
  my ($self,$core_db, %opts) = @_;
 
  my ($spe, $cies, $bioproject) = split "_", $core_db;

  my $species = join "_", $spe, $cies, $bioproject;
   
  my @track_configs;

  $self->{jbrowse_tools}->prepare_sequence(
    core_db => $core_db,
  );
  #push @track_configs, GenomeBrowser::TrackConfig::sequence_track();

  for my $local_track (@$local_tracks){
    $self->{jbrowse_tools}->track_from_annotation(
        %$local_track, 
        core_db => $core_db,
    );

  }
  $self->{jbrowse_tools}->index_names(core_db=>$core_db);
  

  my $assembly = ProductionMysql->staging->meta_value($core_db, "assembly.name");  
  my ($attribute_query_order, $location_per_run_id, @rnaseq_tracks) = $self->{rnaseq_tracks}->get($core_db, $assembly);
  for my $rnaseq_track (@rnaseq_tracks) {
     my $run_id = $rnaseq_track->{run_id};
     my $url = GenomeBrowser::Deployment::sync_ebi_to_sanger($run_id, $location_per_run_id->{$run_id}, %opts);
     push @track_configs, {
       %$TRACK_STANZA,
       urlTemplate => $url,
       key => $rnaseq_track->{label},
       label => "RNASeq/$run_id",
       metadata => $rnaseq_track->{attributes}
     };
  }
  
  
  my %config = %$CONFIG_STANZA;
   $config{trackSelector} =
     $self->track_selector(@$attribute_query_order)
     if @rnaseq_tracks;

  $config{tracks}=\@track_configs; 
  return $self->{jbrowse_tools}->update_config(core_db=>$core_db, new_config=> \%config);
}

sub track_selector {
  my ($self, @as) = @_;
  my %pretty;
  my @ks;
  for my $a (@as){
    (my $k = $a) =~ s/\W+/-/;
    $k = lc($k);
    push @ks, $k;
    $pretty{$k}=$a;
  }
  return {
    type => "Faceted",
    displayColumns => \@ks,
    selectableFacets => [ "category", @ks],
    renameFacets => \%pretty
  }
}
1;
