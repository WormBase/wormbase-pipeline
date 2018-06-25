
package GenomeBrowser::Tracks;

use ProductionMysql;
use GenomeBrowser::ArrayExpressMetadata;
use GenomeBrowser::RnaseqerMetadata;
use GenomeBrowser::Studies;
use GenomeBrowser::Factors;
# This is the convenient representation of tracks
# Fetches its dependent data and puts them on disk

# you can edit the dependent data by hand, like configs

# Convert me to JBrowse format
# Dump me to yaml to see what I'm like
# Use me to move files EBI -> Sanger
# Create hub.txt
sub new {
  my ($class, $root_dir, $core_db) = @_;
  my $assembly = ProductionMysql->staging->meta_value($core_db, "assembly.name");
  my ($spe, $cies, $bioproject) = split "_", $core_db;
  my $rnaseqer_metadata = GenomeBrowser::RnaseqerMetadata->new($root_dir, "${spe}_${cies}");
  my $array_express_metadata = GenomeBrowser::ArrayExpressMetadata->new($root_dir, "${spe}_${cies}");
  my $studies = GenomeBrowser::Studies->new($root_dir, "${spe}_${cies}", $assembly, $rnaseqer_metadata); 
  my $factors = GenomeBrowser::Factors->new($root_dir, "${spe}_${cies}", $assembly, $rnaseqer_metadata, $array_express_metadata);
  my %data;
  for my $study_id (@{$rnaseqer_metadata->access($assembly)}){
    my $study_attributes = $studies->{$study_id};
    for my $run_id (@{$rnaseqer_metadata->access($assembly, $study_id)}){
       my $run_attributes = $rnaseqer_metadata->{$assembly}{$study_id}{$run_id};
       my %attributes = %{{ %$study_attributes, %$run_attributes }};
       $attributes{label} = &_label($run_id,
          map {exists $attributes{$_} ? $attributes{$_}: ()} @$factors
       );
       $attributes{urlTemplate} = "http://ngs.sanger.ac.uk/production/parasites/wormbase/RNASeq_alignments/${spe}_${cies}_${bioproject}/$run_id.bw";
       $data{$run_id}=\%attributes;
    }
  }
  return bless {tracks => \%data }, $class;
}

sub _label {
  my ($run_id, @factor_values) = @_;
  return $run_id unless @factor_values;
  return "$run_id: ". join (", ",
    grep !/[A-Z]+[0-9]+/, @factor_values
  ); 
}

our $TRACK_STANZA = {
  storeClass => "JBrowse/Store/SeqFeature/BigWig",
  type => "JBrowse/View/Track/Wiggle/XYPlot",
  category => "RNASeq",
  autoscale => "local",
  ScalePosition => "right",
};

sub as_jbrowse_track_list {
  my $self = shift;

  my @result;
  for my $track_id (keys $self->{tracks}){
     push @result, {%{$self->{tracks}{$track_id}}, %$TRACK_STANZA};
  }
  sort { scalar (%$b) <=> scalar (%$a) || $a->{label} cmp $b->{label} } @result;
  return { tracks => \@result };
}
sub list_names {
   return keys shift;
}

sub source_path {
  my ($self, $track_id) = @_;
  return "ebi:/arrayexpress_folder/$track_id.bw";
}

sub storage_path {
  my ($self, $track_id) = @_;
  return "sanger:/ngs_folder/$track_id.bw";
}

sub _public_url {
  my ($self, $track_id) = @_;
  return $self->{$track_id}{urlTemplate};
}
1;
