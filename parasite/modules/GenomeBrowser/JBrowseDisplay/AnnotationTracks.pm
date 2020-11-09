use strict;
use warnings;
package GenomeBrowser::JBrowseDisplay::AnnotationTracks;

use SpeciesFtp;
use GenomeBrowser::JBrowseTools;
use File::Path qw(make_path remove_tree);
use IO::Uncompress::Gunzip qw(gunzip);
use File::Basename;
use List::Util qw/pairs/;
use Log::Any qw/$log/;

sub new {
  my ( $class, %args ) = @_;
  die @_
    unless $args{jbrowse_tools}
    and $args{tmp_dir}
    and $args{species_ftp};
  bless \%args, $class;
}

sub tracks_for_species {
  my ($self, $species, $out, %opts) = @_;
  my $tmp_dir = join ("/", $self->{tmp_dir}, $species);
  make_path $tmp_dir;
  my $annotation_path_ftp = $self->{species_ftp}->path_to( $species, "annotations.gff3" );
  my $annotation_path = join ("/", $tmp_dir, basename $annotation_path_ftp);
  $annotation_path =~ s/.gz$//;
  unless (-f $annotation_path) {
    remove_tree $_ for glob("$tmp_dir/*");
    $log->info (__PACKAGE__ . " unzipping $annotation_path_ftp -> $annotation_path.tmp");
    gunzip( $annotation_path_ftp, "$annotation_path.tmp" );
    split_annotation_into_files_by_track_label($tmp_dir, "$annotation_path.tmp");
    rename "$annotation_path.tmp", $annotation_path;
  }
  my %gff_paths_by_track_label = map {(basename $_) => $_ } grep {$_ !~ /annotations.gff3/ } glob("$tmp_dir/*");
  for my $track_label (sort keys %gff_paths_by_track_label){
    $self->{jbrowse_tools}->flatfile_to_json($gff_paths_by_track_label{$track_label}, $out, $track_label, %opts);
  }
  return map {config_for_track_label($_)} keys %gff_paths_by_track_label;
}
sub split_annotation_into_files_by_track_label {
  my ($tmp_dir, $annotation_path) = @_;
  my %out_fhs;
  my %out_paths;
  open (my $annotation_fh, "<", $annotation_path) or die "$!: $annotation_path";
  while(<$annotation_fh>){
    next if /^#/;
    my @F = split "\t";
    my $track_label = track_label_for_source_and_type(@F[1,2]);
    next unless $track_label;
    unless ($out_fhs{$track_label}){
      $log->info($track_label);
      $out_paths{$track_label} = join("/", $tmp_dir, $track_label);
      open ($out_fhs{$track_label}, ">", $out_paths{$track_label}) or die "$!: $out_paths{$track_label} $_";
    }
    print {$out_fhs{$track_label}} $_;
  }
  close $_ for values %out_fhs;
}
sub track_label_for_source_and_type {
  my ($source, $type) = @_;
  if (grep {$type eq $_} qw/gene mRNA exon CDS five_prime_UTR three_prime_UTR tRNA rRNA pseudogene tRNA_pseudogene antisense_RNA lincRNA miRNA miRNA_primary_transcript nc_primary_transcript mRNA piRNA pre_miRNA pseudogenic_rRNA pseudogenic_transcript transposable_element pseudogenic_tRNA scRNA snoRNA snRNA ncRNA nontranslating_CDS/){
    if ($source eq "WormBase" || $source eq "WormBase_imported"){
      return "Gene_Models";
    } elsif ($source eq "WormBase_transposon"){ 
      return "Transposons";
    }else {
      return "";
    }
  }
  return "" if $source eq "WormBase" and $type eq "intron";
  return "" if $source eq "history";
  return "" if $type eq "assembly_component";
  
  return "$source.$type";
}

my $local_tracks_metadata_stanza = {
  "study"=> "(WormBase track)",
  "submitting_centre"=> "WormBase",
  "fraction_of_reads_mapping_uniquely_approximate"=> "(not applicable)",
  "library_size_reads_approximate"=> "(not applicable)"
};
my $genes_track_config = {
  'style' => {
    'className' => 'feature',
    'color'     => '{geneColor}',
    'label'     => '{geneLabel}'
  },
  'key'          => 'Gene Models',
  'storeClass'   => 'JBrowse/Store/SeqFeature/NCList',
  'trackType'    => 'CanvasFeatures',
  'urlTemplate'  => 'tracks/Gene_Models/{refseq}/trackData.jsonz',
  'compress'     => 1,
  'menuTemplate' => [
    {
      'url'    => '/Gene/Summary?g={name}',
      'action' => 'newWindow',
      'label'  => 'View gene in WormBase ParaSite'
    }
  ],
  'metadata' => {
    'category' => 'Genome annotation',
    'track' => 'Gene Models',
    %$local_tracks_metadata_stanza
  },
  'label' => 'Gene_Models'
};
sub feature_track_config {
  my ($category, $track_label, $display_label, $style) = @_;

  return {
    'style' => $style,
    'key'         => $display_label,
    'storeClass'  => 'JBrowse/Store/SeqFeature/NCList',
    'trackType'   => 'CanvasFeatures',
    'urlTemplate' => "tracks/$track_label/{refseq}/trackData.jsonz",
    'compress'    => 1,
    'metadata'    => {
      'category' => $category,
      'track' => $display_label,
      %$local_tracks_metadata_stanza
    },
    'label' => $track_label,
  };
}
my %repeat_labels_by_source = (
  dust => "Low complexity region (Dust)",
  repeatmasker => "Repetitive region",
  RepeatMasker =>  "Repetitive region",
  trf => "Tandem repeat",
  tandem => "Tandem repeat",
);
sub config_for_track_label {
  my ($track_label) = @_;
  my ($source, $type) = split (qr/\./, $track_label, 2);
  $type //="";

  return (
   $track_label eq "Gene_Models" 
     ? $genes_track_config 
     : $repeat_labels_by_source{$source} 
        ? feature_track_config(
            "Repeat regions", $track_label, $repeat_labels_by_source{$source},
            { 
              strandArrow => "",
              color => $source eq "dust" ? "cornflowerblue" : "steelblue",
              showLabels => ($source ne "dust"),
            },
        )
        : $source =~ /rfam|cmscan/i
          ? feature_track_config(
              "Genome annotation", $track_label , "Predicted non-coding RNAs", {
                color => "saddlebrown",
              }
          )
          : $type eq "expressed_sequence_match" || $type eq "protein_match"
            ? feature_track_config(
                "Sequence similarity", $track_label , (do {$source =~ s/_/ /g; $source}) , {
                color => ($source =~ /BEST/ ? "indianred" : $source =~ /OTHER/ ? "lightcoral" : "hotpink"),
              }
            )
            : feature_track_config(
                "Other features", $track_label , (do {
                   my $name = $type =~ /$source/i ? ucfirst $type : $type ? ucfirst "$source: $type": ucfirst $source;
                   $name =~ s/_/ /g;
                   $name 
                }) , {}
              )
  );
}
1;
