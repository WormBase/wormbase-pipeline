use strict;
use warnings;
package GenomeBrowser::JBrowseDisplay::SequenceTrack;

use SpeciesFtp;
use GenomeBrowser::JBrowseTools;
use IO::Uncompress::Gunzip qw(gunzip);
use File::Basename;

sub new {
  my ( $class, %args ) = @_;
  die @_
    unless $args{jbrowse_tools}
    and $args{tmp_dir}
    and $args{species_ftp};
  bless \%args, $class;
}
my $local_tracks_metadata_stanza = {
           "study"=> "(WormBase track)",
        "submitting_centre"=> "WormBase",
        "fraction_of_reads_mapping_uniquely_approximate"=> "(not applicable)",
        "library_size_reads_approximate"=> "(not applicable)"
};
my $sequence_track_config = {
    'seqType'     => 'dna',
    'key'         => 'Reference sequence',
    'chunkSize'   => 80000,
    'storeClass'  => 'JBrowse/Store/Sequence/StaticChunked',
    'urlTemplate' => 'seq/{refseq_dirpath}/{refseq}-',
    'compress'    => 1,
    'label'       => 'DNA',
    'type'        => 'SequenceTrack',
    'metadata'    => {
      'category' => 'Reference sequence',
      'track' => 'Reference sequence',
      %$local_tracks_metadata_stanza
    }
};

sub track_for_species {
  my ($self, $species, $out, %opts) = @_;
  
  my $fasta_path = $self->{species_ftp}->path_to( $species, "genomic.fa" );
  if (-f $fasta_path && $fasta_path =~ /.gz$/){
      my $tmp_path = join ("/", $self->{tmp_dir}, $species, basename $fasta_path);
      $tmp_path =~ s/.gz$//;
      gunzip( $fasta_path, $tmp_path ) unless -f $tmp_path;
      $fasta_path = $tmp_path;
  };
  $self->{jbrowse_tools}->prepare_refseqs($fasta_path, $out, %opts);

  return $sequence_track_config;
}
1;
