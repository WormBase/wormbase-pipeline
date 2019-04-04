package GenomeBrowser::JBrowseTools;
use strict;
use warnings;
use File::Path qw(make_path);
use File::Basename;
use File::Copy qw(copy);
use IO::Uncompress::Gunzip qw(gunzip);
use JSON;
use File::Slurp;
use Log::Any qw($log);

sub new {
  my ( $class, %args ) = @_;
  die @_
    unless $args{install_location}
    and $args{tmp_dir}
    and $args{out_dir}
    and $args{species_ftp};
  bless \%args, $class;
}

sub tool_cmd {
  my ( $self, $name ) = @_;
  return "perl " . join "/", $self->{install_location}, "bin", $name;
}

sub tmp_path {
  my ( $self, $core_db, $name ) = @_;

  $core_db =~ s/_core.*//;
  $name =~ s/.gz//;
  my $dir = join "/", $self->{tmp_dir}, $core_db;
  make_path $dir;
  return join "/", $dir, $name;
}

sub out_dir {
  my ( $self, $core_db ) = @_;
  $core_db =~ s/_core.*//;
  return join "/", $self->{out_dir}, $core_db;
}

sub exec_if_dir_absent {
  my ( $self, $dir, $cmd, %args ) = @_;
  if ( $args{do_jbrowse} // not -d $dir ) {
    $log->info("Executing: $cmd");
    my $output = `$cmd`;
    die $log->fatal("Failed : $cmd, output: $output") if $?;
  }
  else {
    $log->info("Skipping: $cmd");
  }
}

sub filter_gff {
  my ( $self, $path, $processing_path, $features, %args ) = @_;

  if ( $args{do_jbrowse} // not -f $processing_path ) {
    $log->info("filter_gff processing: $processing_path");
    open ANNOTATION, "$path" or die "Could not open: $path";
    open PROCESSING, ">$processing_path"
      or die "Could not write: $processing_path";
    while (<ANNOTATION>) {
      my @gff_elements = split( "\t", $_ );
      print PROCESSING $_ if ( $gff_elements[1] ~~ $features );
    }
    close ANNOTATION;
    close PROCESSING;
  }
  else {
    $log->info("filter_gff skipping: $processing_path");
  }
}

sub track_from_annotation {
  my ( $self, %args ) = @_;

  my $out = $self->out_dir( $args{core_db} );
  ( my $track_label = $args{'trackLabel'} ) =~ s/\s/_/g;
  my $path =
    $self->{species_ftp}->path_to( $args{core_db}, "annotations.gff3" );

  if ( $path =~ /.gz$/ ) {
    my $tmp = $self->tmp_path( $args{core_db}, basename $path);
    gunzip( $path, $tmp ) unless -f $tmp;
    $path = $tmp;
  }
  my $processing_path =
    $self->tmp_path( $args{core_db}, ( join ",", @{ $args{feature} } ) );

  $self->filter_gff( $path, $processing_path, $args{feature}, %args );

  if ( -s $processing_path ) {
    my $cmd = $self->tool_cmd("flatfile-to-json.pl");
    $cmd .= " --gff $processing_path";
    $cmd .= " --type " . join ",", @{ $args{type} };
    $cmd .= " --key '$args{trackLabel}'";
    $cmd .= " --trackLabel '$track_label'";
    $cmd .= " --trackType $args{trackType}";
    $cmd .= ' --nameAttributes "name,id,alias,locus"';
    $cmd .= ' --compress';
    $cmd .= " --out $out";
    $self->exec_if_dir_absent( "$out/tracks/$track_label", $cmd, %args );
  }
  else {
    $log->info("Skipping flatfile-to-json.pl:  $processing_path");
  }
}

sub track_present {
  my ( $self, %args ) = @_;

  my $out = $self->out_dir( $args{core_db} );
  ( my $track_label = $args{'trackLabel'} ) =~ s/\s/_/g;

  return -d "$out/tracks/$track_label";
}

sub prepare_sequence {
  my ( $self, %args ) = @_;

  my $out = $self->out_dir( $args{core_db} );
  my $path = $self->{species_ftp}->path_to( $args{core_db}, "genomic.fa" );

  if ( $path =~ /.gz$/ ) {
    my $tmp = $self->tmp_path( $args{core_db}, basename $path);
    die "Sequence file missing: $path" unless -s $path;
    gunzip( $path, $tmp ) unless -s $tmp;
    $path = $tmp;
  }
  my $cmd = $self->tool_cmd("prepare-refseqs.pl");
  $cmd .= " --fasta $path";
  $cmd .= " --compress";
  $cmd .= " --out $out";

  $self->exec_if_dir_absent( "$out/seq", $cmd, %args );
}

sub index_names {
  my ( $self, %args ) = @_;
  my $out = $self->out_dir( $args{core_db} );
  my $cmd = $self->tool_cmd("generate-names.pl");
  $cmd .= " --compress";
  $cmd .= " --mem 1024000000";    # 1GB is four times the default, 256 MB
  $cmd .= " --out $out";
  $self->exec_if_dir_absent( "$out/names", $cmd, %args );
}

sub add_static_files {
  my ( $self, %args ) = @_;
  my $out = $self->out_dir( $args{core_db} );
  my $doc = <<END_FUNCTIONS_CONF;
geneLabel = function(f) {
  var type = f.get('type');
  var locus = f.get('locus');
  var seq_name = f.get('sequence_name');
  var feature_name = f.get('Name');
  var patt = /RNA|transcript/;
  if(patt.test(type)) { return feature_name; }
  if(typeof seq_name !== 'undefined') {
    if(typeof locus !== 'undefined') {
      return locus + " (" + seq_name + ")";
    } else {
      return seq_name;
    }
  } else {
    if(typeof locus !== 'undefined') {
      return locus;
    } else {
      return feature_name;
    }
  }}

geneColor = function(f) {
  var type = f.get('type');
  if (type.match(/exon/)) {return 'pink';}
  if (type.match(/pseudo/)) {return 'pink';}
  var strand = f.get('strand');
  if (strand == -1) {return 'turquoise';}
  if (strand ==  1) {return 'violet';}
  return 'gray'; }
END_FUNCTIONS_CONF
  make_path $out;
  write_file( "$out/functions.conf", $doc );
}

# JBrowse leaves its own configs in trackList.json.
# They're hard to manipulate programmatically so we overwrite them
# and add canned copies.
# You can JSON->new->pretty if you like pretty, but it will make the files bigger.
sub update_config {
  my ( $self, %args ) = @_;
  my $config_path = join "/", $self->out_dir( $args{core_db} ),
    "trackList.json";
  my $new_config = $args{new_config};
  write_file( $config_path, JSON->new->utf8->encode($new_config) );
}
1;
