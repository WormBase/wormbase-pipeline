package GenomeBrowser::JBrowseTools;
use strict;
use warnings;
use JSON;
use File::Slurp qw/write_file/;
use Log::Any qw($log);

sub new {
  my ( $class, $install_location ) = @_;
  bless {install_location => $install_location}, $class;
}

sub tool_cmd {
  my ( $self, $name ) = @_;
  return "perl " . join "/", $self->{install_location}, "bin", $name;
}

sub exec_if_dir_absent {
  my ($dir, $cmd) = @_;
  return if -d $dir;
  $log->info("Executing: $cmd");
  my $output = `$cmd`;
  die $log->fatal("Failed : $cmd, output: $output") if $?;
}

sub flatfile_to_json {
  my ($self, $source_path, $out, $track_label) = @_;
  my $cmd = $self->tool_cmd("flatfile-to-json.pl");
  $cmd .= " --gff $source_path";
  $cmd .= " --trackLabel '$track_label'";
  $cmd .= ' --nameAttributes "name,id,alias,locus,target,note"';
  $cmd .= ' --compress';
  $cmd .= " --out $out";
  exec_if_dir_absent(
   "$out/tracks/$track_label", $cmd 
  );
  return $cmd;
}

sub prepare_refseqs {
  my ($self, $fasta_path, $out) = @_;
  my $cmd = $self->tool_cmd("prepare-refseqs.pl");
  $cmd .= " --fasta $fasta_path";
  $cmd .= " --compress";
  $cmd .= " --out $out";

  exec_if_dir_absent( "$out/seq", $cmd);
}

sub index_names {
  my ( $self,$out, %args ) = @_;
  my $cmd = $self->tool_cmd("generate-names.pl");
  $cmd .= " --compress";
  $cmd .= " --incremental"; # set experimentally
  $cmd .= " --mem 1024000000";    # 1GB is four times the default, 256 MB
  $cmd .= " --out $out";
  exec_if_dir_absent( "$out/names", $cmd, %args );
}

sub add_static_files {
  my ( $self, $out) = @_;
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
  write_file( join("/", $out, "functions.conf"), $doc );
}

# JBrowse leaves its own configs in trackList.json.
# They're hard to manipulate programmatically so we overwrite them
# and add canned copies.
# You can JSON->new->pretty if you like pretty, but it will make the files bigger.
sub update_config {
  my ( $self, $out, $new_config ) = @_;
  my $path = join("/", $out, "trackList.json");
  $log->info("Updating config: $path");
  write_file( $path, JSON->new->utf8->encode($new_config) );
}
1;
