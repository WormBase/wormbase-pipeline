package GenomeBrowser::JBrowseTools;
use strict;
use File::Path qw(make_path);
use File::Basename;
use File::Copy qw(copy);
use IO::Uncompress::Gunzip qw(gunzip);
use JSON;
use File::Slurp;
use Hash::Merge;
sub new {
    my ( $class, %args ) = @_;
    die @_ unless $args{install_location} and $args{tmp_dir} and $args{out_dir} and $args{species_ftp};
    bless \%args, $class;
}

sub tool_cmd {
    my ( $self, $name ) = @_;
    return "perl " . join "/", $self->{install_location}, "bin", $name;
}
sub tmp_path {
  my ($self, $core_db, $name) = @_;

  $core_db =~ s/_core.*//;
  $name =~s/.gz//;
  my $dir = join "/", $self->{tmp_dir}, $core_db;
  make_path $dir;
  return join "/", $dir, $name;
}
sub out_dir {
  my ($self, $core_db) = @_;
  $core_db =~ s/_core.*//;
  return join "/", $self->{out_dir}, $core_db;
}
sub track_from_annotation {
    my ( $self, %args ) = @_;

    my $path = $self->{species_ftp}->path_to($args{core_db}, "annotations.gff3") ;

    if ( $path =~ /.gz$/ ) {
        my $tmp = $self->tmp_path($args{core_db}, basename $path);
        gunzip( $path, $tmp );
        $path = $tmp;
    }

    my $processing_path = $self->tmp_path($args{core_db}, (basename $path . $args{feature}));
    open ANNOTATION, "$path" or die "Could not open: $path";
    open PROCESSING, ">$processing_path"
      or die "Could not write: $processing_path";
    while (<ANNOTATION>) {
        my @gff_elements = split( "\t", $_ );
        print PROCESSING $_ if ( $gff_elements[1] ~~ $args{type} );
    }
    close ANNOTATION;
    close PROCESSING;

    ( my $track_label = $args{'trackLabel'} ) =~ s/\s/_/g;
    my $cmd = $self->tool_cmd("flatfile-to-json.pl");
    $cmd .= " --gff $processing_path";
    $cmd .= " --type " . join ",", @{ $args{type} };
    $cmd .= " --key '$args{trackLabel}'";
    $cmd .= " --trackLabel '$track_label'";
    $cmd .= " --trackType $args{trackType}";
    $cmd .= ' --nameAttributes "name,id,alias,locus"';
    $cmd .= ' --compress';
    $cmd .= " --out ".$self->out_dir($args{core_db});

#TODO
#--metadata '{ "category": "%s", "menuTemplate" : [{ "label" : "View gene at WormBase ParaSite", "action" : "newWindow", "url" : "/Gene/Summary?g={name}" }] }'
#$gff3_track->{'type'} =~ /^gene/ ? qq(--clientConfig '{ "color" : "{geneColor}", "label" : "{geneLabel}" }') : '',

    print STDERR "Executing: $cmd\n" if $self->{verbose};
    `$cmd`;
    $? and die "Failed: $cmd\n";
}

sub prepare_sequence {
    my ( $self, %args ) = @_;

    my $path = $self->{species_ftp}->path_to($args{core_db}, "genomic.fa") ;

    if ( $path =~ /.gz$/ ) {
        my $tmp = $self->tmp_path($args{core_db}, basename $path);
        gunzip( $path, $tmp );
        $path = $tmp;
    }
    my $cmd = $self->tool_cmd("prepare-refseqs.pl");
    $cmd .= " --fasta $path";
    $cmd .= " --compress";
    $cmd .= " --out ".$self->out_dir($args{core_db});
    print STDERR "Executing: $cmd\n" if $self->{verbose};
    `$cmd`;
    $? and die "Failed: $cmd\n";
}

sub index_names {
    my ( $self, %args ) = @_;
    my $cmd = $self->tool_cmd("generate-names.pl");
    $cmd .= " --compress";
    $cmd .= " --mem 1024000000";         # 1GB is four times the default, 256 MB
    $cmd .= " --out ".$self->out_dir($args{core_db});
    print STDERR "Executing: $cmd\n" if $self->{verbose};
    `$cmd`;
    $? and die "Failed: $cmd\n";
}

sub update_config {
   my ( $self, %args ) = @_;
   my $config_path = join "/", $self->out_dir($args{core_db}), "trackList.json";
   write_file($config_path, $self->merge_configs(from_json(read_file($config_path)), $args{new_config}));
}
sub merge_configs {
  my ($self, $current_config, $new_config) = @_;

# remove rnaseq tracks from the config
  if ($current_config->{tracks}) {
    my @tracks = @{$current_config->{tracks}}; 
    @tracks = grep {$_->{urlTemplate} !~ /rnaseq/i } @tracks;
    $current_config->{tracks} = \@tracks;
  }

  return Hash::Merge::merge($current_config, $new_config);
}
1;
