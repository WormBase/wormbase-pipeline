package GenomeBrowser::JBrowseTools;

use File::Basename;
use File::Copy qw(copy);
use IO::Uncompress::Gunzip qw(gunzip);

sub new {
    my ( $class, %args ) = @_;
    die unless "$args{install_location}";
    die unless "$args{tmp}";
    bless \%args, $class;
}

sub tool_cmd {
    my ( $self, $name ) = shift;
    return "perl " . join "/", $self->{install_location}, "bin", $name;
}

sub track_from_annotation {
    my ( $self, %args ) = @_;

    my $path = $args{input_path};
    my $name = basename $path;

    if ( $path =~ /.gz$/ ) {
        $name =~ s/.gz$//;
        my $tmp = join "/", $self->{tmp}, $name;
        $tmp =~ s/.gz$//;
        gunzip( $path, $tmp );
        $path = $tmp;
    }
    my $processing_path = join "/", $self->{tmp}, ( $name . $args{feature} );
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
    $cmd .= " --type " . join ",", $args{type};
    $cmd .= " --key \"$args{trackLabel}\"";
    $cmd .= " --trackLabel $track_label";
    $cmd .= " --trackType $args{trackType}";
    $cmd .= ' --nameAttributes "name,id,alias,locus"';
    $cmd .= ' --compress';
    $cmd .= " --out $args{output_path}";

#TODO
#--metadata '{ "category": "%s", "menuTemplate" : [{ "label" : "View gene at WormBase ParaSite", "action" : "newWindow", "url" : "/Gene/Summary?g={name}" }] }'
#$gff3_track->{'type'} =~ /^gene/ ? qq(--clientConfig '{ "color" : "{geneColor}", "label" : "{geneLabel}" }') : '',

    `$cmd`;
}

sub prepare_sequence {
    my ( $self, %args ) = @_;

    my $path = $args{input_path};
    my $name = basename $path;

    if ( $path =~ /.gz$/ ) {
        $name =~ s/.gz$//;
        my $tmp = join "/", $self->{tmp}, $name;
        $tmp =~ s/.gz$//;
        gunzip( $path, $tmp );
        $path = $tmp;
    }
    my $cmd = $self->tool_cmd("flatfile-to-json.pl");
    $cmd .= " --fasta $path";
    $cmd .= " --compress";
    $cmd .= " --out $args{output_path}";
    `$cmd`;
}

sub index_names {
  my ($self, %args ) = @_;
  my $cmd = $self->tool_cmd("generate-names.pl");
  $cmd .= " --compress";
  $cmd .= " --out $args{output_path}";
}
1;
