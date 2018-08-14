# Creates hub with this structure:
# hub/genomes.txt
# hub/hub.txt
# hub/rnaseq/
# hub/rnaseq/<species>
# hub/rnaseq/<species>/doc
# hub/rnaseq/<species>/doc/<run>.html
# hub/rnaseq/<species>/doc/<study>.html
# hub/rnaseq/<species>/trackDb.txt

# SuperTracks for studies
# BigWig tracks for runs
# Page with the attributes, for each
package GenomeBrowser::Hub;

use ProductionMysql;
use GenomeBrowser::Resources;
use GenomeBrowser::Deployment;
use File::Path qw(make_path);
use File::Slurp qw(write_file);
use File::Basename qw(dirname);

sub new {

    #in: hub directory
    #params to make RNASeqTracks
    my ( $class, %args ) = @_;
    $args{root_dir} //=
      "$ENV{PARASITE_SCRATCH}/jbrowse/WBPS$ENV{PARASITE_VERSION}";
    return bless {
        dir       => "$args{root_dir}/hub",
        resources => GenomeBrowser::Resources->new("$args{root_dir}/Resources"),
    }, $class;
}

sub path {
    my ( $self, @args ) = @_;
    my $result = join "/", $self->{dir}, @args;
    print "$result\n" if $ENV{HUB_VERBOSE};
    return $result;
}
sub make_all {
    my ( $self, %opts ) = @_;

    make_path $self->path;

    for my $core_db ( ProductionMysql->staging->core_databases ) {
        $self->make_hub($core_db, %opts);
    }
}

sub create_hub_file {
    my ($path, $species) = @_;
    my ($spe, $cies) = split "_", $species;
    my $species_nice= ucfirst($spe)." ".$cies;
    write_file(
        $path,
        "hub $species-RNASeq
shortLabel RNA-Seq Alignments
longLabel $species_nice RNA-Seq Alignments for WormBase ParaSite
genomesFile genomes.txt
email parasite-help\@sanger.ac.uk
"
    );
}

sub make_hub {
    my ( $self, $core_db, %opts ) = @_;

    my ( $spe, $cies, $bioproject ) = split "_", $core_db;
    return if $bioproject eq "core";
    my $species = join "_", $spe, $cies, $bioproject;
    my $_Species = ucfirst($species);
    my $assembly =
      ProductionMysql->staging->meta_value( $core_db, "assembly.name" );
    my ( $attribute_query_order, $location_per_run_id, @studies ) =
      $self->{resources}->get( $core_db, $assembly );
    return unless @studies;
    make_path( $self->path( $_Species, "doc" ) );
    my @study_tracks;
    my @run_tracks;

    for my $study (@studies) {
        my $study_id = $study->{study_id};
        &create_study_doc( $self->path( $_Species, "doc", "$study_id.html" ),
            $study );
        push @study_tracks,
          &study_track( $study_id, $study->{study_description} );
        for my $run ( @{ $study->{runs} } ) {
            my $run_id = $run->{run_id};
            &create_run_doc( $self->path( $_Species, "doc", "$run_id.html" ),
                $study, $run, );
            my $url = GenomeBrowser::Deployment::sync_ebi_to_sanger( $run_id,
                $location_per_run_id->{$run_id}, %opts );
            push @run_tracks,
              &run_track( $study_id, $run_id, $run->{run_description_short}, $run->{run_description_full}, $url );
        }
    }

    &create_trackDb( $self->path( $_Species, "trackDb.txt" ),
        @study_tracks, @run_tracks );

    write_file($self->path($_Species, "genomes.txt"), "genome $assembly\ntrackDb trackDb.txt\n");
    &create_hub_file( $self->path($_Species, "hub.txt") , $species);
}

sub study_track {
    my ( $study_id, $study_description ) = @_;
    return (
        "track $study_id
superTrack on
group $study_id
shortLabel $study_id
longLabel $study_description
html doc/$study_id
"
    );
}

sub run_track {
    my ( $study_id, $run_id, $run_description_short, $run_description_full, $url ) = @_;
    return (
        "track $run_id
parent $study_id
type bigWig
bigDataUrl $url
shortLabel $run_description_short
longLabel $run_description_full
color 0,0,0
html doc/$run_id
visibility hide
"
    );
}

sub create_trackDb {
    my ( $path, @tracks ) = @_;

    write_file( $path,{binmode => ':utf8'}, join( "\n", @tracks ) );
}

sub _study_intro {
    my $study = shift;
    my $si = $study->{study_id};
    my $sd = $study->{study_description};
    return $sd ? "Study $si: $sd\n" : "Study $si\n"; 
}

#TODO I'm not sure where this gets displayed
sub create_study_doc {
    my ( $path, $study ) = @_;
    my $result = &_study_intro($study);
    write_file( $path, {binmode => ':utf8'}, $result );
}

sub create_run_doc {
    my ( $path, $study, $run ) = @_;

    my $rd = $run->{run_description_full};
    my $result = "<table><th>$d</th>";
    $result .= sprintf "<b>Run $rd</b>\n" if $rd; 
    $result .= "\n";
    for my $k ( sort keys $run->{attributes} ) {
        next if grep {$_ eq $k} qw/library_size_reads_approximate fraction_of_reads_mapping_uniquely_approximate sample_name/;
        ( my $property_name = $k ) =~ s/_/ /g;
        $property_name = ucfirst($property_name)
          if $property_name eq lc($property_name);
        $property_name = "Library size (reads)" if $k eq "library_size_reads";
        my $property_value = $run->{attributes}{$k};
        $result .= sprintf "<tr><td>%s</td><td>%s</td></tr>", $property_name, $property_value;
    }
    $result .= "</table>";
    write_file( $path, {binmode => ':utf8'}, $result );
}
1;
