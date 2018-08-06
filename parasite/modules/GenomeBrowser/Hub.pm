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

sub create_hub {
    my ( $self, %opts ) = @_;

    make_path $self->path;
    my @core_dbs_with_tracks;

    my @genomes;
    for my $core_db ( ProductionMysql->staging->core_databases ) {
        my $genome = $self->create_for_species( $core_db, %opts );
        push @genomes, $genome if $genome;
    }

    &create_genomes_file( $self->path("genomes.txt"), @genomes );
    &create_hub_file( $self->path("hub.txt") );
}

sub create_hub_file {
    my ($path) = @_;
    write_file(
        $path,
        'hub WBPS-RNASeq
shortLabel RNA-Seq Alignments
longLabel RNA-Seq Alignments for WormBase ParaSite
genomesFile genomes.txt
email parasite-help@sanger.ac.uk
'
    );
}

sub create_genomes_file {
    my ( $path, @genomes ) = @_;
    open( my $fh, '>', $path ) or die $path;
    for my $genome (@genomes) {
        my $g = $genome->{genome};
        my $t = $genome->{trackDb};
        print $fh "genome $g\ntrackDb $t\n\n";
    }
    close $fh;
}

sub create_for_species {
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
              &run_track( $study_id, $run_id, $run->{run_description}, $url );
        }
    }

    &create_trackDb( $self->path( $_Species, "trackDb.txt" ),
        @study_tracks, @run_tracks );
    return {
        genome  => $assembly,
        trackDb => "$_Species/trackDb.txt"
    };
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
    my ( $study_id, $run_id, $run_description, $url ) = @_;
    return (
        "track $run_id
parent $study_id
type bigWig
bigDataUrl $url
shortLabel $run_id
longLabel $run_description
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
    return sprintf "Study %s: %s\n", $study->{study_id},
      $study->{study_description};
}

#TODO I'm not sure where this gets displayed
sub create_study_doc {
    my ( $path, $study ) = @_;
    my $result = &_study_intro($study);
    write_file( $path, {binmode => ':utf8'}, $result );
}

sub create_run_doc {
    my ( $path, $study, $run ) = @_;

    my $result = &_study_intro($study);
    $result = "<h4>$result</h4>";
    $result .= sprintf "<h5>Run %s</h5>\n", $run->{run_description};
    $result .= "\n";
    for my $k ( sort keys $run->{attributes} ) {
        ( my $property_name = $k ) =~ s/_/ /g;
        $property_name = ucfirst($property_name)
          if $property_name eq lc($property_name);
        $result .= sprintf "<b>%s</b>\t%s\n", $property_name,
          $run->{attributes}{$k};
    }
    $result =~ s/\n/<br>\n/g;
    write_file( $path, {binmode => ':utf8'}, $result );
}
1;
