use strict;
use warnings;
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
package GenomeBrowser::RnaseqHub;

use ProductionMysql;
use GenomeBrowser::Deployment;
use List::MoreUtils qw/uniq/;
use File::Path qw(make_path);
use File::Slurp qw(write_file read_file);
use File::Basename qw(dirname);
use Log::Any qw($log);
use JSON;

sub new {
    my ( $class, $root_dir) = @_;
    return bless {
       out_dir       => "$root_dir/hub",
       wbps_expression_dir => "$root_dir/WbpsExpression",
    }, $class;
}

sub path {
    my ( $self, @args ) = @_;
    my $result = join "/", $self->{out_dir}, @args;
    $log->info(__PACKAGE__ . ": $result");
    return $result;
}

sub make_hubs {
    my ( $self, $core_dbs_pattern, %opts ) = @_;

    make_path $self->path;

    for my $core_db ( ProductionMysql->staging->core_databases($core_dbs_pattern) ) {
        $self->make_hub_for_core_db($core_db, %opts);
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

sub make_hub_for_core_db {
    my ( $self, $core_db, %opts ) = @_;
    my ( $spe, $cies, $bioproject ) = split "_", $core_db;
    my $species = join "_", $spe, $cies, $bioproject;
    my $_Species = ucfirst($species);
    my $path = join("/", $self->{wbps_expression_dir}, $species, "${spe}_${cies}.studies.json");
    return unless -f $path;
    my @studies = @{from_json(read_file($path, { binmode => ':utf8' }))};
    return unless @studies;

    my $assembly =
      ProductionMysql->staging->meta_value( $core_db, "assembly.default" );
    make_path( $self->path( $_Species, "doc" ) );
    my @study_tracks;
    my @run_tracks;
    my $has_multiple_categories = 1 < uniq map {$_->{study_category}} @studies;
    for my $study (sort {$a->{study_category} cmp $b->{study_category} || $a->{study_id} cmp $b->{study_id} } @studies) {
        my $study_id = $study->{study_id};
        &create_study_doc( $self->path( $_Species, "doc", "$study_id.html" ),
            $study );
        push @study_tracks,
          &study_track( $study_id, $study->{study_title},  ($has_multiple_categories ? "$study->{study_category} - $study_id": $study_id));
        for my $run ( @{ $study->{runs} } ) {
            my $run_id = $run->{run_id};
            &create_run_doc( $self->path( $_Species, "doc", "$run_id.html" ),
                $study, $run, );
            my $url = GenomeBrowser::Deployment::sync_ebi_to_sanger(
                $species, $assembly, $run_id,
                $run->{bigwig}, %opts );
            push @run_tracks,
              &run_track( $study_id, $run_id, $run->{condition}, $url );
        }
    }

    &create_trackDb( $self->path( $_Species, "trackDb.txt" ),
        @study_tracks, @run_tracks );

    write_file($self->path($_Species, "genomes.txt"), "genome $assembly\ntrackDb trackDb.txt\n");
    &create_hub_file( $self->path($_Species, "hub.txt") , $species);
}
#Maybe: same format as JBrowse?
sub study_track {
    my ( $study_id, $study_description_full, $short_label ) = @_;
    return (
        "track $study_id
superTrack on
group $study_id
shortLabel $short_label
longLabel $study_description_full
html doc/$study_id
"
    );
}

sub run_track {
    my ( $study_id, $run_id, $condition, $url ) = @_;
    my $track_name = join(": ", grep {$_} $run_id, $condition);
    return (
        "track $run_id
parent $study_id
type bigWig
bigDataUrl $url
shortLabel $track_name
longLabel $track_name
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

#TODO I'm not sure where this gets displayed
sub create_study_doc {
    my ( $path, $study ) = @_;
    my $header = $study->{attributes}{"Study description"} // $study->{study_title};
    my %attributes = %{$study->{attributes}};
    delete $attributes{"Study description"};
    write_file( $path, {binmode => ':utf8'}, attributes_html($header, \%attributes));
}
sub attributes_html {
  my ($header, $attributes) = @_;  

    my $result = "<table>\n";
    $result .= "<th><b>$header</b></th>\n";
    for my $k ( sort keys %{$attributes} ) {
        ( my $property_name = $k ) =~ s/_/ /g;
        $property_name = ucfirst($property_name) if $property_name eq lc($property_name);
        my $property_value = $attributes->{$k};
        $result .= sprintf "<tr>\n  <td>%s</td>\n  <td>%s</td>\n</tr>\n", $property_name, $property_value;
    }
    $result .= "</table>\n";
    return $result;
}
sub create_run_doc {
    my ( $path, $study, $run ) = @_;

    my $track_name = join(": ", grep {$_} $run->{run_id}, $run->{condition});
    write_file( $path, {binmode => ':utf8'}, attributes_html("Run $track_name", $run->{attributes}) );
}
1;
