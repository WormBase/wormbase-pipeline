# Produces a folder of TSV files
# Format:
# study run location  num_reads_approx  mapping_quality_approx  @factors
package GenomeBrowser::Tsv;

use ProductionMysql;
use GenomeBrowser::Resources;
use File::Path qw(make_path);
use File::Slurp qw(write_file);
use File::Basename qw(dirname);

sub new {
    my ( $class, %args ) = @_;
    $args{root_dir} //=
      "$ENV{PARASITE_SCRATCH}/jbrowse/WBPS$ENV{PARASITE_VERSION}";
    return bless {
        dir       => "$args{root_dir}/tsv",
        resources => GenomeBrowser::Resources->new("$args{root_dir}/Resources"),
    }, $class;
}

sub path {
    my ( $self, @args ) = @_;
    my $result = join "/", $self->{dir}, @args;
    print "$result\n" if $ENV{TSV_VERBOSE};
    return $result;
}

sub make_all {
    my ( $self, %opts ) = @_;

    make_path $self->path;

    my @genomes;
    for my $core_db ( ProductionMysql->staging->core_databases ) {
        $self->make_tsv( $core_db, %opts );
    }
}

sub make_tsv {
    my ( $self, $core_db, %opts ) = @_;

    my ( $spe, $cies, $bioproject ) = split "_", $core_db;
    return if $bioproject eq "core";
    my $species = join "_", $spe, $cies, $bioproject;
    my $assembly =
      ProductionMysql->staging->meta_value( $core_db, "assembly.name" );
    my $path = $self->path("$species.tsv");
    unlink $path if -f $path;
    my ( $attribute_query_order, $location_per_run_id, @studies ) =
      $self->{resources}->get( $core_db, $assembly );
    
    return unless @studies;
    open(my $fh, '>:utf8', $path ) or die $core_db;

    print $fh join ("\t", 
       "Study", "Submitting centre", "Track",
       "Library size (reads)", "Mapping quality (reads uniquely mapped)", 
       "Results", 
       map {( my $a = $_) =~ s/[\W_-]+/ /g; ucfirst($a)} @$attribute_query_order
    ) . "\n";
    for my $study (@studies) {
        my $study_id = $study->{study_id};
        for my $run ( @{ $study->{runs} } ) {
            my $run_id = $run->{run_id};
            print $fh join ("\t",
                 $study_id, $study->{attributes}{submitting_centre}, $run->{run_description},
                 $run->{attributes}{library_size_approx}, $run->{attributes}{mapping_quality_approx},
                 dirname($location_per_run_id->{$run_id}),
                 map {$run->{attributes}{$_} || '' } @$attribute_query_order
            ). "\n";     
        }
    }
    close $fh;
}
1;
