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
    for my $core_db ( ProductionMysql->staging->core_databases($opts{core_dbs_pattern}) ) {
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
    my ( $attribute_query_order, $location_per_run_id, @studies ) =
      $self->{resources}->get( $core_db, $assembly );
    @studies = sort {$a->{study_id} cmp $b->{study_id}} @studies;
    my $runs_path = $self->path("$species.runs.tsv");
    unlink $runs_path if -f $runs_path;
    
    return unless @studies;
    open(my $runs_fh, '>:utf8', $runs_path ) or die $runs_path;

    print $runs_fh join ("\t",
       "Study", "Track",
       "Library size (reads)", "Mapping quality (reads uniquely mapped)", 
       "Results", 
       map {( my $a = $_) =~ s/[\W_-]+/ /g; ucfirst($a)} @$attribute_query_order
    ) . "\n";
    for my $study (@studies) {
        my $study_id = $study->{study_id};
        for my $run ( sort {$a->{run_id} cmp $b->{run_id}} @{ $study->{runs} } ) {
            my $run_id = $run->{run_id};
            print $runs_fh join ("\t",
                 $study_id, $run->{run_description_short},
                 $run->{attributes}{library_size_approx}, $run->{attributes}{mapping_quality_approx},
                 dirname($location_per_run_id->{$run_id}),
                 map {$run->{attributes}{$_} || '' } @$attribute_query_order
            ). "\n";     
        }
    }
    close $runs_fh;
    
    my $studies_path = $self->path("$species.studies.tsv");
    unlink $studies_path if -f $studies_path;

    open(my $studies_fh, '>:utf8', $studies_path) or die $studies_path;
    print $studies_fh join ("\t",
      "Study", "Submitting centre", "ENA first public", "ENA last update",
      "Description",
      "PubMed id", "PubMed",
    )."\n";
    for my $study (@studies){
      print $studies_fh join ("\t",
        $study->{study_id},
        $study->{attributes}{submitting_centre}, $study->{attributes}{"ENA first public"}, $study->{attributes}{"ENA last update"},
        $study->{study_description_short},
        join(", " , keys ($study->{pubmed}|| {})), join(", " , map {$_->[1]} values ($study->{pubmed}|| {})),
      )."\n";
    }
    close $studies_fh; 
}
1;
