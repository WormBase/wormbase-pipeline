
package GenomeBrowser::JBrowseDisplay;
use strict;
use Carp;
use File::Path qw(make_path);
use File::Slurp qw(write_file);
use JSON;
use SpeciesFtp;
use GenomeBrowser::JBrowseTools;
use GenomeBrowser::Resources;
use GenomeBrowser::Deployment;
use ProductionMysql;

# This is the data folder for jbrowse consumption
#
# input parameters:
#  - where to construct the folder
#  - corresponding data production location
sub new {
    my ( $class, %args ) = @_;

    $args{jbrowse_install} //=
"/nfs/production/panda/ensemblgenomes/wormbase/software/packages/jbrowse/JBrowse-1.12.5";
    $args{root_dir} //=
      "$ENV{PARASITE_SCRATCH}/jbrowse/WBPS$ENV{PARASITE_VERSION}";

    make_path "$args{root_dir}/out";
    make_path "$args{root_dir}/JBrowseTools";
    make_path "$args{root_dir}/Resources";
    return bless {
        dir           => "$args{root_dir}/out",
        jbrowse_tools => GenomeBrowser::JBrowseTools->new(
            install_location => $args{jbrowse_install},
            tmp_dir          => "$args{root_dir}/JBrowseTools",
            out_dir          => "$args{root_dir}/out",
            species_ftp      => $args{ftp_path}
            ? SpeciesFtp->new( $args{ftp_path} )
            : SpeciesFtp->current_staging,
        ),
        resources =>
          GenomeBrowser::Resources->new("$args{root_dir}/Resources"),
    }, $class;
}

my $CONFIG_STANZA = {
    "names" => {
        "type" => "Hash",
        "url"  => "names/"
    },
    "include" =>
      [ #Gives us the nice gene labels. TODO there's no code to copy them now I think!
        "functions.conf"
      ]
};

my $TRACK_STANZA = {
    storeClass    => "JBrowse/Store/SeqFeature/BigWig",
    type          => "JBrowse/View/Track/Wiggle/XYPlot",
    category      => "RNASeq",
    autoscale     => "local",
    ScalePosition => "right",
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
    'category'    => 'Reference sequence'
};

my $genes_track_config = {
    'style' => {
        'className' => 'feature',
        'color'     => '{geneColor}',
        'label'     => '{geneLabel}'
    },
    'key'          => 'Gene Models',
    'storeClass'   => 'JBrowse/Store/SeqFeature/NCList',
    'trackType'    => 'CanvasFeatures',
    'urlTemplate'  => 'tracks/Gene_Models/{refseq}/trackData.jsonz',
    'compress'     => 1,
    'menuTemplate' => [
        {
            'url'    => '/Gene/Summary?g={name}',
            'action' => 'newWindow',
            'label'  => 'View gene at WormBase ParaSite'
        }
    ],
    'metadata' => {
        'category' => 'Genome Annotation',
    },
    'type'  => 'CanvasFeatures',
    'label' => 'Gene_Models'
};

sub feature_track_config {
    my $niceLabelForReading = shift;
    (my $label = $niceLabelForReading ) =~ s/ /_/g;

    return {
        'style' => {
            'className' => 'feature'
        },
        'key'         => $niceLabelForReading,
        'storeClass'  => 'JBrowse/Store/SeqFeature/NCList',
        'trackType'   => 'FeatureTrack',
        'urlTemplate' => "tracks/$label/{refseq}/trackData.jsonz",
        'compress'    => 1,
        'metadata'    => {
            'category' => 'Repeat Regions',    #All our features are repeats now
        },
        'type'  => 'FeatureTrack',
        'label' => $label,
    };
}
my $genes_track = {
    'feature'    => [qw/WormBase WormBase_imported/],
    'trackLabel' => 'Gene Models',
    'trackType'  => 'CanvasFeatures',
    'type'       => [
        qw/gene mRNA exon CDS five_prime_UTR three_prime_UTR tRNA rRNA pseudogene tRNA_pseudogene antisense_RNA lincRNA miRNA miRNA_primary_transcript mRNA piRNA pre_miRNA pseudogenic_rRNA pseudogenic_transcript pseudogenic_tRNA scRNA snoRNA snRNA ncRNA/
    ]
};

my $feature_tracks = [
    {
        'feature'    => ['ncrnas_predicted'],
        'trackLabel' => 'Predicted non-coding RNA (ncRNA)',
        'trackType'  => 'FeatureTrack',
        'type'       => ['nucleotide_match']
    },
    {
        'feature'    => ['RepeatMasker'],
        'trackLabel' => 'Repeat Region',
        'trackType'  => 'FeatureTrack',
        'type'       => ['repeat_region']
    },
    {
        'feature'    => ['dust'],
        'trackLabel' => 'Low Complexity Region (Dust)',
        'trackType'  => 'FeatureTrack',
        'type'       => ['low_complexity_region']
    },
    {
        'feature'    => ['tandem'],
        'trackLabel' => 'Tandem Repeat (TRFs)',
        'trackType'  => 'FeatureTrack',
        'type'       => ['tandem_repeat']
    }
];

sub make_all {
    my ( $self, $core_db, %opts ) = @_;

    my ( $spe, $cies, $bioproject ) = split "_", $core_db;

    my $species = join "_", $spe, $cies, $bioproject;

    $self->{jbrowse_tools}->prepare_sequence(
        core_db => $core_db,
        %opts
    );
    $self->{jbrowse_tools}->track_from_annotation(
        %$genes_track,
        core_db => $core_db,
        %opts
    );

    my @feature_track_configs;
    for my $feature_track (@$feature_tracks) {
        my $args = {
            %$feature_track,
            core_db => $core_db,
            %opts
        };
        $self->{jbrowse_tools}->track_from_annotation(%$args);

        push @feature_track_configs,
          feature_track_config( $feature_track->{trackLabel} )
          if $self->{jbrowse_tools}->track_present(%$args);
    }

    $self->{jbrowse_tools}->index_names( core_db => $core_db, %opts );
    $self->{jbrowse_tools}->add_static_files( core_db => $core_db, %opts );

    my $assembly =
      ProductionMysql->staging->meta_value( $core_db, "assembly.name" );
    my ( $attribute_query_order, $location_per_run_id, @studies ) =
      $self->{resources}->get( $core_db, $assembly );
    my @rnaseq_track_configs;
    for my $study (@studies) {
        for my $run (@{$study->{runs}}) {
          my $run_id = $run->{run_id};
          my $url    = GenomeBrowser::Deployment::sync_ebi_to_sanger( $run_id,
              $location_per_run_id->{$run_id}, %opts );
          push @rnaseq_track_configs,
            {
              %$TRACK_STANZA,
              urlTemplate => $url,
              key         => $run->{run_description},
              label       => "RNASeq/$run_id",
              metadata    => {%{$study->{attributes}}, %{$run->{attributes}}}
            };
        }
    }

    my %config = %$CONFIG_STANZA;
    if (@rnaseq_track_configs) {
        $config{trackSelector} = $self->track_selector(@$attribute_query_order);
        $config{defaultTracks} = "DNA,Gene_Models";
    }
    else {
        #Default track selector
        #All local tracks on
        $config{defaultTracks} = join ",", "DNA",
          map { my $m = $_->{'trackLabel'}; $m =~ s/\s/_/g; $m }
          @$feature_tracks;
    }

    $config{tracks} = [
        $sequence_track_config, $genes_track_config,
        @feature_track_configs, @rnaseq_track_configs
    ];
    $config{containerID} =
      "WBPS$ENV{PARASITE_VERSION}_${species}_" . scalar( @{ $config{tracks} } );
    return $self->{jbrowse_tools}
      ->update_config( core_db => $core_db, new_config => \%config );
}

sub track_selector {
    my ( $self, @as ) = @_;
    my %pretty;
    for my $a (@as) {
        ( my $p = $a ) =~ s/[\W_-]+/ /g;
        $pretty{$a} = ucfirst($p);
    }
    return {
        type             => "Faceted",
        displayColumns   => [ "key", @as ],
        selectableFacets => [
            "category",            "study",
            "submitting_centre",
            "library_size_approx", "mapping_quality_approx",
            @as
        ],
        renameFacets => {
            study                  => "Study",
            submitting_centre      => "Submitting centre",
            key                    => "Track",
            library_size_approx    => "Library size (reads)",
            mapping_quality_approx => "Mapping quality (reads uniquely mapped)",
            %pretty
        }
    };
}
1;
