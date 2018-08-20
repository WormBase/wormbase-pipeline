
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
my $help_menu_html= <<EOF;
<div class="jbrowse help_dialog">
<div class="main">
<dl>
<dt>Welcome to the secret tips menu!</dt>
<dd><ul>
Here are a couple of less discoverable features you might want to try if you haven't used them yet.
<br>
For an overview of usual JBrowse features please refer to the <a href="http://jbrowse.org">JBrowse main page</a>.
<br>
</ul></dd>
<dt>Advanced navigation</dt>
<dd><ul>
<li>You can use the search box to jump between regions and scaffolds: <br>
e.g. inputting <i>chr4:79,500,000..80,000,000</i> jumps to the region on chromosome 4 between 79.5Mb and 80Mb</li>
<li>You can also search for gene identifiers there:<br>
e.g. <i>WBGene0001324</i></li>
<li>Right clicking on a gene model lets you see the corresponding gene page in WormBase ParaSite </li>
</ul></dd>
<dt>"About this track"</dt>
This option in the track label menu is full of useful information about the RNASeq tracks:
<dd><ul>
<li>All metadata available for the track</li>
<li>Links to ENA study page and other resources where available</li>
<li>Link to FTP location with data files: QC, quantification, and more</li>
</ul></dd>
<dt>Track selector</dt>
There are frequently too many tracks to view at once. <br>
Select runs in the study you need, or pick several runs with good library size across several studies with good metadata. <br>
Searching for keywords may help you narrow down the results.</br>
<br>
If you have any questions, suggestions, or comments, write to us at <a href="mailto:parasite-help\@sanger.ac.uk">parasite-help\@sanger.ac.uk</a>!
</dl>
</div>
</div>
EOF
$help_menu_html =~ s/\n//g;
my $CONFIG_STANZA = {
    "names" => {
        "type" => "Hash",
        "url"  => "names/"
    },
    "include" =>
      [ 
        "functions.conf"
      ],
  "css" => ".detail .field_container .field.track {display: none} .detail .field_container .value_container.track {display: none}",
  "aboutThisBrowser" =>  {
      "title"=> "WormBase ParaSite:",
      "description" => "<div class=\"default_about\"><img src='/i/parasite.png' /><br>Browser for WormBase ParaSite genomes and annotation features <br>RNASeq tracks have been imported from <a href=\"https://www.ebi.ac.uk/fg/rnaseq/api/\">the RNASeq-er project</a><br></div>",
  },
  "quickHelp" => {
     "title"=> "Usage tips", 
     "content" => $help_menu_html
  }
};

my $TRACK_STANZA = {
    storeClass    => "JBrowse/Store/SeqFeature/BigWig",
    type          => "JBrowse/View/Track/Wiggle/XYPlot",
    category      => "RNASeq",
    autoscale     => "local",
    ScalePosition => "right",
};
my $local_tracks_metadata_stanza = {
           "study"=> "(WormBase track)",
        "submitting_centre"=> "WormBase",
        "fraction_of_reads_mapping_uniquely_approximate"=> "(not applicable)",
        "library_size_reads_approximate"=> "(not applicable)"
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
    'metadata'    => {
      'category' => 'Reference sequence',
      'track' => 'Reference sequence',
      %$local_tracks_metadata_stanza
    }
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
            'label'  => 'View gene in WormBase ParaSite'
        }
    ],
    'metadata' => {
        'category' => 'Genome Annotation',
        'track' => 'Gene Models',
        %$local_tracks_metadata_stanza
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
        'track' => $niceLabelForReading,
        %$local_tracks_metadata_stanza
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
  my ($self,  %opts ) = @_;
  for my $core_db (ProductionMysql->staging->core_databases){
      next unless $core_db =~ /_core_$ENV{PARASITE_VERSION}/;
      print "Starting: $core_db\n";
      $self->make_tracks($core_db, %opts);
  }
}

sub make_tracks {
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
          my $attributes = {
             %{$study->{attributes}},
             %{$run->{attributes}},
             track => $run->{run_description_short}
          };
# We don't want both exact and approximate values to show, but we need the approximate values for facets
# So, delete exact values ( I don't know how to stop JBrowse from displaying some values) 
          delete $attributes->{library_size_reads};
          delete $attributes->{fraction_of_reads_mapping_uniquely};
          push @rnaseq_track_configs,
            {
              %$TRACK_STANZA,
              urlTemplate => $url,
              key         => $run->{run_description_full},
              label       => "RNASeq/$run_id",
              metadata    => $attributes,
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
        displayColumns   => [ "track", @as ],
        selectableFacets => [
            "category",            "study",
            "submitting_centre",
            "library_size_reads_approximate",
            "fraction_of_reads_mapping_uniquely_approximate",
            @as
        ],
        renameFacets => {
            study                  => "Study",
            submitting_centre      => "Submitting centre",
            track  => "Track",
            library_size_reads_approximate    => "Library size (reads)",
            fraction_of_reads_mapping_uniquely_approximate  => "Fraction of reads mapping uniquely",
            %pretty
        }
    };
}
1;
