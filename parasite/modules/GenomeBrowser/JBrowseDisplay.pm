
package GenomeBrowser::JBrowseDisplay;
use strict;
use warnings;
use File::Path qw(make_path);
use SpeciesFtp;
use ProductionMysql;
use Log::Any qw/$log/;

use GenomeBrowser::Deployment;
use GenomeBrowser::JBrowseTools;
use GenomeBrowser::JBrowseDisplay::AnnotationTracks;
use GenomeBrowser::JBrowseDisplay::SequenceTrack;
use GenomeBrowser::JBrowseDisplay::RnaseqTracks;

# This is the data folder for jbrowse consumption
#
# input parameters:
#  - where to construct the folder
#  - corresponding data production location
sub new {
    my ( $class, %args ) = @_;
    die "$class: need root dir" unless $args{root_dir};
    $args{jbrowse_install} //=
"/nfs/production/panda/ensemblgenomes/wormbase/software/packages/jbrowse/JBrowse-1.12.5";

    make_path "$args{root_dir}/out";
    make_path "$args{root_dir}/JBrowseTools";
    make_path "$args{root_dir}/Resources";
    my $species_ftp = $args{ftp_path} ? SpeciesFtp->new( $args{ftp_path} ) : SpeciesFtp->dot_next;
    die "No FTP dir at " . $species_ftp->{root} unless $species_ftp->root_exists;
    my $jbrowse_tools = GenomeBrowser::JBrowseTools->new( $args{jbrowse_install});
    return bless {
        dir           => "$args{root_dir}/out",
        jbrowse_tools => $jbrowse_tools,
        sequence_track =>
          GenomeBrowser::JBrowseDisplay::SequenceTrack->new(
           jbrowse_tools => $jbrowse_tools,
           species_ftp => $species_ftp,
           tmp_dir => "$args{root_dir}/SequenceTrack",
          ),
        annotation_tracks => 
          GenomeBrowser::JBrowseDisplay::AnnotationTracks->new(
           jbrowse_tools => $jbrowse_tools,
           species_ftp => $species_ftp,
           tmp_dir => "$args{root_dir}/AnnotationTracks",
        ),
        rnaseq_tracks => 
          GenomeBrowser::JBrowseDisplay::RnaseqTracks->new("$args{root_dir}/WbpsExpression"),
    }, $class;
}

sub make_displays {
  my ($self, $core_db_pattern,  %opts ) = @_;
  my @core_dbs = ProductionMysql->staging->core_databases($core_db_pattern);
  die "No core dbs for: $core_db_pattern" unless @core_dbs;
  for my $core_db (@core_dbs){
      $self->make_display_for_core_db($core_db, %opts);
  }
}

sub make_display_for_core_db {
    my ( $self, $core_db, %opts ) = @_;
    die unless $core_db =~ /_core_/;
    return unless $core_db =~ /_core_$ENV{PARASITE_VERSION}/;
    $log->info("make_tracks: $core_db");

    my ( $spe, $cies, $bioproject ) = split "_", $core_db;

    my $species = join "_", $spe, $cies, $bioproject;
    my $out = join ("/", $self->{dir}, $species);
    make_path $out;
    my $sequence_track_config = $self->{sequence_track}->track_for_species(
       $species, $out, %opts
    );
    my @annotation_track_configs = $self->{annotation_tracks}->tracks_for_species(
       $species, $out, %opts
    );

    $self->{jbrowse_tools}->index_names( $out, %opts );
    $self->{jbrowse_tools}->add_static_files( $out, %opts );

    my ($track_selector, @rnaseq_track_configs) = $self->{rnaseq_tracks}->track_selector_and_tracks_for_species_and_assembly(
      $species, ProductionMysql->staging->meta_value( $core_db, "assembly.name" ), %opts
    );

    my %config = %{config_stanza()};
    if (@rnaseq_track_configs) {
        $config{trackSelector} = $track_selector;
        $config{defaultTracks} = "DNA,Gene_Models";
    }
    else {
        #Default track selector
        #All local tracks on
        $config{defaultTracks} = join ",", "DNA",
          map { $_->{label} }
          @annotation_track_configs;
    }

    $config{tracks} = [
        $sequence_track_config, @annotation_track_configs, @rnaseq_track_configs
    ];
    $config{containerID} =
      "WBPS$ENV{PARASITE_VERSION}_${species}_" . scalar( @{ $config{tracks} } );
    return $self->{jbrowse_tools}
      ->update_config( $out, \%config );
}
sub config_stanza {
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
  return {
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
}
1;
