use Test::More;

use File::Temp qw/tempdir/;
use GenomeBrowser::JBrowseTools;
use SpeciesFtp;
use JSON;

my $canned_json = do {local $/; <DATA>};
my $subject = GenomeBrowser::JBrowseTools->new(
  install_location => "/invalid/JBrowse_install_path",
  tmp_dir => "/invalid/tmp",
  out_dir => "/invalid/out",
  species_ftp => SpeciesFtp->new("/invalid/ftp")
);

is_deeply($subject->merge_configs({},{}), {}, "null case");
is_deeply($subject->merge_configs(from_json($canned_json),{}), from_json($canned_json), "merge empty hash");

my $track = {
  "key" => "test track",
  "urlTemplate"=>"./rnaseq"
};

my $other_track = {
  "key" => "other test track",
  "urlTemplate"=>"./rnaseq"
};
my $out = $subject->merge_configs(from_json($canned_json),{tracks=>[$track]});
ok(scalar(grep {$_->{key} =~/test/ } @{$out->{tracks}}) == 1, "New track gets added");
my @entries = grep {$_->{key} !~/test/} @{$out->{tracks}};
is_deeply(\@entries, from_json($canned_json)->{tracks}, "Old tracks unchanged");

my $out_2 = $subject->merge_configs($out, {tracks=>[$track, $other_track]});
ok(scalar(grep {$_->{key} =~/test/ } @{$out_2->{tracks}}) == 2, "New track gets added 2");

my $out_3 = $subject->merge_configs($out_2, {tracks=>[$other_track]});
ok(scalar(grep {$_->{key} =~/test/ } @{$out_3->{tracks}}) == 1, "New track gets removed");

done_testing();
__DATA__
{
  "tracks": [
    {
      "key": "Reference sequence",
      "storeClass": "JBrowse/Store/Sequence/StaticChunked",
      "chunkSize": 80000,
      "urlTemplate": "seq/{refseq_dirpath}/{refseq}-",
      "compress": 1,
      "label": "DNA",
      "type": "SequenceTrack",
      "category": "Reference sequence"
    },
    {
      "style": {
        "className": "feature",
        "color": "{geneColor}",
        "label": "{geneLabel}"
      },
      "key": "Gene Models",
      "storeClass": "JBrowse/Store/SeqFeature/NCList",
      "trackType": "CanvasFeatures",
      "urlTemplate": "tracks/Gene_Models/{refseq}/trackData.jsonz",
      "compress": 1,
      "metadata": {
        "category": "Genome Annotation",
        "menuTemplate": [
          {
            "action": "newWindow",
            "url": "/Gene/Summary?g={name}",
            "label": "View gene at WormBase ParaSite"
          }
        ]
      },
      "type": "CanvasFeatures",
      "label": "Gene_Models"
    },
    {
      "style": {
        "className": "feature"
      },
      "key": "Predicted non-coding RNA (ncRNA)",
      "storeClass": "JBrowse/Store/SeqFeature/NCList",
      "trackType": "FeatureTrack",
      "urlTemplate": "tracks/Predicted_non-coding_RNA_(ncRNA)/{refseq}/trackData.jsonz",
      "compress": 1,
      "type": "FeatureTrack",
      "label": "Predicted_non-coding_RNA_(ncRNA)",
      "metadata": {
        "category": "Genome Annotation",
        "menuTemplate": [
          {
            "url": "/Gene/Summary?g={name}",
            "action": "newWindow",
            "label": "View gene at WormBase ParaSite"
          }
        ]
      }
    },
    {
      "style": {
        "className": "feature"
      },
      "key": "Repeat Region",
      "storeClass": "JBrowse/Store/SeqFeature/NCList",
      "trackType": "FeatureTrack",
      "urlTemplate": "tracks/Repeat_Region/{refseq}/trackData.jsonz",
      "compress": 1,
      "label": "Repeat_Region",
      "metadata": {
        "category": "Repeat Regions",
        "menuTemplate": [
          {
            "action": "newWindow",
            "url": "/Gene/Summary?g={name}",
            "label": "View gene at WormBase ParaSite"
          }
        ]
      },
      "type": "FeatureTrack"
    },
    {
      "style": {
        "className": "feature"
      },
      "key": "Low Complexity Region (Dust)",
      "storeClass": "JBrowse/Store/SeqFeature/NCList",
      "trackType": "FeatureTrack",
      "urlTemplate": "tracks/Low_Complexity_Region_(Dust)/{refseq}/trackData.jsonz",
      "compress": 1,
      "metadata": {
        "category": "Repeat Regions",
        "menuTemplate": [
          {
            "url": "/Gene/Summary?g={name}",
            "action": "newWindow",
            "label": "View gene at WormBase ParaSite"
          }
        ]
      },
      "type": "FeatureTrack",
      "label": "Low_Complexity_Region_(Dust)"
    },
    {
      "style": {
        "className": "feature"
      },
      "key": "Tandem Repeat (TRFs)",
      "storeClass": "JBrowse/Store/SeqFeature/NCList",
      "trackType": "FeatureTrack",
      "urlTemplate": "tracks/Tandem_Repeat_(TRFs)/{refseq}/trackData.jsonz",
      "compress": 1,
      "type": "FeatureTrack",
      "label": "Tandem_Repeat_(TRFs)",
      "metadata": {
        "category": "Repeat Regions",
        "menuTemplate": [
          {
            "action": "newWindow",
            "url": "/Gene/Summary?g={name}",
            "label": "View gene at WormBase ParaSite"
          }
        ]
      }
    }
  ],
  "names": {
    "url": "names/",
    "type": "Hash"
  },
  "include": [
    "functions.conf"
  ]
}
