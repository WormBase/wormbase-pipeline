use Test::More;
use Test::MockModule;

use File::Temp qw/tempdir/;
use JSON;
use GenomeBrowser::RnaseqerMetadata;

my $module = new Test::MockModule('GenomeBrowser::LocallyCachedResource');
our $json = from_json(do {local $/; <DATA>});
$module->mock('get_json', sub {
  my ($class, $url) = @_;
  if ($url =~ /getRunsByOrganism\/schistosoma_mansoni/){
    return $json->[0]; #cropped payload to one small study
  } elsif($url=~/getSampleAttributesPerRunByStudy\/SRP026308/) {
     return $json->[1]; #sample attributes for this study
  } else {
    return undef;
  }
});

my $root_dir = tempdir( CLEANUP => 1 );
my $species = "schistosoma_mansoni";

my $subject = GenomeBrowser::RnaseqerMetadata->new($root_dir, $species);

is_deeply(
  $subject->{location_per_run_id},
  {
    SRR922067 => "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR922/SRR922067/SRR922067.nospliced.bw",
    SRR922068=> "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR922/SRR922068/SRR922068.nospliced.bw",
  },
  "Get location from RNASeq-er"
);
is_deeply(
  $subject->access,
  ["ASM23792v2"],
  "Assemblies",
) or diag explain $subject;
is_deeply(
   $subject->access("ASM23792v2"), 
   ["SRP026308"], 
   "Studies for assembly"
) or diag explain $subject;
is_deeply(
  $subject->access("wrong assembly"),
  [],
  "No studies for wrong assembly"
);

is_deeply(
   $subject->access("ASM23792v2","SRP026308"),
   ["SRR922067", "SRR922068"],
   "Runs for assembly and study"
) or diag explain $subject;

is_deeply(
   $subject->access("ASM23792v2","wrong study"), 
   [], 
   "No runs for assembly and wrong study"
) or diag explain $subject;
is_deeply(
   $subject->access("wrong assembly","SRP026308"), 
   [], 
   "No runs for wrong assembly and study"
) or diag explain $subject;

is_deeply(
   $subject->access("ASM23792v2","SRP026308", "SRR922067"), 
   ['developmental_stage','organism_part','sample_name','source_name','strain'],
   "Property types for assembly, study, and run. Normalise values. Skip dummy type with bad value of 'not applicable'"
) or diag explain $subject;
is_deeply(
   $subject->access("ASM23792v2","SRP026308", "SRR922067","strain"), 
   "NMRI",
   "Single property value for assembly, study, run, and study"
) or diag explain $subject;
done_testing();
__DATA__
[
  [
    {
      "STUDY_ID": "SRP026308",
      "SAMPLE_IDS": "SAMN02213661",
      "BIOREP_ID": "SRR922068",
      "RUN_IDS": "SRR922068",
      "ORGANISM": "schistosoma_mansoni",
      "REFERENCE_ORGANISM": "schistosoma_mansoni",
      "STATUS": "Complete",
      "ASSEMBLY_USED": "ASM23792v2",
      "ENA_LAST_UPDATED": "Tue Apr 15 2014 07:18:23",
      "LAST_PROCESSED_DATE": "Mon Dec 25 2017 09:25:33",
      "CRAM_LOCATION": "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR922/SRR922068/SRR922068.cram",
      "BEDGRAPH_LOCATION": "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR922/SRR922068/SRR922068.bedgraph",
      "BIGWIG_LOCATION": "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR922/SRR922068/SRR922068.bw",
      "MAPPING_QUALITY": 95
    },
    {
      "STUDY_ID": "SRP026308",
      "SAMPLE_IDS": "SAMN02213660",
      "BIOREP_ID": "SRR922067",
      "RUN_IDS": "SRR922067",
      "ORGANISM": "schistosoma_mansoni",
      "REFERENCE_ORGANISM": "schistosoma_mansoni",
      "STATUS": "Complete",
      "ASSEMBLY_USED": "ASM23792v2",
      "ENA_LAST_UPDATED": "Tue Apr 15 2014 07:18:23",
      "LAST_PROCESSED_DATE": "Mon Dec 25 2017 09:24:43",
      "CRAM_LOCATION": "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR922/SRR922067/SRR922067.cram",
      "BEDGRAPH_LOCATION": "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR922/SRR922067/SRR922067.bedgraph",
      "BIGWIG_LOCATION": "ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR922/SRR922067/SRR922067.bw",
      "MAPPING_QUALITY": 95
    }
  ],
  [
    {
      "STUDY_ID": "SRP026308",
      "RUN_ID": "SRR922067",
      "TYPE": "Sample Name",
      "VALUE": "Miracidia",
      "EFO_URL": "NA"
    },
    {
      "STUDY_ID": "SRP026308",
      "RUN_ID": "SRR922067",
      "TYPE": "development stage",
      "VALUE": "miracidia",
      "EFO_URL": "NA"
    },
    {
      "STUDY_ID": "SRP026308",
      "RUN_ID": "SRR922067",
      "TYPE": "organism part",
      "VALUE": "whole body",
      "EFO_URL": "NA"
    },
    {
      "STUDY_ID": "SRP026308",
      "RUN_ID": "SRR922067",
      "TYPE": "source name",
      "VALUE": "Miracidia",
      "EFO_URL": "NA"
    },
    {
      "STUDY_ID": "SRP026308",
      "RUN_ID": "SRR922067",
      "TYPE": "strain",
      "VALUE": "NMRI",
      "EFO_URL": "NA"
    },
    {
      "STUDY_ID": "SRP026308",
      "RUN_ID": "SRR922068",
      "TYPE": "Sample Name",
      "VALUE": "Sporocyst",
      "EFO_URL": "NA"
    },
    {
      "STUDY_ID": "SRP026308",
      "RUN_ID": "SRR922068",
      "TYPE": "development stage",
      "VALUE": "48 h post-transformation sporocyst",
      "EFO_URL": "NA"
    },
    {
      "STUDY_ID": "SRP026308",
      "RUN_ID": "SRR922068",
      "TYPE": "organism part",
      "VALUE": "whole body",
      "EFO_URL": "NA"
    },
    {
      "STUDY_ID": "SRP026308",
      "RUN_ID": "SRR922068",
      "TYPE": "source name",
      "VALUE": "Sporocyst",
      "EFO_URL": "NA"
    },
    {
      "STUDY_ID": "SRP026308",
      "RUN_ID": "SRR922068",
      "TYPE": "strain",
      "VALUE": "NMRI",
      "EFO_URL": "NA"
    },
        {
      "STUDY_ID": "SRP026308",
      "RUN_ID": "SRR922067",
      "TYPE": "Dummy",
      "VALUE": "not applicable",
      "EFO_URL": "NA"
    }
  ]
]
