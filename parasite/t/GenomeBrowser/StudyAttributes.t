use Test::More;
use Test::MockModule;

use File::Temp qw/tempdir/;
use XML::Simple;
use GenomeBrowser::RnaseqerMetadata;
use GenomeBrowser::StudyAttributes;
my $study_id = "ERP006623";

my $module = new Test::MockModule('GenomeBrowser::LocallyCachedResource');
my $xml    = XMLin(
    do { local $/; <DATA> }
);

my $imported_study_id  = "SRP124650";
my $imported_study_str = <<EOF;
<?xml version="1.0" encoding="UTF-8"?>
<ROOT request="SRP124650&amp;display=xml">
<STUDY accession="SRP124650" alias="GSE106693" center_name="GEO" broker_name="NCBI">
     <IDENTIFIERS>
          <PRIMARY_ID>SRP124650</PRIMARY_ID>
          <SECONDARY_ID>PRJNA417697</SECONDARY_ID>
          <EXTERNAL_ID label="primary" namespace="BioProject">PRJNA417697</EXTERNAL_ID>
          <EXTERNAL_ID namespace="GEO">GSE106693</EXTERNAL_ID>
          <SUBMITTER_ID namespace="GEO">GSE106693</SUBMITTER_ID>
     </IDENTIFIERS>
     <STUDY_ATTRIBUTES>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-FIRST-PUBLIC</TAG>
               <VALUE>2017-11-30</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-LAST-UPDATE</TAG>
               <VALUE>2017-11-30</VALUE>
          </STUDY_ATTRIBUTE>
     </STUDY_ATTRIBUTES>
</STUDY>
</ROOT>
EOF
my $imported_study_xml = XMLin($imported_study_str);

$module->mock(
    'get_xml',
    sub {
        my ( $class, $url ) = @_;
        if ( $url =~ /$study_id/ ) {
            return $xml;
        }
        elsif ( $url =~ /$imported_study_id/ ) {
            return $imported_study_xml;
        }
        else {
            return undef;
        }
    }
);

my $root_dir          = tempdir( CLEANUP => 1 );
my $species           = "fasciola_hepatica";
my $assembly          = "Fasciola_10x_pilon";
my $rnaseqer_metadata = bless {
    metadata => {
        $assembly => {
            $study_id          => {},
            $imported_study_id => {},
        }
    }
  },
  'GenomeBrowser::RnaseqerMetadata';

my $subject = GenomeBrowser::StudyAttributes->new( $root_dir, $species,
    $rnaseqer_metadata );
is_deeply(
    $subject,
    (
        bless {
            $assembly => {
                $study_id => {

                    #lowercase because it's a facet
                    "study" =>
"ERP006623: Some RNA-seq reads form different developmental stages of the liver fluke Fasciola hepatica",
                    "Study description" =>
"RNA was prepared from various stages of the liver fluke Fasciola hepatica by John Dalton's group and sequenced by Genome Quebec.",
                    "PubMed" =>
'<a href="https://www.ncbi.nlm.nih.gov/pubmed/25887684">25887684</a>',
                    "ArrayExpress" =>
'<a href="http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-451">E-MTAB-451 in ArrayExpress</a>',
                    "ENA first public"  => "2014-12-31",
                    "ENA last update"   => "2016-04-19",
                    'submitting_centre' => 'University of Liverpool'
                },
                'SRP124650' => {
                    'ENA first public' => '2017-11-30',
                    'ENA last update'  => '2017-11-30',
                    'study'            => 'SRP124650'
                }
            }
        },
        'GenomeBrowser::StudyAttributes'
    ),
    "Get properties from ENA xml"
) or diag explain $subject;

is_deeply(
    GenomeBrowser::StudyAttributes::_clean_messy_text(
        "Globodera_pallida_transcriptomics"),
    "Globodera pallida transcriptomics",
    "Names get better"
);

is_deeply( GenomeBrowser::StudyAttributes::_clean_messy_text("MIYAZAKI"),
    "Miyazaki", "Names get better 2" );
done_testing();
__DATA__
<?xml version="1.0" encoding="UTF-8"?>
<ROOT request="ERP006623&amp;display=xml">
<STUDY alias="ena-STUDY-LIV-10-08-2014-19:12:12:707-88" center_name="University of Liverpool" accession="ERP006623">
     <IDENTIFIERS>
          <PRIMARY_ID>ERP006623</PRIMARY_ID>
          <SECONDARY_ID>PRJEB6948</SECONDARY_ID>
          <SUBMITTER_ID namespace="LIV">ena-STUDY-LIV-10-08-2014-19:12:12:707-88</SUBMITTER_ID>
     </IDENTIFIERS>
     <DESCRIPTOR>
          <STUDY_TITLE>Some RNA-seq reads form different developmental stages of the liver fluke Fasciola hepatica</STUDY_TITLE>
          <STUDY_ABSTRACT>RNA was prepared from various stages of the liver fluke Fasciola hepatica by John Dalton's group and sequenced by Genome Quebec.</STUDY_ABSTRACT>
          <STUDY_DESCRIPTION>RNA was prepared from various stages of the liver fluke Fasciola hepatica by John Dalton's group and sequenced by Genome Quebec.</STUDY_DESCRIPTION>
          <CENTER_PROJECT_NAME>Transcriptome of Fasciola hepatica</CENTER_PROJECT_NAME>
          <STUDY_TYPE existing_study_type="Other"/>
     </DESCRIPTOR>
     <STUDY_LINKS>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>PUBMED</DB>
                    <ID>25887684</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <URL_LINK>
                    <LABEL>E-MTAB-451 in ArrayExpress</LABEL>
                    <URL>http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-451</URL>
               </URL_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-SAMPLE</DB>
                    <ID>ERS524684,ERS524692,ERS524694-ERS524695,ERS524697</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-EXPERIMENT</DB>
                    <ID>ERX535559-ERX535563</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-RUN</DB>
                    <ID>ERR577156-ERR577160</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMISSION</DB>
                    <ID>ERA345947</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-FASTQ-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=ERP006623&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes]]></ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMITTED-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=ERP006623&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format]]></ID>
               </XREF_LINK>
          </STUDY_LINK>
     </STUDY_LINKS>
     <STUDY_ATTRIBUTES>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-SPOT-COUNT</TAG>
               <VALUE>232639591</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-BASE-COUNT</TAG>
               <VALUE>46527918200</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-FIRST-PUBLIC</TAG>
               <VALUE>2014-12-31</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-LAST-UPDATE</TAG>
               <VALUE>2016-04-19</VALUE>
          </STUDY_ATTRIBUTE>
     </STUDY_ATTRIBUTES>
</STUDY>
</ROOT>
