use Test::More;
use Test::MockModule;

use File::Temp qw/tempdir/;
use JSON;
use GenomeBrowser::Resources::PubMed;

my $module = new Test::MockModule('GenomeBrowser::Resources::LocallyCachedResource');
$module->mock('get_xml', do {local $/; <DATA>});

my $species = "schistosoma_mansoni";
my $pubmed_id = "29069413";

sub assert_pubmed {
  my ($metadata, $expected, $description) = @_;
  my $result = GenomeBrowser::Resources::PubMed->new(tempdir( CLEANUP => 1 ), $species, $metadata);
  is_deeply(
     $result, bless ($expected, 'GenomeBrowser::Resources::PubMed'),
  ,$description) or diag explain $result;
}

assert_pubmed({
  rnaseqer => bless({}, 'GenomeBrowser::Resources::RnaseqerMetadata'),
  array_express => bless({}, 'GenomeBrowser::Resources::ArrayExpressMetadata'),
},{}, "Null case");
sub pick_up_paper {
  my ($metadata, $description) = @_;
  assert_pubmed($metadata, {assembly_id => {study_id => {
     $pubmed_id => ["Lee & Howe, 2018", "WormBase 2017: molting into a new stage. (Lee & Howe, 2018)"]
  }}}, $description);
}
pick_up_paper({
  rnaseqer => bless ({metadata =>{assembly_id => {study_id=>{}}}},'GenomeBrowser::Resources::RnaseqerMetadata'),
  ena => {assembly_id => {study_id=>{pubmed => [$pubmed_id]}}},
  array_express => bless({}, 'GenomeBrowser::Resources::ArrayExpressMetadata'),
}, "ENA");

pick_up_paper({
  rnaseqer => bless ({metadata => {assembly_id => {study_id=>{}}}},'GenomeBrowser::Resources::RnaseqerMetadata'),
  geo => {assembly_id => {study_id=>{pubmed => [$pubmed_id]}}},
  array_express => bless({}, 'GenomeBrowser::Resources::ArrayExpressMetadata'),
}, "GEO");

pick_up_paper({
  rnaseqer => bless ({metadata => {assembly_id => {study_id=>{}}}},'GenomeBrowser::Resources::RnaseqerMetadata'),
  array_express => bless ({
    secondary_to_primary_accession=>{study_id=>"E-MOCK-1"},
    primary_accession_to_pubmed=>{"E-MOCK-1", [$pubmed_id]}
  }, 'GenomeBrowser::Resources::ArrayExpressMetadata'),
}, "AE");
done_testing();
__DATA__
<PubmedArticle>
    <MedlineCitation Status="In-Data-Review" Owner="NLM">
        <PMID Version="1">29069413</PMID>
        <DateRevised>
            <Year>2018</Year>
            <Month>06</Month>
            <Day>07</Day>
        </DateRevised>
        <Article PubModel="Print">
            <Journal>
                <ISSN IssnType="Electronic">1362-4962</ISSN>
                <JournalIssue CitedMedium="Internet">
                    <Volume>46</Volume>
                    <Issue>D1</Issue>
                    <PubDate>
                        <Year>2018</Year>
                        <Month>Jan</Month>
                        <Day>04</Day>
                    </PubDate>
                </JournalIssue>
                <Title>Nucleic acids research</Title>
                <ISOAbbreviation>Nucleic Acids Res.</ISOAbbreviation>
            </Journal>
            <ArticleTitle>WormBase 2017: molting into a new stage.</ArticleTitle>
            <Pagination>
                <MedlinePgn>D869-D874</MedlinePgn>
            </Pagination>
            <ELocationID EIdType="doi" ValidYN="Y">10.1093/nar/gkx998</ELocationID>
            <Abstract>
                <AbstractText>WormBase (http://www.wormbase.org) is an important knowledge resource for biomedical researchers worldwide. To accommodate the ever increasing amount and complexity of research data, WormBase continues to advance its practices on data acquisition, curation and retrieval to most effectively deliver comprehensive knowledge about Caenorhabditis elegans, and genomic information about other nematodes and parasitic flatworms. Recent notable enhancements include user-directed submission of data, such as micropublication; genomic data curation and presentation, including additional genomes and JBrowse, respectively; new query tools, such as SimpleMine, Gene Enrichment Analysis; new data displays, such as the Person Lineage browser and the Summary of Ontology-based Annotations. Anticipating more rapid data growth ahead, WormBase continues the process of migrating to a cutting-edge database technology to achieve better stability, scalability, reproducibility and a faster response time. To better serve the broader research community, WormBase, with five other Model Organism Databases and The Gene Ontology project, have begun to collaborate formally as the Alliance of Genome Resources.</AbstractText>
                <CopyrightInformation>© The Author(s) 2017. Published by Oxford University Press on behalf of Nucleic Acids Research.</CopyrightInformation>
            </Abstract>
            <AuthorList CompleteYN="Y">
                <Author ValidYN="Y">
                    <LastName>Lee</LastName>
                    <ForeName>Raymond Y N</ForeName>
                    <Initials>RYN</Initials>
                    <AffiliationInfo>
                        <Affiliation>Division of Biology and Biological Engineering 156-29, California Institute of Technology, Pasadena, CA 91125, USA.</Affiliation>
                    </AffiliationInfo>
                </Author>
                <Author ValidYN="Y">
                    <LastName>Howe</LastName>
                    <ForeName>Kevin L</ForeName>
                    <Initials>KL</Initials>
                    <AffiliationInfo>
                        <Affiliation>European Molecular Biology Laboratory, European Bioinformatics Institute, Wellcome Trust Genome Campus, Hinxton, Cambridge CB10 1SD, UK.</Affiliation>
                    </AffiliationInfo>
                </Author>
            </AuthorList>
        </Article>
    </MedlineCitation>
    <PubmedData>
       ...
    </PubmedData>
</PubmedArticle>
