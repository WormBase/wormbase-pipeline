#!/usr/bin/env perl
use ProductionMysql;
sub read_templates {
  my (%templates, $current);
  while(<DATA>) {
    if (/^BEGIN_(\S+)_TEMPLATE/) {
      $current = $1;
      next;
    }

    if (/^END_\S+\_TEMPLATE/) {
      $current = undef;
      next;
    }

    $templates{$current} .= $_ if defined $current;
  }

  return \%templates;
}

my $templates = &read_templates;

my @all_core_dbs = ProductionMysql->staging->core_databases;
@ARGV = @all_core_dbs unless @ARGV;

my @core_dbs;
for my $core_db (@all_core_dbs){
  my $include;
  for my $pat (@ARGV){
    $include = 1 if $core_db =~ /$pat/;
  }
  push @core_dbs, $core_db if $include;
}
die "Usage: $0 <core_dbs_pattern" unless @core_dbs;

print $templates->{ENSEMBL_PARSERS};
my $has_wormbase_parsers;
for my $core_db (@core_dbs){
   my ($spe, $cies, $bioproject) = split "_", $core_db;
   my $species = "${spe}_${cies}";
   my $wormbase_species = substr ($spe, 0, 1 ) . "_" . $cies;
   my $alias = "${spe}_${cies}_${bioproject}";
   my $taxon = ProductionMysql->staging->meta_value($core_db, "species.taxonomy_id");
   my $wormbase_annotation_path = join ("/",   
     "/ebi/ftp/pub/databases/wormbase/releases",
     "WS$ENV{WORMBASE_VERSION}",
     "species",
     $wormbase_species, 
     uc($bioproject),
     "annotation",
     join(".", $wormbase_species, uc($bioproject), "WS$ENV{WORMBASE_VERSION}", "xrefs.txt.gz"),
   );
   print "[species $species]\n";
   print "alias           = $alias\n";
   print "taxonomy_id     = $taxon\n";
   print $templates->{STANDARD_SOURCES};
   if (-f $wormbase_annotation_path){
       $has_wormbase_parsers++;
       (my $ftp_path = $wormbase_annotation_path) =~ s/\/ebi\/ftp/ftp:\/\/ftp.ebi.ac.uk/;
       my $spe_1 = substr ($spe, 0, 1 );
       my $source =  "wormbase::$spe_1$cies";
       print "source = $source\n";
       print "\n";
       print "[source $source]\n";
       print "name            = wormbase_$spe_1$cies\n";
       print $templates->{WORMBASE_SOURCE};
       print "data_uri = $ftp_path\n";
   }
   print "\n";
}
print $templates->{WORMBASE_FAKE_SOURCES} if $has_wormbase_parsers;

1;
__DATA__

BEGIN_STANDARD_SOURCES_TEMPLATE
source          = EntrezGene::MULTI
source          = GO::MULTI
source          = RefSeq_dna::MULTI-invertebrate
source          = RefSeq_peptide::MULTI-invertebrate
source          = Uniprot/SPTREMBL::MULTI-invertebrate
source          = Uniprot/SWISSPROT::MULTI-invertebrate
source          = UniParc::MULTI
END_STANDARD_SOURCES_TEMPLATE

BEGIN_WORMBASE_SOURCE_TEMPLATE
download        = Y
order           = 50
priority        = 1
prio_descr      =
parser          = WormbaseDirectParser
release_uri     =
END_WORMBASE_SOURCE_TEMPLATE

BEGIN_ENSEMBL_PARSERS_TEMPLATE

[source EntrezGene::MULTI]
name            = EntrezGene
download        = Y
order           = 10
priority        = 1
prio_descr      =
parser          = EntrezGeneParser
release_uri     =
data_uri        = ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz

[source WikiGene::MULTI]
name            = WikiGene
download        = N
order           = 100
priority        = 1
prio_descr      =
parser          = EntrezGeneParser
release_uri     =
data_uri        = comes via EntrezGene

[source GO::MULTI]
name            = GO
download        = Y
order           = 80
priority        = 1
prio_descr      = main
parser          = GOParser
dependent_on    = Uniprot/SPTREMBL,Uniprot/SWISSPROT,RefSeq_dna,RefSeq_peptide,SGD
release_uri     = ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/current_release_numbers.txt
data_uri        = ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz

[source RefSeq_dna::MULTI-invertebrate]
name            = RefSeq_dna
download        = Y
order           = 20
priority        = 1
prio_descr      = refseq
parser          = RefSeqParser
release_uri     = ftp://ftp.ncbi.nih.gov/refseq/release/release-notes/RefSeq-release*.txt
data_uri        = ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate*.rna.fna.gz

[source RefSeq_peptide::MULTI-invertebrate]
name            = RefSeq_peptide
download        = Y
order           = 30
priority        = 2
prio_descr      =
parser          = RefSeqGPFFParser
release_uri     = ftp://ftp.ncbi.nih.gov/refseq/release/release-notes/RefSeq-release*.txt
data_uri        = ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/invertebrate*.protein.gpff.gz

[source Uniprot/SPTREMBL::MULTI-invertebrate]
name            = Uniprot/SPTREMBL
download        = Y
order           = 20
priority        = 3
prio_descr      =
parser          = UniProtParser
dependent_on    = MIM
release_uri     = ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/reldate.txt
data_uri        = ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_invertebrates.dat.gz

[source Uniprot/SWISSPROT::MULTI-invertebrate]
name            = Uniprot/SWISSPROT
download        = Y
order           = 20
priority        = 3
prio_descr      = sequence_mapped
parser          = UniProtParser
dependent_on    = MIM
release_uri     = ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/reldate.txt
data_uri        = ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_invertebrates.dat.gz

[source UniParc::MULTI]
name        = UniParc
download    = Y
order       = 20
priority    = 1
prio_descr  =
parser      = ChecksumParser
release_uri =
data_uri    = ftp://ftp.ebi.ac.uk/pub/contrib/uniparc/upidump.lis

END_ENSEMBL_PARSERS_TEMPLATE

BEGIN_WORMBASE_FAKE_SOURCES_TEMPLATE
[source wormpep_id::wormbase]
name            = wormpep_id
download        = N
order           = 50
priority        = 1
prio_descr      =
parser          = comes from WormbaseDirectParser
release_uri     =
data_uri        =

[source wormbase_gene::wormbase]
name            = wormbase_gene
download        = N
order           = 50
priority        = 1
prio_descr      =
parser          = comes from WormbaseDirectParser
release_uri     =
data_uri        =

[source wormbase_locus::wormbase]
name            = wormbase_locus
download        = N
order           = 50
priority        = 1
prio_descr      =
parser          = comes from WormbaseDirectParser
release_uri     =
data_uri        =

[source wormbase_gseqname::wormbase]
name            = wormbase_gseqname
download        = N
order           = 50
priority        = 1
prio_descr      =
parser          = comes from WormbaseDirectParser
release_uri     =
data_uri        =

[source wormbase_transcript::wormbase]
name            = wormbase_transcript
download        = N
order           = 50
priority        = 1
prio_descr      =
parser          = comes from WormbaseDirectParser
release_uri     =
data_uri        =
END_WORMBASE_FAKE_SOURCES_TEMPLATE
