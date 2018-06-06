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
print $templates->{ENSEMBL_FAKE_SOURCES};
my %core_dbs_per_species;
for my $core_db (@core_dbs){
   my ($spe, $cies, $bioproject) = split "_", $core_db;
   next if $bioproject eq "core";
   my $species = "${spe}_${cies}";
   $core_dbs_per_species{$species} //= [];
   push $core_dbs_per_species{$species}, $core_db;
}

my $has_wormbase_parsers;
for my $species (keys %core_dbs_per_species){
   my @core_dbs = @{$core_dbs_per_species{$species}};
   my ($spe, $cies) = split "_", $species;
   my $wormbase_species = substr ($spe, 0, 1 ) . "_" . $cies;

   my %taxons;
   my %wormbase_annotation_paths;
   my %aliases;
   for my $core_db (@core_dbs) {
      my ($_spe, $_cies, $bioproject) = split "_", $core_db;
      $taxons{ProductionMysql->staging->meta_value($core_db, "species.taxonomy_id")}++;
      my $possible_wormbase_annotation_path = 
       join ("/",
         "/ebi/ftp/pub/databases/wormbase/releases",
         "WS$ENV{WORMBASE_VERSION}",
         "species",
         $wormbase_species,
         uc($bioproject),
         "annotation",
         join(".", $wormbase_species, uc($bioproject), "WS$ENV{WORMBASE_VERSION}", "xrefs.txt.gz")
       );
      $wormbase_annotation_paths{$possible_wormbase_annotation_path} ++ if -f $possible_wormbase_annotation_path;
     $aliases{
       "${spe}_${cies}_${bioproject}"
     }++;
   }
   my ($taxon, @other_taxons) = keys %taxons;
   die @core_dbs if @other_taxons;
   my $aliases= join ", ", keys %aliases;

   print "[species $species]\n";
   print "aliases         = $aliases\n";
   print "taxonomy_id     = $taxon\n";
   print $templates->{STANDARD_SOURCES};
   if (%wormbase_annotation_paths){
       my ($wormbase_annotation_path, @other_wormbase_annotation_paths) = keys %wormbase_annotation_paths;
       die %wormbase_annotation_paths if @other_wormbase_annotation_paths;
       $has_wormbase_parsers++;
       (my $ftp_path = $wormbase_annotation_path) =~ s/\/ebi\/ftp/ftp:\/\/ftp.ebi.ac.uk/;
       my $spe_1 = substr ($spe, 0, 1 );
       my $source =  "wormbase::$spe_1$cies";
       print "taxon = wormbase\n";
       print "source = $source\n";
       print "\n";
       print "[source $source]\n";
       print "name            = wormbase_$spe_1$cies\n";
       print $templates->{WORMBASE_SOURCE};
       print "data_uri = $ftp_path\n";
   } else {
       print "taxon = parasite\n";
   }
   print "\n";
}
print $templates->{WORMBASE_FAKE_SOURCES} if $has_wormbase_parsers;

1;
__DATA__

BEGIN_STANDARD_SOURCES_TEMPLATE
source          = EntrezGene::MULTI
source          = RNACentral::MULTI
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

[source RNACentral::MULTI]
name        = RNAcentral
download    = Y
order       = 1
priority    = 1
prio_descr  =
parser      = ChecksumParser
release_uri =
data_uri    = ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/md5/md5.tsv.gz

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

BEGIN_ENSEMBL_FAKE_SOURCES_TEMPLATE
[source RefSeq_dna::MULTI-predicted]
name            = RefSeq_dna_predicted
download        = N
order           = 20
priority        = 1
prio_descr      = refseq
parser          = RefSeqParser
release_uri     =

[source RefSeq_mRNA::MULTI]
name            = RefSeq_mRNA
download        = N
order           = 20
priority        = 3
prio_descr      = refseq
parser          = RefSeqParser
release_uri     =

[source RefSeq_mRNA::CCDS]
name            = RefSeq_mRNA
download        = N
order           = 20
priority        = 2
prio_descr      = ccds
parser          = RefSeqParser
release_uri     =

[source RefSeq_mRNA::otherfeatures]
name            = RefSeq_mRNA
download        = N
order           = 20
priority        = 1
prio_descr      = otherfeatures
parser          = RefSeqCoordinateParser
release_uri     =

[source RefSeq_peptide::otherfeatures]
name            = RefSeq_peptide
download        = N
order           = 20
priority        = 1
prio_descr      = otherfeatures
parser          = RefSeqCoordinateParser
release_uri     =

[source RefSeq_mRNA_predicted::otherfeatures]
name            = RefSeq_mRNA_predicted
download        = N
order           = 20
priority        = 1
prio_descr      = otherfeatures
parser          = RefSeqCoordinateParser
release_uri     =

[source RefSeq_peptide_predicted::otherfeatures]
name            = RefSeq_peptide_predicted
download        = N
order           = 20
priority        = 1
prio_descr      = otherfeatures
parser          = RefSeqCoordinateParser
release_uri     =

[source RefSeq_ncRNA::MULTI]
name            = RefSeq_ncRNA
download        = N
order           = 20
priority        = 2
prio_descr      =
parser          = RefSeqParser
release_uri     =

[source RefSeq_ncRNA::otherfeatures]
name            = RefSeq_ncRNA
download        = N
order           = 20
priority        = 1
prio_descr      = otherfeatures
parser          = RefSeqCoordinateParser
release_uri     =

[source RefSeq_ncRNA_predicted::otherfeatures]
name            = RefSeq_ncRNA_predicted
download        = N
order           = 20
priority        = 1
prio_descr      = otherfeatures
parser          = RefSeqCoordinateParser
release_uri     =

[source RefSeq_mRNA_predicted::MULTI]
name            = RefSeq_mRNA_predicted
download        = N
order           = 20
priority        = 2
prio_descr      = refseq
parser          = RefSeqParser
release_uri     =

[source RefSeq_mRNA_predicted::CCDS]
name            = RefSeq_mRNA_predicted
download        = N
order           = 20
priority        = 1
prio_descr      = ccds
parser          = RefSeqParser
release_uri     =

[source RefSeq_ncRNA_predicted::MULTI]
name            = RefSeq_ncRNA_predicted
download        = N
order           = 20
priority        = 1
prio_descr      =
parser          = RefSeqParser
release_uri     =

[source RefSeq_peptide::MULTI-predicted]
name            = RefSeq_peptide_predicted
download        = N
order           = 30
priority        = 2
prio_descr      =
parser          = RefSeqGPFFParser
release_uri     =

END_ENSEMBL_FAKE_SOURCES_TEMPLATE
