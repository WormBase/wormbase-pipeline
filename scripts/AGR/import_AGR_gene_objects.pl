#!/bin/env perl
# trawl throught he rest endpoint and import from the BGI genes into WormBase
# potential issues:
#  * the REST endpoint will change regularly (or its payload)
#  * BGI schema will also change regulary (https://github.com/alliance-genome/agr_schemas/tree/master/ingest/gene)
#
# URLs:
# * https://fms.alliancegenome.org/api/datafile/by/BGI?latest=true
# * https://download.alliancegenome.org


use LWP::Simple;
use JSON;
use Wormbase;
use Log_files;
use Storable;
use IO::File;
use Getopt::Long;
use Compress::Zlib;
use strict;

my ($debug,$test,$store,$wormbase,$load,$outfile);
GetOptions( 'debug=s'      => \$debug,
            'test'         => \$test,
            'store=s'      => \$store,
	    'load'         => \$load,
            'outfile=s'    => \$outfile,
            )||die("invalid commandline options\n");

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             );
}

my $log = Log_files->make_build_log($wormbase); # establish log file.

$outfile||=$wormbase->acefiles . '/AGR_genes.ace';
my $outfh = IO::File->new($outfile,'w');
$outfh->binmode(':utf8'); # that might bite us when importing into acedb
$log->write_to("writing to $outfile\n\n");

my %taxon2name = (
	'NCBITaxon:9606'   => 'Homo sapiens',            # human
	'NCBITaxon:7227'   => 'Drosophila melanogaster', # fly
	'NCBITaxon:10090'  => 'Mus musculus',            # mouse
	'NCBITaxon:10116'  => 'Rattus norvegicus',       # rat
	'NCBITaxon:7955'   => 'Danio rerio',             # fish
	'NCBITaxon:559292' => 'Saccharomyces cerevisiae',# mushroom
);

my %so2term = ( # probably shouldn't be hardcoded. We could altermatively pull it from ACeDB
	'SO:0000704' => 'gene',
	'SO:0001263' => 'ncRNA_gene',
	'SO:0001637' => 'rRNA_gene',
	'SO:0001217' => 'protein_coding_gene',
	'SO:0000336' => 'pseudogene',
	'SO:0001268' => 'snRNA_gene',
	'SO:0001272' => 'tRNA_gene',
	'SO:0001267' => 'snoRNA_gene',
	'SO:0001265' => 'miRNA_gene',
	'SO:0001500' => 'heritable_phenotypic_marker',
	'SO:0001269' => 'SRP_RNA_gene',
	'SO:0002182' => 'antisense_lncRNA_gene',
	'SO:0001841' => 'polymorphic_pseudogene',
	'SO:0002127' => 'lncRNA_gene',
	'SO:3000000' => 'gene_segment',
	'SO:0001741' => 'pseudogenic_gene_segment',
	'SO:0001639' => 'RNase_P_RNA_gene',
	'SO:0001266' => 'scRNA_gene',
	'SO:0001643' => 'telomerase_RNA_gene',
        'SO:0002185' => 'bidirectional_promoter_lncRNA',
	'SO:0002184' => 'sense_intronic_ncRNA_gene',
	'SO:0002183' => 'sense_overlap_ncRNA_gene',
	'SO:0001641' => 'lincRNA_gene',
	'SO:0002181' => 'ribozyme_gene',
	'SO:0001640' => 'RNase_MRP_RNA_gene',
	'SO:0000111' => 'transposable_element_gene',
	'SO:0000718' => 'blocked_reading_frame',
	'SO:0001411' => 'biological_region',
);

my %prefix2db = (
	HUMAN => 'HGNC',
	FB    => 'FlyBase', # the only one that is different
	RGD   => 'RGD',
	MGI   => 'MGI',
	ZFIN  => 'ZFIN',
	SGD   => 'SGD',
);

my $content = get('https://fms.alliancegenome.org/api/datafile/by/BGI?latest=true');
die "Couldn't query the REST endpoint\n" unless $content;

my $latest = decode_json($content);
# datasubType->name [FB,MGI,RGD,SGD,ZFIN,HUMAN]
foreach my $db (@$latest){
	if ($db->{dataSubType}->{name} =~/^(FB|MGI|RGD|SGD|ZFIN|HUMAN)$/){
		my $downloadPath='https://download.alliancegenome.org/'.$db->{s3Path};
		$log->write_to("$1 => $downloadPath\n");
		get_bgi($downloadPath, $db->{dataSubType}->{name});
	}
}

$wormbase->load_to_database($wormbase->autoace,$outfile,'AGR_import',$log) if $load;
$log->mail();

# that below creates the ACE object
sub get_bgi{
	my ($url,$prefix)=@_;
	my $payload = get($url);

	$log->log_and_die("couldn't get data for $prefix\n") unless $payload;

	$payload = Compress::Zlib::memGunzip($payload) if $url =~ /\.gz$/; # If clause only required in interim while some AGR files are still unzipped

	my $json = decode_json($payload);
	foreach my $gene (@{$json->{data}}){
	        my $bgE=$gene->{basicGeneticEntity};

		print $outfh "Gene : \"${\$bgE->{primaryId}}\"\n";
		print $outfh "Database $prefix2db{$prefix} id \"".(split(":",$bgE->{primaryId}))[1]."\"\n";
		print $outfh "Database AGR cURI \"${\$bgE->{primaryId}}\"\n";

		foreach my $x (@{$bgE->{crossReferences}}){
			print $outfh "Database \"EnsEMBL\" \"ENSEMBL_geneID\" \"$1\"\n" if $x->{id}=~/ENSEMBL:(\w+)/;
			print $outfh "Database \"UniProt\" \"UniProt_AC\" \"$1\"\n" if $x->{id}=~/UniProtKB:(\w+)/;
	                print $outfh "Database \"EntrezGene\" \"EntrezGene_id\" \"$1\"\n" if $x->{id}=~/NCBI_gene:(\w+)/;
		}

		foreach my $oId (@{$bgE->{secondaryIds}}){
			print $outfh "Other_name \"$oId\"\n";
		}
		
		map {print $outfh "Other_name \"$_\"\n"} @{$bgE->{synonyms}} if $bgE->{taxonId} eq 'NCBITaxon:9606'; # only for human

		print $outfh "Live\n";
		$log->log_and_die("cannot find species ${\$bgE->{taxonId}}") unless($taxon2name{$bgE->{taxonId}});

		print $outfh 'Species "'.$taxon2name{$bgE->{taxonId}}."\"\n";

		if ($so2term{$gene->{soTermId}}){
			print $outfh 'Biotype "'.$so2term{$gene->{soTermId}}."\"\n";
		}else{
			$log->log_and_die("cannot find term for ${\$gene->{soTermId}}\n")
		}

	        print $outfh 'CGC_name "'.$gene->{symbol}."\" Inferred_automatically \"AGR_import\"\n" if $gene->{symbol};

		if ($gene->{symbol}){
		  	print $outfh 'Public_name "'.$gene->{symbol} . "\"\n";
                }else{
                        print $outfh 'Public_name "'.$bgE->{name} ."\"\n";
		}

		# $gene->{geneSynopsis}=~s/[^\x00-\x7f]//g; # strip all non-ASCII characters out of the string
		$gene->{geneSynopsis}=~s/"/\\"/g;    # FB:FBgn0001150 - escape quotation marks
		$gene->{geneSynopsis}=~s/[\n\r]/ /g; # SGD:S000001091 - replace line-breaks with white-spaces
		$gene->{geneSynopsis}=~s/\s+$//g;    # SGD:S000001091 - remove trailing whitespaces

		print $outfh 'Automated_description "'.$gene->{geneSynopsis}."\" Inferred_automatically \"AGR_import\"\n" if $gene->{geneSynopsis};
                print $outfh "\n";
	}
}
