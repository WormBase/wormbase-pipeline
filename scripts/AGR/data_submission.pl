#!/usr/bin/env perl
use strict;
use warnings;

use Log_files;
use Getopt::Long;
use File::Path 'make_path';
use JSON;
use Wormbase;

my ($fasta, $gff, $bgi, $disease, $allele, $phenotype, $expression, $ltp_variations, $htp_variations, $agm,
    $construct, $interactions, $genetic_interactions, $hts, $reference, $molecule);
my ($all, $test, $logfile, $bgi_file, $allele_file, $disease_file, $fasta_file, $gff_file, $debug);


GetOptions(
    "fasta"                => \$fasta,
    "gff"                  => \$gff,
    "bgi"                  => \$bgi,
    "disease"              => \$disease,
    "allele"               => \$allele,
    "phenotype"            => \$phenotype,
    "expression"           => \$expression,
    "ltp_variations"       => \$ltp_variations,
    "htp_variations"       => \$htp_variations,
    "agm"                  => \$agm,
    "construct"            => \$construct,
    "interactions"         => \$interactions,
    "genetic_interactions" => \$genetic_interactions,
    "hts"                  => \$hts,
    "reference"            => \$reference,
    "molecule"             => \$molecule,
    "all"                  => \$all,
    "test"                 => \$test,
    "logfile=s"            => \$logfile,
    "bgi_file=s"           => \$bgi_file,
    "allele_file=s"        => \$allele_file,
    "disease_file=s"       => \$disease_file,
    "fasta_file=s"         => \$fasta_file,
    "gff_file=s"           => \$gff_file,
    "debug=s"              => \$debug
    );

my $agr_resources   = $ENV{'AGR_DIR'} . '/resources';
my $agr_uploads     = $ENV{'AGR_UPLOADS'};
my $build_home      = $ENV{'BUILD_HOME'};
my $agr_release     = $ENV{'AGR_RELEASE'};
my $agr_schema      = $ENV{'AGR_SCHEMA'};
my $ws_release      = $ENV{'WS_RELEASE'};
my $upload_date     = $ENV{'UPLOAD_DATE'};
my $cvs_dir         = $ENV{'CVS_DIR'};
my $agr_schema_repo = $ENV{'AGR_SCHEMA_REPO'};
my $sub_dir = "${agr_uploads}/${agr_release}/${upload_date}";

$bgi_file = "${sub_dir}/WB_${agr_schema}_BGI.json" unless $bgi_file;
if (!$bgi and ($gff or $allele or $phenotype or $expression)) {
    die "BGI file ${bgi_file} doesn't exist, please specify using --bgi_file\n" unless -e $bgi_file;
}


$allele_file = "${sub_dir}/WB_${agr_schema}_allele.json" unless $allele_file;
if (!$allele and $agm) {
    die "Alleles file ${allele_file} doesn't exist, please specify using --allele_file\n" unless -e $allele_file;
}

$disease_file = "${sub_dir}/WB_${agr_schema}_disease.json" unless $disease_file;
if ($disease and $agm) {
    die "Disease file ${disease_file} doesn't exist, please specify using --disease_file\n" unless -e $disease_file;
}

$fasta_file = "${sub_dir}/WB_${agr_schema}_FASTA.fa" unless $fasta_file;
if (!$fasta and ($ltp_variations or $htp_variations)) {
    die "FASTA file ${fasta_file} doesn't exist, please specify using --fasta_file\n" unless -e $fasta_file;
}

$gff_file = "${sub_dir}/WB_${agr_schema}_gff3.gff3" unless $gff_file;
if (!$gff and $htp_variations) {
    die "GFF file ${gff_file} doesn't exist, please specify using --gff_file\n" unless -e $gff_file;
}

make_path($sub_dir);    

$logfile = "${sub_dir}/data_submission.log" unless $logfile;
my $log = Log_files->make_log($logfile, $debug);
$log->write_to("Starting AGR submission of $ws_release data\n");

if ($fasta or $all) {
    my @cmds = (
	"cp ${build_home}/DATABASES/${ws_release}/SEQUENCES/elegans.genome.fa ${sub_dir}/WB_${agr_schema}_FASTA.fa",
	"perl -i -pne 's/CHROMOSOME_//' ${fasta_file}"
	);
    my $datatype_processed = process_datatype('FASTA', \@cmds, $log);
    submit_data('FASTA', $fasta_file, $log)
	if $datatype_processed and !$test;
}

if ($bgi or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_basic_gene_json.pl -database ${build_home}/DATABASES/${ws_release} " .
	"-gtf ${build_home}/DATABASES/${ws_release}/SEQUENCES/elegans.gtf -wsversion ${ws_release} " .
	"-outfile ${bgi_file}",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/gene/geneMetaData.json " .
	"-d ${bgi_file}"
	);
    my $datatype_processed = process_datatype('BGI', \@cmds, $log);
    submit_data('BGI', $bgi_file, $log)
	if $datatype_processed and !$test;
}

if ($disease or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_disease_json.pl -database ${build_home}/DATABASES/${ws_release} " . 
	"-wsversion ${ws_release} -outfile ${disease_file} -AGRwhitelist",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/disease/diseaseMetaDataDefinition.json " .
	"-d ${disease_file}"
	);
    my $datatype_processed = process_datatype('disease', \@cmds, $log);
    submit_data('DAF', $disease_file, $log)
	if $datatype_processed and !$test;
}

if ($gff or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_gff.pl  -wsversion ${ws_release} -bgijson ${bgi_file} " . 
	"-gffin ${build_home}/DATABASES/${ws_release}/SEQUENCES/elegans.processed.gff3.gz " .
	"-gffout ${gff_file}"
	);
    my $datatype_processed = process_datatype('GFF', \@cmds, $log);
    submit_data('GFF', $gff_file, $log)
	if $datatype_processed and !$test;
}

if ($allele or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_allele_json.pl -database ${build_home}/DATABASES/${ws_release} " .
	"-wsversion ${ws_release} -bgijson ${bgi_file} -outfile ${allele_file}",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/allele/alleleMetaData.json " .
	"-d ${allele_file}"
	);
    my $datatype_processed = process_datatype('allele', \@cmds, $log);
    submit_data('ALLELE', $allele_file, $log)
	if $datatype_processed and !$test;
}

if ($phenotype or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_phenotype_json.pl -database ${build_home}/DATABASES/${ws_release} " .
	"-wsversion ${ws_release} -bgijson ${bgi_file} -outfile ${sub_dir}/WB_${agr_schema}_phenotype.json",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/phenotype/phenotypeMetaDataDefinition.json " .
	"-d ${sub_dir}/WB_${agr_schema}_phenotype.json"
	);
    my $datatype_processed = process_datatype('phenotype', \@cmds, $log);
    submit_data('PHENOTYPE', "${sub_dir}/WB_${agr_schema}_phenotype.json", $log)
	if $datatype_processed and !$test;
}

if ($expression or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_expression_json.pl -database ${build_home}/DATABASES/${ws_release} " .
	"-wsversion ${ws_release} -bgijson ${bgi_file}  -wb2uberon ${agr_resources}/wormbase_to_uberon.agr_2_2.txt " .
	"-outfile ${sub_dir}/WB_${agr_schema}_expression.json",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/expression/wildtypeExpressionMetaDataDefinition.json " .
	"-d ${sub_dir}/WB_${agr_schema}_expression.json"
	);
    my $datatype_processed = process_datatype('expression', \@cmds, $log);
    submit_data('EXPRESSION', "${sub_dir}/WB_${agr_schema}_expression.json", $log)
	if $datatype_processed and !$test;
}

if ($ltp_variations or $all) {
    my @cmds = (
	"ruby ${cvs_dir}/AGR/make_agr_variations.rb -f ${fasta_file} " .
	"-g ${build_home}/DATABASES/${ws_release}/SEQUENCES/elegans.processed.gff3.gz -d ${build_home}/DATABASES/${ws_release} " .
	"-w ${ws_release} -q ${cvs_dir}/AGR/agr_variations.def -o ${sub_dir}/WB_${agr_schema}_variations.json",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/allele/variantMetaData.json " .
	"-d ${sub_dir}/WB_${agr_schema}_variations.json" 
	);
    my $datatype_processed = process_datatype('LTP variation', \@cmds, $log);
    submit_data('VARIATION', "${sub_dir}/WB_${agr_schema}_variations.json", $log)
	if $datatype_processed and !$test;
}

if ($htp_variations or $all) {
    my @cmds = (
	"ruby ${cvs_dir}/AGR/make_agr_variations.rb -a -f ${fasta_file} " .
	"-g ${build_home}/DATABASES/${ws_release}/SEQUENCES/elegans.processed.gff3.gz -d ${build_home}/DATABASES/${ws_release} " .
	"-w ${ws_release} -q ${cvs_dir}/AGR/agr_variations.def -o ${sub_dir}/WB_${agr_schema}_htp_variations.json",
	"python3 ${cvs_dir}/AGR/agr_variations_json2vcf.py -j ${sub_dir}/WB_${agr_schema}_htp_variations.json -g ${gff_file} " .
	"-m WB -o ${sub_dir}/WB_${agr_schema}_htp_variations.vcf -f ${fasta_file}"
	);
    my $datatype_processed = process_datatype('HTP variation', \@cmds, $log);
    submit_data('HTVCF', "${sub_dir}/WB_${agr_schema}_htp_variations.vcf", $log)
	if $datatype_processed and !$test;
}

if ($agm or $all) {
    my @cmds = (
	"perl  ${cvs_dir}/AGR/make_agr_AGM.pl -alleles ${allele_file} -database ${build_home}/DATABASES/${ws_release} " .
	"-wsversion ${ws_release} -outfile ${sub_dir}/WB_${agr_schema}_AGM.json -diseases ${disease_file}",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/affectedGenomicModel/affectedGenomicModelMetaData.json " .
	"-d ${sub_dir}/WB_${agr_schema}_AGM.json"
	);
    my $datatype_processed = process_datatype('AGM', \@cmds, $log);
    submit_data('AGM', "${sub_dir}/WB_${agr_schema}_AGM.json", $log)
	if $datatype_processed and !$test;
}

if ($construct or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_construct_json.pl -database ${build_home}/DATABASES/${ws_release} " .
	"-wsversion ${ws_release} -outfile ${sub_dir}/WB_${agr_schema}_construct.json",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/construct/constructMetaData.json " .
	"-d ${sub_dir}/WB_${agr_schema}_construct.json"
	);
    my $datatype_processed = process_datatype('construct', \@cmds, $log);
    submit_data('CONSTRUCT', "${sub_dir}/WB_${agr_schema}_construct.json", $log)
	if $datatype_processed and !$test;
}

if ($interactions or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_interactions_psimitab.pl  -database ${build_home}/DATABASES/${ws_release} " .
	"-outfile ${sub_dir}/WB_${agr_schema}_interactions.psi-mi-tab"
	);
    my $datatype_processed = process_datatype('molecular interactions', \@cmds, $log);
    submit_data('INTERACTION-SOURCE', "${sub_dir}/WB_${agr_schema}_interactions.psi-mi-tab", $log, 'MOL')
	if $datatype_processed and !$test;
}

if ($genetic_interactions or $all) {
    my @cmds = (
	"perl  ${cvs_dir}/AGR/make_agr_genetic_interaction_psimitab.pl -database ${build_home}/DATABASES/${ws_release} " . 
	"-outfile ${sub_dir}/WB_${agr_schema}_genetic_interactions.psi-mi-tab"
	);
    my $datatype_processed = process_datatype('genetic interactions', \@cmds, $log);
    submit_data('INTERACTION-SOURCE', "${sub_dir}/WB_${agr_schema}_genetic_interactions.psi-mi-tab", $log, 'GEN')
	if $datatype_processed and !$test;
}

if ($hts or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_HTS.pl  -database ${build_home}/DATABASES/${ws_release} " .
	"-sample ${sub_dir}/WB_${agr_schema}_sample.json -dataset ${sub_dir}/WB_${agr_schema}_dataset.json " .
	"-wsversion ${ws_release} -wb2uberon ${build_home}/AGR/resources/wormbase_to_uberon.agr_2_2.txt",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/htp/datasetSample/datasetSampleMetaDataDefinition.json " .
	"-d ${sub_dir}/WB_${agr_schema}_sample.json",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/htp/dataset/datasetMetaDataDefinition.json " .
	"-d ${build_home}/AGR/uploads/${agr_release}/${upload_date}/WB_${agr_schema}_dataset.json"
	);
    my $datatype_processed = process_datatype('HTS dataset and sample', \@cmds, $log);
    if ($datatype_processed and !$test) {
	submit_data('HTPDATASAMPLE', "${sub_dir}/WB_${agr_schema}_sample.json", $log);
	submit_data('HTPDATASET', "${sub_dir}/WB_${agr_schema}_dataset.json", $log);
    }
}

if ($reference or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_reference_json.pl -database ${build_home}/DATABASES/${ws_release} " .
	"-wsversion ${ws_release} -ref ${sub_dir}/WB_${agr_schema}_reference.json " .
	"-refex ${sub_dir}/WB_${agr_schema}_reference_exchange.json",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/resourcesAndReferences/referenceMetaData.json " .
	"-d ${sub_dir}/WB_${agr_schema}_reference.json",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/resourcesAndReferences/referenceExchangeMetaData.json " .
	"-d ${sub_dir}/WB_${agr_schema}_reference_exchange.json"
	);
    my $datatype_processed = process_datatype('reference', \@cmds, $log);
    if ($datatype_processed and !$test) {
	submit_data('REFERENCE', "${sub_dir}/WB_${agr_schema}_reference.json", $log);
	submit_data('REF-EXCHANGE', "${sub_dir}/WB_${agr_schema}_reference_exchange.json", $log);
    }
}

if ($molecule or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_molecules.pl -database ${build_home}/DATABASES/${ws_release} -wsversion ${ws_release} " .
	"-outfile ${sub_dir}/WB_${agr_schema}_molecule.json",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/molecules/moleculeMetaData.json " .
	"-d ${sub_dir}/WB_${agr_schema}_molecule.json"
	);
    my $datatype_processed = process_datatype('molecule', \@cmds, $log);
    submit_data('MOLECULE', "${sub_dir}/WB_${agr_schema}_molecule.json", $log)
	if $datatype_processed and !$test;
}

cleanup_unzipped($sub_dir);

$log->mail;

exit(0);


sub cleanup_unzipped {
    my $sub_dir = shift;

    my @files = glob($sub_dir . '/*');
    for my $file (@files) {
	unlink $file if -e $file . '.gz';
    }

    return;
}


sub process_datatype {
    my ($datatype, $cmds, $log) = @_;

    my @cmd_output;
    for my $cmd(@$cmds) {
	$log->write_to("Running $datatype command: $cmd\n\n");
	my $output = `$cmd 2>&1`;
	my $exit_code = $? >> 8;
	if ($exit_code != 0) {
	    $log->log_and_die("Processing of $datatype data failed with exit code ${exit_code}:\n$cmd\n$output\n\n");
	    return 0;
	}
	push @cmd_output, $output unless $cmd =~ /agr_validate.py/;
    }
    my $all_output = join("\n", @cmd_output);
    $log->write_to("Processing of $datatype data succeeded\n${all_output}\n\n");
    return 1;
}


sub submit_data {
    my ($fms_datatype, $file, $log, $fms_subtype) = @_;

    my $exit_code = system("gzip -9 -c $file > $file.gz");
    if ($exit_code != 0) {
	$log->error("Compression of $file failed with exit code $exit_code, $fms_datatype submission not completed\n\n");
	return;
    }

    my $fms_file_link = $fms_subtype ? '_WB-' . $fms_subtype . '=@' : '_WB=@';

    my $cmd = 'curl -H "Authorization: Bearer ' . $ENV{'TOKEN'} . '" -X POST ' .
	'"https://fms.alliancegenome.org/api/data/submit" -F "' . $ENV{'AGR_RELEASE'} . '_' .
	$fms_datatype . $fms_file_link . $file . '.gz"';

    my $response_json = `$cmd`;
    my $response = decode_json($response_json);
    if ($response->{status} eq 'failed') {
	$log->error("Upload of $fms_datatype failed:\n$response_json\n\n");
    }
    else {
	$log->write_to("Upload of ${fms_datatype} succeeded\n\n");
    }

    return;
}    
    

