#!/usr/bin/env perl
use strict;
use warnings;

use Const::Fast;
use Log_files;
use Getopt::Long;
use File::Path 'make_path';
use JSON;
use Wormbase;
use Modules::WormSlurm;

const my @DATATYPES => ('FASTA', 'BGI', 'DAF', 'GFF', 'ALLELE', 'PHENOTYPE', 'EXPRESSION',
			'VARIATION', 'HTVCF', 'AGM', 'CONSTRUCT', 'INTERACTION-SOURCE_MOL',
			'INTERACTION-SOURCE_GEN', 'HTPDATASAMPLE', 'HTPDATASET', 'REFERENCE',
			'REF-EXCHANGE', 'MOLECULE');
const my $LOW_MEM => '1g';
const my $DEFAULT_MEM => '6g'; # Default amount of memory in Gb for Slurm jobs
const my $HIGH_MEM => '12g';    # Amount of memory in Gb for jobs requiring extra

my ($fasta, $gff, $bgi, $disease, $allele, $phenotype, $expression, $ltp_variations, $htp_variations, $agm,
    $construct, $interactions, $genetic_interactions, $hts, $reference, $molecule, $all, $test, $logfile,
    $bgi_file, $allele_file, $disease_file, $fasta_file, $gff_file, $chebi_map_file, $debug);


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
    "chebi_map_file=s"     => \$chebi_map_file,
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

# FASTA                               1
# BGI                                 1
# DISEASE                             1
# GFF - BGI                           1  
# ALLELE - BGI                        1
# PHENOTYPE - BGI                     1
# EXPRESSION - BGI                    1
# LTP VARIATIONS - FASTA              1
# HTP VARIATIONS - GFF, FASTA .. BGI  1
# AGM - DISEASE, ALLELE .. BGI        1
# CONSTRUCT                           2
# INTERACTIONS                        3
# GENETIC INTERACTIONS                4
# HTS                                 5
# REFERENCE                           6
# MOLECULE

my %datatypes_processed = map {$_ => 0} @DATATYPES;

$bgi_file = "${sub_dir}/WB_${agr_schema}_BGI.json" unless $bgi_file;
if (!$bgi and ($gff or $allele or $phenotype or $expression)) {
    die "BGI file ${bgi_file} doesn't exist, please specify using --bgi_file\n" unless -e $bgi_file;
}


$allele_file = "${sub_dir}/WB_${agr_schema}_allele.json" unless $allele_file;
if (!$allele and $agm) {
    die "Alleles file ${allele_file} doesn't exist, please specify using --allele_file\n" unless -e $allele_file;
}

$disease_file = "${sub_dir}/WB_${agr_schema}_disease.json" unless $disease_file;
if (!$disease and $agm) {
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

$chebi_map_file = "${sub_dir}/chebi_id_to_name_map.txt" unless $chebi_map_file;

make_path("${sub_dir}/slurm_logs");

my %files_to_submit = (
	'FASTA'                  => $fasta_file,
	'BGI'                    => $bgi_file,
	'DAF'                    => $disease_file,
	'GFF'                    => $gff_file,
	'ALLELE'                 => $allele_file,
	'PHENOTYPE'              => "${sub_dir}/WB_${agr_schema}_phenotype.json",
	'EXPRESSION'             => "${sub_dir}/WB_${agr_schema}_expression.json",
	'VARIATION'              => "${sub_dir}/WB_${agr_schema}_variations.json",
	'HTVCF'                  => "${sub_dir}/WB_${agr_schema}_htp_variations.vcf",
	'AGM'                    => "${sub_dir}/WB_${agr_schema}_AGM.json",
	'CONSTRUCT'              => "${sub_dir}/WB_${agr_schema}_construct.json",
	'INTERACTION-SOURCE_MOL' => "${sub_dir}/WB_${agr_schema}_interactions.psi-mi-tab",
	'INTERACTION-SOURCE_GEN' => "${sub_dir}/WB_${agr_schema}_genetic_interactions.psi-mi-tab",
	'HTPDATASAMPLE'          => "${sub_dir}/WB_${agr_schema}_sample.json",
	'HTPDATASET'             => "${sub_dir}/WB_${agr_schema}_dataset.json",
	'REFERENCE'              => "${sub_dir}/WB_${agr_schema}_reference.json",
	'REF-EXCHANGE'           => "${sub_dir}/WB_${agr_schema}_reference_exchange.json",
	'MOLECULE'               => "${sub_dir}/WB_${agr_schema}_molecule.json" 
    );

my ($fasta_id, $bgi_id, $disease_id, $gff_id, $allele_id, $pheno_id, $exp_id, $ltp_id, $htp_id, $agm_id,
    $construct_id, $int_id, $gen_int_id, $hts_id, $ref_id, $mol_id);

$logfile = "${sub_dir}/data_submission.log" unless $logfile;
my $log = Log_files->make_log($logfile, $debug);
$log->write_to("Starting AGR submission of $ws_release data\n");

if ($fasta or $all) {
    my @cmds = (
	"cp ${build_home}/DATABASES/${ws_release}/SEQUENCES/elegans.genome.fa ${sub_dir}/WB_${agr_schema}_FASTA.fa",
	"perl -i -pne 's/CHROMOSOME_//' ${fasta_file}"
	);
    $datatypes_processed{'FASTA'} = process_datatype('FASTA', \@cmds, $sub_dir, $log, $LOW_MEM, '00:30:00');
}

if ($bgi or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_basic_gene_json.pl -database ${build_home}/DATABASES/${ws_release} " .
	"-gtf ${build_home}/DATABASES/${ws_release}/SEQUENCES/elegans.gtf -wsversion ${ws_release} " .
	"-outfile ${bgi_file}",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/gene/geneMetaData.json " .
	"-d ${bgi_file}"
	);
    $datatypes_processed{'BGI'} = process_datatype('BGI', \@cmds, $sub_dir, $log, $DEFAULT_MEM, '1:00:00');
}

if ($disease or $phenotype or $all) {
    my @cmds = ("perl ${cvs_dir}/AGR/get_chebi_name_map.pl -o $chebi_map_file") unless -e $chebi_map_file;
    $datatypes_processed{'CHEBI'} = process_datatype('chebi_map', \@cmds, $sub_dir, $log, $HIGH_MEM, '1:00:00');
}
    

if ($disease or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_disease_json.pl -database ${build_home}/DATABASES/${ws_release} " . 
	"-wsversion ${ws_release} -outfile ${disease_file} -AGRwhitelist -chebi ${chebi_map_file}",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/disease/diseaseMetaDataDefinition.json " .
	"-d ${disease_file}"
	);
    $datatypes_processed{'DAF'} = process_datatype('disease', \@cmds, $sub_dir, $log, $DEFAULT_MEM, '1:00:00', $datatypes_processed{'CHEBI'});
}

if ($gff or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_gff.pl  -wsversion ${ws_release} -bgijson ${bgi_file} " . 
	"-gffin ${build_home}/DATABASES/${ws_release}/SEQUENCES/elegans.processed.gff3.gz " .
	"-gffout ${gff_file}"
	);
    $datatypes_processed{'GFF'} = process_datatype('GFF', \@cmds, $sub_dir, $log, $DEFAULT_MEM, '01:00:00', $datatypes_processed{'BGI'});
}

if ($allele or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_allele_json.pl -database ${build_home}/DATABASES/${ws_release} " .
	"-wsversion ${ws_release} -bgijson ${bgi_file} -outfile ${allele_file}",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/allele/alleleMetaData.json " .
	"-d ${allele_file}"
	);
    $datatypes_processed{'ALLELE'} = process_datatype('allele', \@cmds, $sub_dir, $log, $DEFAULT_MEM, '1:00:00', $datatypes_processed{'BGI'});
}

if ($phenotype or $all) {
    my @cmds = ("perl ${cvs_dir}/AGR/get_chebi_name_map.pl -o $chebi_map_file") unless -e $chebi_map_file;
    push @cmds, (
	"perl ${cvs_dir}/AGR/make_agr_phenotype_json.pl -database ${build_home}/DATABASES/${ws_release} " .
	"-wsversion ${ws_release} -bgijson ${bgi_file} -outfile ${sub_dir}/WB_${agr_schema}_phenotype.json -chebi ${chebi_map_file}",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/phenotype/phenotypeMetaDataDefinition.json " .
	"-d ${sub_dir}/WB_${agr_schema}_phenotype.json"
	);
    $datatypes_processed{'PHENOTYPE'} = process_datatype('phenotype', \@cmds, $sub_dir, $log, $DEFAULT_MEM, '1:00:00',
							 $datatypes_processed{'CHEBI'}, $datatypes_processed{'BGI'});
}

if ($expression or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_expression_json.pl -database ${build_home}/DATABASES/${ws_release} " .
	"-wsversion ${ws_release} -bgijson ${bgi_file}  -wb2uberon ${agr_resources}/wormbase_to_uberon.agr_2_2.txt " .
	"-outfile ${sub_dir}/WB_${agr_schema}_expression.json",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/expression/wildtypeExpressionMetaDataDefinition.json " .
	"-d ${sub_dir}/WB_${agr_schema}_expression.json"
	);
    $datatypes_processed{'EXPRESSION'} = process_datatype('expression', \@cmds, $sub_dir, $log, $DEFAULT_MEM, '1:00:00', $datatypes_processed{'BGI'});
}

if ($ltp_variations or $all) {
    my @cmds = (
	"ruby ${cvs_dir}/AGR/make_agr_variations.rb -f ${fasta_file} " .
	"-g ${build_home}/DATABASES/${ws_release}/SEQUENCES/elegans.processed.gff3.gz -d ${build_home}/DATABASES/${ws_release} " .
	"-w ${ws_release} -q ${cvs_dir}/AGR/agr_variations.def -o ${sub_dir}/WB_${agr_schema}_variations.json",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/allele/variantMetaData.json " .
	"-d ${sub_dir}/WB_${agr_schema}_variations.json" 
	);
    $datatypes_processed{'VARIATION'} = process_datatype('LTP variation', \@cmds, $sub_dir, $log, $DEFAULT_MEM, '1:00:00', $datatypes_processed{'FASTA'});
}

if ($htp_variations or $all) {
    my @cmds = (
	"ruby ${cvs_dir}/AGR/make_agr_variations.rb -a -f ${fasta_file} " .
	"-g ${build_home}/DATABASES/${ws_release}/SEQUENCES/elegans.processed.gff3.gz -d ${build_home}/DATABASES/${ws_release} " .
	"-w ${ws_release} -q ${cvs_dir}/AGR/agr_variations.def -o ${sub_dir}/WB_${agr_schema}_htp_variations.json",
	"python3 ${cvs_dir}/AGR/agr_variations_json2vcf.py -j ${sub_dir}/WB_${agr_schema}_htp_variations.json -g ${gff_file} " .
	"-m WB -o ${sub_dir}/WB_${agr_schema}_htp_variations.vcf -f ${fasta_file} -w"
	);
    $datatypes_processed{'HTVCF'} = process_datatype('HTP variation', \@cmds, $sub_dir, $log, $HIGH_MEM, '2:00:00',
						    $datatypes_processed{'BGI'}, $datatypes_processed{'GFF'},
						    $datatypes_processed{'FASTA'});
}

if ($agm or $all) {
    my @cmds = (
	"perl  ${cvs_dir}/AGR/make_agr_AGM.pl -alleles ${allele_file} -database ${build_home}/DATABASES/${ws_release} " .
	"-wsversion ${ws_release} -outfile ${sub_dir}/WB_${agr_schema}_AGM.json",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/affectedGenomicModel/affectedGenomicModelMetaData.json " .
	"-d ${sub_dir}/WB_${agr_schema}_AGM.json"
	);
    $datatypes_processed{'AGM'} = process_datatype('AGM', \@cmds, $sub_dir, $log, $DEFAULT_MEM, '1:00:00', $datatypes_processed{'BGI'},
						   $datatypes_processed{'ALLELE'});
}

if ($construct or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_construct_json.pl -database ${build_home}/DATABASES/${ws_release} " .
	"-wsversion ${ws_release} -outfile ${sub_dir}/WB_${agr_schema}_construct.json",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/construct/constructMetaData.json " .
	"-d ${sub_dir}/WB_${agr_schema}_construct.json"
	);
    $datatypes_processed{'CONSTRUCT'} = process_datatype('construct', \@cmds, $sub_dir, $log, $DEFAULT_MEM, '1:00:00');
}

if ($interactions or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_interactions_psimitab.pl  -database ${build_home}/DATABASES/${ws_release} " .
	"-outfile ${sub_dir}/WB_${agr_schema}_interactions.psi-mi-tab"
	);
    $datatypes_processed{'INTERACTION-SOURCE_MOL'} = process_datatype('molecular interactions', \@cmds, $sub_dir, $log,
								      $DEFAULT_MEM, '1:00:00');
}

if ($genetic_interactions or $all) {
    my @cmds = (
	"perl  ${cvs_dir}/AGR/make_agr_genetic_interaction_psimitab.pl -database ${build_home}/DATABASES/${ws_release} " . 
	"-outfile ${sub_dir}/WB_${agr_schema}_genetic_interactions.psi-mi-tab"
	);
    $datatypes_processed{'INTERACTION-SOURCE_GEN'} = process_datatype('genetic interactions', \@cmds, $sub_dir, $log,
								      $DEFAULT_MEM, '1:00:00');
}

if ($hts or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_HTS.pl  -database ${build_home}/DATABASES/${ws_release} " .
	"-sample ${sub_dir}/WB_${agr_schema}_sample.json -dataset ${sub_dir}/WB_${agr_schema}_dataset.json " .
	"-wsversion ${ws_release} -wb2uberon ${agr_resources}/wormbase_to_uberon.agr_2_2.txt",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/htp/datasetSample/datasetSampleMetaDataDefinition.json " .
	"-d ${sub_dir}/WB_${agr_schema}_sample.json",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/htp/dataset/datasetMetaDataDefinition.json " .
	"-d ${agr_uploads}/${agr_release}/${upload_date}/WB_${agr_schema}_dataset.json"
	);
    $datatypes_processed{'HTPDATASAMPLE'} = $datatypes_processed{'HTPDATASET'}  = process_datatype('HTS dataset and sample', \@cmds,
												   $sub_dir, $log, $DEFAULT_MEM, '1:00:00');
}

# Reference files now submitted directly by Caltech

#if ($reference or $all) {
#    my @cmds = (
#	"perl ${cvs_dir}/AGR/make_agr_reference_json.pl -database ${build_home}/DATABASES/${ws_release} " .
#	"-wsversion ${ws_release} -ref ${sub_dir}/WB_${agr_schema}_reference.json " .
#	"-refex ${sub_dir}/WB_${agr_schema}_reference_exchange.json",
#	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/resourcesAndReferences/referenceMetaData.json " .
#	"-d ${sub_dir}/WB_${agr_schema}_reference.json",
#	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/resourcesAndReferences/referenceExchangeMetaData.json " .
#	"-d ${sub_dir}/WB_${agr_schema}_reference_exchange.json"
#	);
 #   $datatypes_processed{'REFERENCE'} = $datatypes_processed{'REF-EXCHANGE'} = process_datatype('reference', \@cmds, $sub_dir, $log,
#												$lsf_manager, $HIGH_MEM);
#}

if ($molecule or $all) {
    my @cmds = (
	"perl ${cvs_dir}/AGR/make_agr_molecules.pl -database ${build_home}/DATABASES/${ws_release} -wsversion ${ws_release} " .
	"-outfile ${sub_dir}/WB_${agr_schema}_molecule.json",
	"${agr_schema_repo}/bin/agr_validate.py -s ${agr_schema_repo}/ingest/molecules/moleculeMetaData.json " .
	"-d ${sub_dir}/WB_${agr_schema}_molecule.json"
	);
    $datatypes_processed{'MOLECULE'} = process_datatype('molecule', \@cmds, $sub_dir, $log, $DEFAULT_MEM, '1:00:00');
}

my @job_ids_to_check;
for my $datatype (keys %datatypes_processed) {
    push @job_ids_to_check, $datatypes_processed{$datatype} unless $datatypes_processed{$datatype} == 0;
}

WormSlurm::wait_for_jobs(@job_ids_to_check);
check_exit_statuses(\@job_ids_to_check, $log);

submit_data(\%datatypes_processed, \%files_to_submit, $log) if !$test;
cleanup_unzipped($sub_dir);

$log->mail;

exit(0);


sub check_exit_statuses{
    my ($job_ids, $log) = @_;

    for my $job_id (@$job_ids) {
	my $exit_code = WormSlurm::get_exit_code($job_id);
	if ($exit_code == 0) {
	    $log->write_to("Slurm job $job_id completed successfully\n");
	}
	else {
	    $log->error("ERROR: Slurm job $job_id exited with exit code $exit_code\n");
	}
    }

    return;
}


sub cleanup_unzipped {
    my $sub_dir = shift;

    my @files = glob($sub_dir . '/*');
    for my $file (@files) {
	unlink $file if -e $file . '.gz';
    }

    return;
}


sub process_datatype {
    my ($datatype, $cmds, $sub_dir, $log, $mem, $time, @dependencies) = @_;

    my $datatype_cmd_nr = 0;
    my $slurm_job_id;
	
    for my $cmd(@$cmds) {
	$datatype_cmd_nr++;
	my $job_name = $datatype . '_' . $datatype_cmd_nr;
	$job_name =~ s/\s/_/g;
	my $outfile = "${sub_dir}/slurm_logs/${job_name}.out";
	my $errfile = "${sub_dir}/slurm_logs/${job_name}.err";

	$time = '00:05:00' if $cmd =~ /agr_validate/;
	if (@dependencies) {
	    $slurm_job_id = WormSlurm::submit_job_with_name_and_dependencies($cmd, 'production', $mem, $time, $errfile, $outfile, $job_name, \@dependencies);
	} else {
	    $slurm_job_id = WormSlurm::submit_job_with_name($cmd, 'production', $mem, $time, $errfile, $outfile, $job_name);
	}

	if ($slurm_job_id) {
	    $log->write_to("$cmd submitted to Slurm with job ID " . $slurm_job_id . "\n\n");
	}
	else {
	    $log->log_and_die("Slurm submission of $datatype data failed: $?, $@\n");
	}
	@dependencies = ($slurm_job_id);
    }
    
    
    return $slurm_job_id;
}


sub submit_data {
    my ($datatypes_processed, $files_to_submit, $log) = @_;


    for my $type (@DATATYPES) {	
	next unless $datatypes_processed->{$type};
	my $file = $files_to_submit->{$type};

	my $fms_subtype = '';
	my $fms_datatype = $type;
	if ($type =~ /^(.+)_(.+)$/) {
	    $fms_datatype = $1;
	    $fms_subtype = $2;
	}

	
	my $exit_code = system("gzip -9 -c $file > $file.gz");
	if ($exit_code != 0) {
	    $log->error("Compression of $file failed with exit code $exit_code, $fms_datatype submission not completed\n\n");
	    return;
	}

	my $fms_file_link;
	if ($fms_datatype eq 'FASTA') {
	    $fms_file_link = '_WBcel235=@';
	} else {
	    $fms_file_link = $fms_subtype eq '' ? '_WB=@' : $fms_datatype . '_WB-' . $fms_subtype . '=@';
	}

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
    }

    return;
}    
