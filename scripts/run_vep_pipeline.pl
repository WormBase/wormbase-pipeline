#!/usr/bin/perl
use strict;

use Getopt::Long;
use Log_files;
use Wormbase;

my ($password, $species, $database, $fasta, $gff_splits, $log_dir, $ace_dir, $test, $load, $debug, $help);

$species = 'Elegans';
GetOptions(
    "password=s"   => \$password,
    "species=s"    => \$species,
    "database=s"   => \$database,
    "gff_splits=s" => \$gff_splits,
    "fasta=s"      => \$fasta,
    "log_dir=s"    => \$log_dir,
    "ace_dir=s"    => \$ace_dir,
    "test"         => \$test,
    "load"         => \$load,
    "debug=s"      => \$debug,
    "help"         => \$help
    ) or print_usage();

print_usage() if $help;

$species = ucfirst(lc($species));
die "Species must be Elegans or Briggsae\n" unless $species eq 'Elegans' or $species eq 'Briggsae';

my $test_run = $test ? 1 : 0;
my $debug_mode = $debug ? $debug : 0;

if (!defined $database) {
    my $basedir = $test ? $ENV{'WORMPUB'} . "/TEST/BUILD" : $ENV{'WORMPUB'} . "/BUILD";
    $database = $species eq 'Elegans' ? "$basedir/autoace" : "$basedir/" . lc($species);
}

if (!defined $log_dir) {
    $log_dir = "$database/logs";
}

my $wb = Wormbase->new(
    -autoace  => $database,
    -organism => $species,
    -debug    => $debug,
    -test     => $test,
    );

$ace_dir = $wb->acefiles unless defined $ace_dir;


if ($load) {
    # Load collated file into AceDB
    my $log = Log_files->make_build_log($wb);
    
    my $out_file = $ace_dir . '/mapped_alleles.' . $wb->get_wormbase_version_name . '.ace';
    $wb->load_to_database($wb->autoace, $out_file, 'WB_VEP_pipeline', $log) 

}
else {
    # Initialise and run ProteinConsequence eHive pipeline
    if (!defined $gff_splits) {
	$gff_splits = "$database/GFF_SPLITS";
    }

    if (!defined $fasta) {
	$fasta = "$database/SEQUENCES/" . lc($species) . '.genome.fa';
    }

    my $url = 'mysql://' . $ENV{'WORM_DBUSER'} . ':' . $password . '@' . $ENV{'WORM_DBHOST'} . ':' . 
    $ENV{'WORM_DBPORT'} . '/wb_vep_' . lc($species) . '_ehive';

    $ENV{EHIVE_URL} = $url;

    my $init_cmd = "init_pipeline.pl ProteinConsequence::ProteinConsequence_conf -password $password " . 
	"-species $species -database $database -fasta $fasta -gff_dir $gff_splits " . 
	"-test $test_run -log_dir $log_dir -ace_dir $ace_dir -debug $debug_mode";
    system($init_cmd);
    system("beekeeper.pl -url $url -loop");
}


sub print_usage{
print  <<USAGE;
run_vep_pipeline.pl options:
    -password PASSWORD      password for eHive pipeline database
    -species SPECIES        Elegans or Briggsae
    -database DATABASE_DIR  specify a different db for testing
    -gff_splits DIR NAME    specify a different dir containing split GFF files for testing
    -fasta FILE_NAME        specify a different FASTA file containing genome sequence for testing
    -log_dir DIR NAME       specify a different directory for log files for testing
    -ace_dir DIR NAME       specify a different directory for final .ace file for testing
    -test                   use the test database
    -load                   load data into AceDB
    -help                   print this message
    -debug EMAIL            eamil for sending log messages to
USAGE

exit 1;
}
