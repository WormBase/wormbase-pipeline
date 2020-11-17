#!/usr/bin/perl
use strict;

use Getopt::Long;

my ($password, $species, $database, $fasta, $gff_splits, $test, $no_load, $debug);

$species = 'Elegans';
GetOptions(
    "password=s"   => \$password,
    "species=s"    => \$species,
    "database=s"   => \$database,
    "gff_splits=s" => \$gff_splits,
    "fasta=s"      => \$fasta,
    "test"         => \$test,
    "noload"       => \$no_load,
    "debug"        => \$debug,
    ) or print_usage();

$species = ucfirst(lc($species));
die "Species must be Elegans or Briggsae\n" unless $species eq 'Elegans' or $species eq 'Briggsae';

my $url = 'mysql://' . $ENV{'WORM_DBUSER'} . ':' . $password . '@' . $ENV{'WORM_DBHOST'} . ':' . 
    $ENV{'WORM_DBPORT'} . '/wb_vep_' . lc($species) . '_ehive';
my $load_data = $no_load ? 0 : 1;
my $test_run = $test ? 1 : 0;
my $debug_mode = $debug ? 1 : 0;

if (!defined $database) {
    my $basedir = $test ? $ENV{'WORMPUB'} . "/TEST/BUILD" : $ENV{'WORMPUB'} . "/BUILD";
    $database = $species eq 'Elegans' ? "$basedir/autoace" : "$basedir/" . lc($species);
}

if (!defined $gff_splits) {
    $gff_splits = "$database/GFF_SPLITS";
}

if (!defined $fasta) {
    $fasta = "$database/SEQUENCES/" . lc($species) . '.genome.fa';
}

my $init_cmd = "init_pipeline.pl ProteinConsequence::ProteinConsequence_conf -password $password -ws_release $ws_release " . 
    "-species $species -database $database -fasta $fasta -gff_dir $gff_splits -load_to_ace $load_data -test $test_run " .
    "-debug $debug_mode";

for my $cmd ($init_cmd, "setenv EHIVE_URL \"$url\"", "beekeeper.pl -url $url -loop") {
    system($cmd);
};


sub print_usage{
print  <<USAGE;
run_vep_pipeline.pl options:
    -password PASSWORD      password for eHive pipeline database
    -species SPECIES        Elegans or Briggsae
    -database DATABASE_DIR  specify a different db for testing
    -gff_splits DIR NAME    specify a different dir containing split GFF files for testing
    -fasta FILE_NAME        specify a different FASTA file containing genome sequence for testing
    -test                   use the test database
    -noload                 do no update AceDB
    -debug                  add debugging messages
USAGE

exit 1;
}
