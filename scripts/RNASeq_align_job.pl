#!/usr/bin/env perl
#
# RNASeq_align_job.pl
# 
# by Gary Williams                        
#
# Runs one alignment and then runs cufflinks and gets the junctions ace file.
# This script is usually run under LSF.
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2014-02-06 14:04:00 $      

use strict;
use lib $ENV{'CVS_DIR'};
use Carp;
use Modules::RNASeq;

use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species, $new_genome, $check, $experiment_acession, $database, $threads, $notbuild, $test_set);
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
	    "test_set"   => \$test_set, # use the limited set of experiments in Studies_test.ini
            "verbose"    => \$verbose,
            "store:s"    => \$store,
            "species:s"  => \$species, # the default is elegans
	    "new_genome" => \$new_genome, # redo all of the alignments against a fresh copy of the genome
	    "check"      => \$check, # test to see if any cufflinks etc. results are missing
	    "expt:s"     => \$experiment_acession, # the experiment accession to align
	    "database:s" => \$database,
	    "threads:i"  => \$threads, # number of threads to use when running STAR
	    "notbuild"   => \$notbuild, # don't try to run cufflinks (for when it is run before doing the Build) 
);


$species = 'elegans' unless $species;

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
                             );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


######################################
# variables and command-line options # 
######################################

my $data;

if (!defined $database) {$database = $wormbase->autoace}

my $RNASeq = RNASeq->new($wormbase, $log, $new_genome, $check, $test_set);

$log->write_to("Get experiment details from config\n");
my $data = $RNASeq->get_transcribed_long_experiments();

$data = $data->{$experiment_acession};
    
my $accession=$data->{'experiment_accession'};
my $library_strategy=$data->{'library_strategy'};
my $library_selection=$data->{'library_selection'};
my $library_layout=$data->{'library_layout'};
my $library_source=$data->{'library_source'};
$log->write_to("Running alignment job with experiment=$accession $library_source $library_strategy $library_selection $library_layout\n");

$RNASeq->align_star($data, $threads, $notbuild);



$log->mail();
print "Finished.\n" if ($verbose);
exit(0);
