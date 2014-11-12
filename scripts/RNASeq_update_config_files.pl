#!/usr/bin/env perl
#
# test_RNASeq.pl
# 
# by Gary Williams                        
#
# This searches the ENA SRA for Studies from the current species
# It updates the MISC_DYNAMIC/SHORT_READS/species/Studies.ini files
# withe details of new Studies and then makes an INI file of the details 
# of Experiments in that Study.
#
# This will miss Experiments that have a missing Study entry or which have
# a Study that has a missing or inappropriate species field.
# To catch these, a search of the ENA SRA for all Experiments in the
# current species is then made and any new Experiments have their Study details
# stored and are added to their Study file.
#
# This appears to catch all new details.
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2014-11-12 15:38:21 $      

use strict;
use lib $ENV{'CVS_DIR'};
use Carp;
use Modules::RNASeq;

use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species);
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
            "verbose"    => \$verbose,
            "store:s"    => \$store,
            "species:s"  => \$species, # the default is elegans
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

my %data;
my $data;

my $RNASeq = RNASeq->new($wormbase, $log);

# get studies and then get the experiment details

#my @studies_changed = $RNASeq->update_study_config_data;

# we may have missed some experiments if we didn't find their study in
# the first search, so do the search the other way round as well -
# find Experiments and then get their Study details

my @new_experiments = $RNASeq->update_experiment_config_data;

#$log->write_to("\nStudies changed\n");
#foreach my $study (@studies_changed) {
#  $log->write_to("$study\n");
#}

$log->write_to("\nNew experiments\n");
foreach my $experiment (@new_experiments) {
  $log->write_to("$experiment\n");
}

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);
