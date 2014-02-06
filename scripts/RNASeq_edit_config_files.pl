#!/usr/bin/env perl
#
# test_RNASeq.pl
# 
# by Gary Williams                        
#
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2014-02-06 15:33:09 $      

use strict;
use lib $ENV{'CVS_DIR'};
use Carp;
use Modules::RNASeq;

use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species, $infile);
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
            "verbose"    => \$verbose,
            "store:s"    => \$store,
            "species:s"  => \$species, # the default is elegans
	    "infile:s"   => \$infile, # the file of editing commands
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

my $RNASeq = RNASeq->new($wormbase, $log, 0, 0);

my $experiments = $RNASeq->get_all_experiments();

# the input file should be in the form of
# condition thing_to_add
# where condition is in the form of 
# parameter=value
# e.g. experiment_accession=SRX037199
# and the thing to add in in the form of 
# parameter=value
# e.g. strandedness=unstranded

open (IN, "< $infile") || $log->log_and_die("Can't open $infile\n");
while (my $line = <IN>) {
  chomp $line;
  if ($line =~ /^#/) {next}
  if ($line eq '') {next}
  if ($line =~ /(.+?)=(.+?)\s+(.+?)=(.+)/) {
    my $param1 = $1;
    my $val1 = $2;
    my $param2 = $3;
    my $val2 = $4;

    # find all experiments satisfying the condition
    my @experiment_accessions;
    foreach my $expt (keys %{$experiments}) {
      foreach my $param (keys %{$experiments->{$expt}}) {
	if ($param1 eq $param && $val1 eq $experiments->{$expt}{$param}) {
	  push @experiment_accessions, $expt;
	}
      }
    }

    # find the ini file to use
    if ($#experiment_accessions >= 0) {
      foreach my $experiment_accession ( @experiment_accessions) {
	my $study_accession = $experiments->{$experiment_accession}{study_accession};
	my $experiment_ini = $RNASeq->read_experiments_from_study_config($study_accession);
	$experiment_ini->newval($experiment_accession, $param2, $val2);
	$experiment_ini->RewriteConfig;
	$log->write_to("SET $param2 = $val2 WHERE $param1 = $val1 (expt $experiment_accession) in Study $study_accession\n");
      }
    } else {
      $log->write_to("Can't find matching experiments for line: $line\n");
    }
  } else {
    $log->write_to("Unrecognised format: $line\n");
  }

}
close(IN);

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);
