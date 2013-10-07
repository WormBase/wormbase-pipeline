#!/usr/bin/env perl
#####################################################
# put interpolate snps jobs on queue and waits for them to exit
#
####################################################
use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;

##############################
# Script variables (run)     #
##############################

my ( $help, $debug, $test, $store, $wb, $flags, $noload, $nopseudo );

##############################
# command-line options       #
##############################

GetOptions(
    "help"    => \$help,
    "debug=s" => \$debug,
    "test"    => \$test,
    "noload"  => \$noload,
    "store:s" => \$store,
    "nopseudo"=> \$nopseudo
);


if ($store) { 
  $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n"); 
  $flags = "-store";
}
else { 
  $wb = Wormbase->new(-debug => $debug, 
                      -test => $test,
      );
  $flags .= "-debug $debug " if $debug;
  $flags .= "-test "         if $test;
}

my $log = Log_files->make_build_log($wb);

####################################

my $m      = LSF::JobManager->new();
my @bsub_opts = (-M => 4000,
                 -R => 'select[mem>=4000] rusage[mem=4000]',
                 -o => '/dev/null',
    );
my $cmd = $wb->build_cmd_line("interpolate_gff.pl -prep $flags");
my $mother = $m->submit(@bsub_opts, $cmd);
my $myid   = $mother->id;

push @bsub_opts, (-w => "ended($myid)");
foreach my $i ($wb->get_chromosome_names(-prefix => 1)) {
  my $cmd2 = $wb->build_cmd_line("interpolate_gff.pl -chrom $i -all $flags");
  $m->submit( @bsub_opts, $cmd2);
}

$m->wait_all_children( history => 1 );
print "All children have completed!\n";

############################
for my $job ( $m->jobs ) {    # much quicker if history is pre-cached
  $log->write_to("$job exited non zero\n") if $job->history->exit_status != 0;
}
$m->clear;                    # clear out the job manager to reuse.

#################
unless ( $noload ) {
  my $acedir = $wb->acefiles;
  foreach my $file ( glob("$acedir/interpolated_*_*.ace") ) {
    $wb->load_to_database( $wb->autoace, $file, "interpolated_positions", $log );
  }

  my $gmap_fixes =  "$acedir/genetic_map_fixes.ace";
  if (-e $gmap_fixes and -s $gmap_fixes) {
    $wb->load_to_database($wb->autoace,
                          "$acedir/genetic_map_fixes.ace",
                          "genetic_map_corrections",
                          $log, 
                          0, # no backup
                          1 #accept large differences
        );
  }
}

$wb->run_script("make_pseudo_map_positions.pl", $log) unless $nopseudo;
$log->mail();

exit(0);
