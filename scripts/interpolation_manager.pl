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
use Modules::WormSlurm;

##############################
# Script variables (run)     #
##############################

my ( $help, $debug, $test, $store, $wb, $flags, $noload, $pseudo, $fix );

##############################
# command-line options       #
##############################

GetOptions(
  "help"    => \$help,
  "debug=s" => \$debug,
  "test"    => \$test,
  "noload"  => \$noload,
  "store:s" => \$store,
  "pseudo"  => \$pseudo,
  "fix"     => \$fix,
);

$flags = '';

if ($store) { 
  $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n"); 
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

my $acedir = $wb->acefiles;
my $fix_acefile    = "$acedir/genetic_map_fixes.ace";
my $pseudo_acefile = "$acedir/pseudo_map_positions.ace";

my $prep_flags = ($fix) ? "$flags -preparefix -fixacefile $fix_acefile" : "$flags -preparenofix";
my $cmd = $wb->build_cmd_line("interpolate_gff.pl $prep_flags",$store);

WormSlurm::submit_job_and_wait($cmd, 'production', '4000m', '1-00:00:00', '/dev/null', '/dev/null');

my %slurm_jobs;
foreach my $i ($wb->get_chromosome_names(-prefix => 1)) {
    my $cmd2 = $wb->build_cmd_line("interpolate_gff.pl -chrom $i -all $flags",$store);
    my $job_id = WormSlurm::submit_job($cmd2, 'production', '4000m', '1-00:00:00', '/dev/null', '/dev/null');
    $slurm_jobs{$job_id} = $cmd2;
}
WormSlurm::wait_for_jobs(keys %slurm_jobs);

for my $job_id (keys %slurm_jobs) {
  $log->write_to("Slurm job $job_id (" . $slurm_jobs{$job_id} . ") exited non zero\n") if WormSlurm::get_exit_code($job_id) != 0;
}

if ($pseudo) {
  $wb->run_script("make_pseudo_map_positions.pl -acefile $pseudo_acefile -noload", $log);
}

#################
unless ( $noload ) {
  foreach my $file ( glob("$acedir/interpolated_*_*.ace") ) {
    $wb->load_to_database( $wb->autoace, $file, "interpolated_positions", $log );
  }

  if ($fix and -e $fix_acefile) {
    $wb->load_to_database($wb->autoace,
                          $fix_acefile,
                          "genetic_map_corrections",
                          $log, 
                          0, # no backup
                          1 #accept large differences
        );
  }
  if ($pseudo and -e $pseudo_acefile) {
    $wb->load_to_database($wb->autoace,
                          $pseudo_acefile,
                          "pseudo_map_posn",
                          $log, 
                          0, # no backup
                          1 #accept large differences
        );
  }
}

$log->mail();

exit(0);
