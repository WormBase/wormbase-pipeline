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
my @bsub_opts = (-M => 4000,
                 -R => 'select[mem>=4000] rusage[mem=4000]',
                 -o => '/dev/null');
my $lsf = LSF::JobManager->new();

my $acedir = $wb->acefiles;
my $fix_acefile    = "$acedir/genetic_map_fixes.ace";
my $pseudo_acefile = "$acedir/pseudo_map_positions.ace";

my $prep_flags = ($fix) ? "$flags -preparefix -fixacefile $fix_acefile" : "$flags -preparenofix";
my $cmd = $wb->build_cmd_line("interpolate_gff.pl $prep_flags",$store);

my $mother = $lsf->submit(@bsub_opts, $cmd);
my $myid   = $mother->id;

push @bsub_opts, (-w => "ended($myid)");
foreach my $i ($wb->get_chromosome_names(-prefix => 1)) {
  my $cmd2 = $wb->build_cmd_line("interpolate_gff.pl -chrom $i -all $flags",$store);
  $lsf->submit( @bsub_opts, $cmd2);
}

$lsf->wait_all_children( history => 1 );

for my $job ( $lsf->jobs ) {
  $log->write_to("$job exited non zero\n") if $job->history->exit_status != 0;
}
$lsf->clear;                    # clear out the job manager to reuse.

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
