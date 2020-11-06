#!/bin/env perl
# that should fire off 10 dump_promoter jobs with 16GB each and concatenate the results

use lib $ENV{CVS_DIR};
use LSF::JobManager;
use Storable;

my ($species,$store,$debug,$test,$database,$outfile,$chunk);
GetOptions(
     'species=s'  => \$species,
     'store=s'    => \$store,
     'debug=s'    => \$debug,
     'test'       => \$test,
     'database=s' => \$database,
     'outfile=s'  => \$outfile,
     'chunk=s'    => \$chunk,
)||die(@!);

my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")
}else {
  $wormbase = Wormbase->new( -debug => $debug, 
                             -test => $test,
                             -organism => $species);
}
$database = $wormbase->autoace if not defined $database;

my $log = Log_files->make_build_log($wormbase);

my $defaultMem = 16000;
my $lsf_out = $wormbase->build_lsfout;
my $lsf = LSF::JobManager->new();


for my $chunk (0..9){
  my $cmd = "preFTPdumps/perSpecies/dump_promoters.pl -outfile $outfile.$chunk -chunk $c";
  $cmd = $wb->build_cmd($cmd);

  my @lsf_opts = (-M => $mem, 
                  -R => "select[mem>=$mem] rusage[mem=$mem]",
                  -J => 'perSpeciesDumps', 
                  -o => "${lsf_out}/dump_promoters.pl.$chunk.lsfout");
  $lsf->submit(@lsf_opts, $cmd);
}

$log->write_to("Waiting for LSF jobs to finish.\n");
$lsf->wait_all_children( history => 1 );
for my $job ( $lsf->jobs ) {
  if ($job->history->exit_status != 0) {
    $log->write_to("WARNING: Job $job (" . $job->history->command . ") exited non zero: " . $job->history->exit_status . "\n");
  }
}
$lsf->clear;

# cleanup
`cat $outfile.* > $outfile`;
for my $chunk(0..9){
	unlink "$outfile.$chunk" if -e "$outfile.$chunk";
}

$log->mail
