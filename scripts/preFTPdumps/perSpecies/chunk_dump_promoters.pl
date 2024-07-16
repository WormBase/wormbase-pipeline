#!/bin/env perl
# that should fire off 10 dump_promoter.pl jobs with 16GB each and concatenate the results

use lib $ENV{CVS_DIR};
use Modules::WormSlurm;
use Storable;
use Getopt::Long;
use Wormbase;
use Log_files;

my ($species,$store,$debug,$test,$database,$outfile);
GetOptions(
     'species=s'  => \$species,
     'store=s'    => \$store,
     'debug=s'    => \$debug,
     'test'       => \$test,
     'database=s' => \$database,
     'outfile=s'  => \$outfile,
)||die(@!);

my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")
}else {
  $wormbase = Wormbase->new( -debug => $debug, 
                             -test => $test,
                             -organism => $species);
}
$database ||= $wormbase->autoace;
$outfile ||= $wormbase->reports . '/potential_promoters.fa';

my $log = Log_files->make_build_log($wormbase);

my $mem = '16g';
my $slurm_out = $wormbase->build_lsfout;

my %slurm_jobs;
for my $chunk (1..10){
    my $cmd = "preFTPdumps/perSpecies/dump_promoters.pl -outfile $outfile.$chunk -chunk $chunk -database $database";
    $cmd = $wormbase->build_cmd($cmd);
    ;
    my $job_id = WormSlurm::submit_job_with_name($cmd, 'production', $mem, '5:00:00', "$slurm_out/dump_promoters.pl.$chunk.slurm_out",
						 "$slurm_out/dump_promoters.pl.$chunk.slurm_err", 'perSpeciesPromoterDumps');
    $slurm_jobs{$job_id} = $cmd;
}

$log->write_to("Waiting for Slurm jobs to finish.\n");
WormSlurm::wait_for_jobs(keys %slurm_jobs);
for my $job_id (keys %slurm_jobs) {
    my $exit_code = WormSlurm::get_exit_code($job_id);
    if ($exit_code != 0) {
	$log->write_to("WARNING: Slurm job $job_id (" . $slurm_jobs{$job_id} . ") exited non zero: " . $exit_code . "\n");
    }
}

# cleanup
`cat $outfile.* > $outfile`;
for my $chunk(1..10){
	unlink "$outfile.$chunk" if -e "$outfile.$chunk";
}

$log->mail
