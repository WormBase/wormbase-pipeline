#!/usr/bin/env perl
#
# DESCRIPTION:
# script to submit blastx dumping scripts onto the farm
# and concatenate them at the end
#
# Last edited by: $Author: gw3 $
# Last edited on: $Date: 2015-05-18 10:06:10 $
#

my $usage = <<USAGE;
dump_blastx.pl options:
  -debug USER_NAME sets email address and debug mode
  -store FILE_NAME use a Storable wormbase configuration file
  -database DB_NAME Ensembl datbase to dump from. If not used deafults to the database of the storable
  -species SPECIES_NAME which species you want to dump (a.e. elegans, briggsae,....)
  -test if you want to use TEST_BUILD instead of BUILD
  -dumpdir DIRECTORY_NAME id you want to dump it to a different directory
  -rerun INT rerun jobs if output files (e.g. brugia_remapepx.50.ace, brugia_brugpepx.208.ace) failed and need to be repeated with '-rerun 50,208'
USAGE
 
use strict;
use Getopt::Long;

use lib $ENV{CVS_DIR};
use Wormbase;
use Coords_converter;
use Modules::WormSlurm;

my ($database,$store,$debug,$species,$test,$dumpdir, $slurm_mem, $rerun);
GetOptions(
	   'database=s' => \$database,
	   'store=s'    => \$store,
	   'debug=s'    => \$debug,
	   'species=s'  => \$species,
	   'test'       => \$test,
	   'dumpdir=s'  => \$dumpdir,
	   'slurmmem=s' => \$slurm_mem,
	   'rerun=s'    => \$rerun, 
	  ) || die($usage);

my $wormbase;

if ($store) {
  $wormbase = retrieve($store) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new(
			    -debug => $debug,
			    -test => $test,
			    -organism => $species,
			   );
}

my $log = Log_files->make_build_log($wormbase);

$species = $wormbase->species;
my @rerun = split /,/, $rerun if (defined $rerun);

# that might work until we change logic_names
my %logic2type = (
		  wormpepx => 'elegans',
		  ppapepx => 'pristionchus',
		  jappepx => 'japonica',
		  brugpepx => 'brugia',
		  brigpepx => 'briggsae',
		  remapepx => 'remanei',
		  brepepx => 'brenneri',
		  ovolpepx => 'ovolvulus',
                  srapepx => 'sratti',
                  tmupepx => 'tmuris',
		  gadflyx => '1',
		  ipi_humanx => '1',
		  yeastx => '1',
		  slimswissprotx => '1',
		 );

$slurm_mem = "4g" if not defined $slurm_mem;

$dumpdir ||= "$ENV{PIPELINE}/dumps";
$database ||= "worm_ensembl_${species}";

my $job_count;
my (@outfiles, @cmds);

foreach my $db (keys %logic2type){
  my @chroms = @{$wormbase->get_binned_chroms(-bin_size => 20)};
  $log->write_to("bsub commands . . . . \n\n");
  foreach my $chrom ( @chroms ) {
    $job_count++;
    my $outfile = "$dumpdir/${species}_$db.$job_count.ace";
    push @outfiles,$outfile;
    my $options="-database $database -logicname $db -outfile $outfile -sequence $chrom";
    $options .= ' -self' if $logic2type{$db} eq lc(ref $wormbase);

    my $cmd;
    if ($store) {
      $cmd = $wormbase->build_cmd_line("BLAST_scripts/blastx_dump.pl $options", $store);
    } else {
      $cmd = $wormbase->build_cmd("BLAST_scripts/blastx_dump.pl $options");
    }

    if (!defined $rerun || grep /^$job_count$/, @rerun) {
      push @cmds, $cmd;
    }
  }
}

my $slurm_jobs = WormSlurm::submit_jobs_and_wait_for_all(\@cmds, $ENV{LSB_DEFAULTQUEUE}, $slurm_mem, '4:00:00', '/dev/null', '/dev/null');
$log->write_to("All children have completed!\n");
for my $job_id (keys %$slurm_jobs) {
    my $exit_code = WormSlurm::get_exit_code($job_id);
    $log->error("Slurm job $job_id (" . $slurm_jobs->{$job_id} . ") exited non zero ($exit_code)\n") if $exit_code != 0;
}

# check that all files end with a blank line,
# otherwise the job that created them was probably terminated prematurely by LSF
$log->write_to("\nTesting output ace files to see if they completed\n");
foreach my $file (@outfiles) {
  $log->write_to("$file ");
  my $endline = `tail -1 $file`;
  if (-e $file && $endline =~ /^\s*\n$/) {
    $log->write_to("- looks OK\n");
  } else {
    $log->write_to("- appears to be prematurely terminated.\nlast line in $file is:\n$endline\n");
    $log->error;
  }
}

my $outfile="$dumpdir/${species}_blastx.ace";
$wormbase->run_command("rm -f $outfile", $log); # ensure we don't have a results file left over from previous Builds

if (! $log->report_errors ) {
  # concatenate the ace files into a big blob for later parsing with ensembl/ipi scripts
  $log->write_to("Concatenating the ace files to create $outfile\n");
  
  # in case of Elegans do something else
  if ($wormbase->species eq 'elegans' or
      $wormbase->species eq 'briggsae'){
    my @files = glob("$dumpdir/$species*x.*.ace");
    # the following will trigger the generation of SL_coords and clone_coords files if they do not exists
    Coords_converter->invoke($wormbase->orgdb, undef, $wormbase);
    open(my $out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");
    foreach my $file (@files ){
      $log->write_to("\tcat $file\n");
      open(my $conv_fh, "perl $ENV{CVS_DIR}/BLAST_scripts/convert_chromblast2clone.pl -species $species $file |")
	or $log->log_and_die("Could not open convert_chromblast2clone command\n");
      while(<$conv_fh>) {
	print $out_fh $_;
      }
    }
    close($out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
  } else {
    $wormbase->run_command("cat $dumpdir/$species*x*.ace > $outfile", $log);
  }
} else {
  $log->write_to("\n\nERROR FOUND WHEN PRODUCING OUTPUT FILES.\nNO RESULTS FILE WILL BE PRODUCED\n\n");
}

$log->mail();
