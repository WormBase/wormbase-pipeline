#!/usr/bin/env perl 

use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use DBI;
use strict;
use Modules::WormSlurm;
use File::Which qw(which);

my ($debug,
    $store,
    $test,
    $wb,
    $rdb_user,
    $rdb_pass,
    $rdb_host, 
    $rdb_port,
    $rdb_name,
    $aligner_dir,
    $aligner_exe,
    @run_species,
    $run_clustal,
    $dump_clustal,
    $ace_database,
    $dontclean);

GetOptions(
  'debug=s'          => \$debug,
  'store=s'          => \$store,
  'runspecies=s@'    => \@run_species,
  'test=s'           => \$test,
  'host=s'           => \$rdb_host,
  'port=s'           => \$rdb_port,
  'user=s'           => \$rdb_user,
  'pass=s'           => \$rdb_pass,
  'acedb=s'          => \$ace_database,
  'dontclean'        => \$dontclean,
  'alignerdir'       => \$aligner_dir,
  'aligneerexe'      => \$aligner_exe,
  'run'              => \$run_clustal,
  'dump'             => \$dump_clustal,
) ||die(@!);

if ($store) {
  $wb = Storable::retrieve($store) or die('cannot load from storable');
}
else { 
  $wb = Wormbase->new( 
    -debug => $debug, 
    -test => $test);
}

my $log = Log_files->make_build_log($wb);


if (not $run_clustal and not $dump_clustal) {
  $run_clustal = $dump_clustal = 1;
}

$rdb_name = "worm_clw";
$rdb_host = "ebiworm-db" if not defined $rdb_host;
$rdb_port = "3478" if not defined $rdb_port;
$rdb_user = "wormadmin" if not defined $rdb_user;
$rdb_pass = "worms" if not defined $rdb_pass;
$ace_database = $wb->autoace if not defined $ace_database;


$aligner_exe = "muscle" if not defined $aligner_exe;
if (not defined $aligner_dir) {
  $aligner_dir = which($aligner_exe);
  if (not defined $aligner_dir) {
    $log->die("Could not find '$aligner_exe' on PATH. Exiting\n");
  }
  $aligner_dir =~ s/$aligner_exe$//;
}

my $run_jobs_failed = 0;
if ($run_clustal) {
  my %accessors = $wb->species_accessors;
  $accessors{elegans} = $wb;
  
  my $dbconn = DBI->connect("DBI:mysql:dbname=${rdb_name};host=${rdb_host};port=${rdb_port}" ,$rdb_user,$rdb_pass) or $log->log_and_die($DBI::errstr);
  
  my $scratch_dir = $wb->build_lsfout;
  my %job_info;

  my %slurm_jobs;
  foreach my $species (sort keys %accessors) {
      # if no runspecies were supplied at command line, we use the 
      # database versions to work out which species to run

      $log->write_to("Considering running $species...\n");
      if (not @run_species) {
	  if ($accessors{$species}->version != $wb->version) {
	      $log->write_to(" Looks like $species was not (re)built. Will not re-run alignments\n");
	      next;
	  }
      } else {
	  if (not grep { $species eq $_ } @run_species) {
	      $log->write_to(" Skipping $species because explicit list defined, and not part of that list\n");
	      next;
	  } 
      }
      
      $log->write_to(" Running $species...\n");
      
      my $infile = sprintf("%s/%spep%s", 
			   $accessors{$species}->wormpep,
			   $accessors{$species}->pepdir_prefix,
			   $accessors{$species}->version);
      
      $log->log_and_die("Could not find $infile\n") if not -e $infile;
      $log->write_to("Will repopulate for $species using $infile\n");
      
      my $prefix = $accessors{$species}->pep_prefix;        
      $dbconn->do("DELETE FROM clustal WHERE peptide_id LIKE \'$prefix\%\'") unless $dontclean;
      
      my $cmd_prefix = "$scratch_dir/clustal.$species.$$.";
      
      my $batch_total = 20;
      for(my $batch_idx = 1; $batch_idx <= $batch_total; $batch_idx++) {
	  my $cmd_file = "${cmd_prefix}.${batch_idx}.cmd.csh";
	  my $cmd_out  = "${cmd_prefix}.${batch_idx}.slurmout";
	  my $cmd_err  = "${cmd_prefix}.${batch_idx}.slurmerr";
	  my $job_name = "worm_clustal";
	  
	  my $cmd = "clustal_runner.pl" 
	      . " -batchid $batch_idx"
	      . " -batchtotal $batch_total"
	      . " -user $rdb_user"
	      . " -pass $rdb_pass"
	      . " -host $rdb_host"
	      . " -port $rdb_port"
	      . " -dbname $rdb_name"
	      . " -pepfile $infile"
	      . " -database $ace_database"
	      . " -alignerdir $aligner_dir";
	  $cmd = $accessors{$species}->build_cmd($cmd);
	  
	  # it is useful to still make the command files so that if any batch jobs fail, the command files are available to be run by hand
	  open(my $cmd_fh, ">$cmd_file") or $log->log_and_die("Could not open $cmd_file for writing\n");
	  print $cmd_fh "#!/bin/csh\n";
	  print $cmd_fh "$cmd\n";
	  close($cmd_fh);
	  chmod 0777, $cmd_file;
	  
	  my $job_id = WormSlurm::submit_job_with_name($cmd, 'production', '6g', '2:00:00', $cmd_out, $cmd_err, $job_name);
	  $slurm_jobs{$job_id} = $cmd;
	  if (defined $job_id) {
	      $job_info{$job_id} = {
		  output_file => $cmd_out,
		  command_file => $cmd,
	      }
	  } else {
	      $log->log_and_die("Could not submit job $cmd\n");
	  }
      }
  }
  WormSlurm::wait_for_jobs(keys %slurm_jobs);
  
  $log->write_to("All clustal jobs have completed!\n");
  for my $job_id (keys %slurm_jobs) {    # much quicker if history is pre-cached
    my $job_cmd =  $job_info{$job_id}->{command_file};
    my $job_out =  $job_info{$job_id}->{output_file};
    
    if (WormSlurm::get_exit_code($job_id) == 0) {
	unlink $job_out;   
    } else {    
	$log->error("$job exited non zero: check $job_out and re-run '$job_cmd' before dumping\n");
	$run_jobs_failed++;
    }
  }
}

if ($run_jobs_failed) {
  $log->log_and_die("Some clustal jobs failed; you will need to re-run these manually before dumping\n");
}


if ($dump_clustal) {
  my $outdir = $wb->misc_output;
  my $outfile = "$outdir/wormpep_clw.sql.gz";
  my $cmd = "mysqldump -u $rdb_user -p$rdb_pass -h $rdb_host -P $rdb_port $rdb_name | gzip > $outfile";
  $wb->run_command($cmd, $log);
}

$log->mail;
exit(0);
