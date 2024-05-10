#!/usr/bin/env perl

use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use strict;
use Log_files;
use Storable;

my ($debug, $test, $database,$species, $verbose, $store );
my ($giface, $giface_server, $giface_client, $port);
my ($gff3, $gff, $dump_dir, $rerun_if_failed, $methods, $chrom_choice, $mem, $time);

my $dumpGFFscript = "GFF_method_dump.pl";

use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;

GetOptions (
    "debug:s"        => \$debug,
    "test"           => \$test,
    "verbose"        => \$verbose,
    "database:s"     => \$database,
    "dump_dir:s"     => \$dump_dir,
    "methods:s"      => \$methods,
    "chromosomes:s"  => \$chrom_choice,
    "store:s"        => \$store,
    "giface:s"       => \$giface,
    "gifaceserver:s" => \$giface_server, 
    "gifaceclient:s" => \$giface_client, 
    "gff3"           => \$gff3,
    "gff"            => \$gff,
    "rerunfail"      => \$rerun_if_failed,
    "species:s"      => \$species,
    "port:s"         => \$port,
    "mem:s"          => \$mem,
    "time:s"         => \$time,
    );
my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			     -organism => $species
			   );
  $store = $wormbase->autoace . "/".ref($wormbase).".store";
}

$species||=lc(ref($wormbase));

my $scratch_dir = $wormbase->build_lsfout;
my $host = qx('hostname');chomp $host;
$port = 23100 if not $port;

my @chroms = $wormbase->get_chromosome_names(-prefix => 1,-mito => 1);
my $chrom_lengths = $wormbase->get_chromosome_lengths(-prefix => 1,-mito => 1);
$wormbase->check_slurm;

if ($gff) {
    print "-gff is not a necessary command line option as the script defaulty to gff2 but this catches those who like editing the command lines instead of taking them from the build guide.\n";
}
my @methods     = split(/,/,join(',',$methods)) if $methods;
my @chromosomes = $chrom_choice ? split(/,/,join(',',$chrom_choice)):@chroms;

$database = $wormbase->autoace    unless $database;
if (not defined $dump_dir) {
  $dump_dir = (@methods) ? $wormbase->gff_splits : $wormbase->gff;
}

$giface = $wormbase->giface if not defined $giface;
$giface_server = $wormbase->giface_server if not defined $giface_server;
$giface_client = $wormbase->giface_client if not defined $giface_client;

my $log = Log_files->make_build_log($wormbase);
my $store_file = $wormbase->build_store;

$log->write_to("Dumping from DATABASE : $database\n\tto $dump_dir\n\n");
$log->write_to("\t chromosomes ".@chromosomes."\n");
$log->write_to("\tmethods ".@methods."\n\n") if @methods;
$log->write_to("\tno method specified\n\n") if not @methods;
$log->write_to("Slurm commands . . . . \n\n");

my (@individual_chrs, @batch_chrs);
if ($species eq 'elegans') {
  @individual_chrs = @chromosomes;
} else {
  foreach my $chr (@chromosomes) {
    $log->log_and_die("Could not find length for $chr - aborting\n") if not exists $chrom_lengths->{$chr};

    if ($chrom_lengths->{$chr} > 5000000) {
      push @individual_chrs, $chr;
    } else {
      push @batch_chrs, $chr;
    }
  }
}

if (@batch_chrs){
  &start_giface_server();
}

my $cmd_dir = $wormbase->autoace . "/TMP";
mkdir $cmd_dir, 0777;
my $cmd_file_root = "${cmd_dir}/dump_gff_batch_$$";
my $cmd_base = "$dumpGFFscript -database $database -species $species -dump_dir $dump_dir";
my $cmd_number = 0;
my $this_mem;
if (defined $mem){
    unless ($mem =~ /^\d+[m|g|M|G]$/) {
	$log->log_and_die("Expecting memory parameter in format /^\d+[m|g|M|G]$/\n");
    }
    $this_mem = $mem;
} else {
    $this_mem = '4500m';
}


unless (defined $time) {
    $time = '1:00:00';
}
unless ($time =~ /^(\d+\-[0-2]\d:[0-5]\d:[0-5]\d|[0-2]?\d:[0-5]\d:[0-5]\d)$/) {
    $log->log_and_die("Expecting time parameter in format D-HH:MM:SS\n");
}

my %submitted_jobs;
foreach my $chrom (@individual_chrs) {
  my $common_additional_params = "-chromosome $chrom -giface $giface";
  if ( @methods ) {
    foreach my $method ( @methods ) {
      my $this_cmd_num = ++$cmd_number;

      my $gff_out = sprintf("%s/%s.%s.gff%s", $dump_dir, $chrom, $method, ($gff3) ? "3" : "");
      my $slurm_out = "$scratch_dir/wormpubGFFdump.$chrom.$method.$this_cmd_num.slurmout";
      my $slurm_err = "$scratch_dir/wormpubGFFdump.$chrom.$method.$this_cmd_num.slurmerr";
      my $job_name = "worm_".$wormbase->species."_gffbatch.$this_cmd_num";

      my $cmd = "$cmd_base $common_additional_params";
      $cmd .= " -method $method";
      $cmd .= " -debug $debug" if $debug;
      $cmd .= " -gff3" if $gff3;
      $cmd = $wormbase->build_cmd_line($cmd, $store_file);      
      $log->write_to("Command: $cmd\n") if ($verbose);
      print "Command: $cmd\n" if ($verbose);

      # write out the command file to be executed
      my $cmd_file = "$cmd_file_root.$this_cmd_num";
      open (my $cmd_fh, ">$cmd_file") || die "Can't open file $cmd_file: $!";
      print $cmd_fh "#!/bin/csh\n";
      print $cmd_fh "$cmd\n";
      close ($cmd_fh);
      chmod 0777, $cmd_file;

      my $job_id = WormSlurm::submit_job_with_name($cmd_file, 'production', $this_mem, $time, $slurm_out, $slurm_err, $job_name);
      $submitted_jobs{$job_id} = $cmd_file;
    }
  } else {
    my $this_cmd_num = ++$cmd_number;

    my $gff_out = sprintf("%s/%s.gff%s", $dump_dir, $chrom, ($gff3) ? "3" : "");
    my $slurm_out = "$scratch_dir/wormpubGFFdump.$chrom.$this_cmd_num.lsfout";
    my $slurm_err = "$scratch_dir/wormpubGFFdump.$chrom.$this_cmd_num.lsferr";
    my $job_name = "worm_".$wormbase->species."_gffbatch.$this_cmd_num";

    my $cmd = "$cmd_base $common_additional_params";
    $cmd .= " -debug $debug" if $debug;
    $cmd .= " -gff3" if $gff3;
    $cmd = $wormbase->build_cmd_line($cmd, $store_file);
    $log->write_to("Command: $cmd\n") if ($verbose);
    print "Command: $cmd\n" if ($verbose);

    # write out the command file to be executed
    my $cmd_file = "$cmd_file_root.$this_cmd_num";
    open (my $cmd_fh, ">$cmd_file") || die "Can't open file $cmd_file: $!";
    print $cmd_fh "#!/bin/csh\n";
    print $cmd_fh "$cmd\n";
    close ($cmd_fh);
    chmod 0777, $cmd_file;

    my $job_id = WormSlurm::submit_job_with_name($cmd_file, 'production', $this_mem, $time, $slurm_out, $slurm_err, $job_name);
    $submitted_jobs{$job_id} = $cmd_file;
  }
}

if (@batch_chrs) {
  my $seq_list_file = "$cmd_dir/batch_seq_list.txt";
  open(my $batch_fh, ">$seq_list_file") or 
      $log->log_and_die("Could not open $seq_list_file for writing\n");
  foreach my $seq (@batch_chrs) {
    print $batch_fh "$seq\n";
  }
  close($batch_fh);
  my $this_mem;
  if (defined $mem){$this_mem = $mem} else {$this_mem = '100m'}

  my $common_additional_params = "-host $host -port $port -gifaceclient $giface_client -list $seq_list_file";
  
  if ( @methods ) {
    foreach my $method ( @methods ) {
      my $this_cmd_num = ++$cmd_number;

      my $gff_out = sprintf("%s/%s.gff%s", $dump_dir, $method, ($gff3) ? "3" : "");
      my $slurm_out = "$scratch_dir/wormpubGFFdump.$method.$this_cmd_num.slurmout";
      my $slurm_err = "$scratch_dir/wormpubGFFdump.$method.$this_cmd_num.slurmerr";
      my $job_name = "worm_".$wormbase->species."_gffbatch.$this_cmd_num";

      my $cmd = "$cmd_base $common_additional_params";
      $cmd .= " -method $method";
      $cmd .= " -debug $debug" if $debug;
      $cmd .= " -gff3" if $gff3;

      $cmd = $wormbase->build_cmd_line($cmd, $store_file);      
      $log->write_to("Command: $cmd\n") if ($verbose);
      print "Command: $cmd\n" if ($verbose);

      my $cmd_file = "$cmd_file_root.$this_cmd_num";
      open (my $cmd_fh, ">$cmd_file") || die "Can't open file $cmd_file: $!";
      print $cmd_fh "#!/bin/csh\n";
      print $cmd_fh "$cmd\n";
      close($cmd_fh);
      chmod 0777, $cmd_file;

      my $job_id = WormSlurm::submit_job_with_name($cmd_file, 'production', $this_mem, $time, $slurm_out, $slurm_err, $job_name);
      $submitted_jobs{$job_id} = $cmd_file;
    }
  }
  else {
    my $this_cmd_num = ++$cmd_number;

    my $gff_out = sprintf("%s/%s.gff%s", $dump_dir, $species, ($gff3) ? "3" : "");
    my $slurm_out = "$scratch_dir/wormpubGFFdump.$this_cmd_num.slurmout";
    my $slurm_err = "$scratch_dir/wormpubGFFdump.$this_cmd_num.slurmerr";
    my $job_name = "worm_".$wormbase->species."_gffbatch.$this_cmd_num";

    my $cmd = "$cmd_base $common_additional_params";
    $cmd .= " -debug $debug" if $debug;
    $cmd .= " -gff3" if $gff3;

    $cmd = $wormbase->build_cmd_line($cmd, $store_file);
    $log->write_to("Command: $cmd\n") if ($verbose);
    print "Command: $cmd\n" if ($verbose);

    # write out the command file to be executed
    my $cmd_file = "$cmd_file_root.$this_cmd_num";
    open (my $cmd_fh, ">$cmd_file") || die "Can't open file $cmd_file: $!";
    print $cmd_fh "#!/bin/csh\n";
    print $cmd_fh "$cmd\n";
    close ($cmd);
    chmod 0777, $cmd_file;

    my $job_id = WormSlurm::submit_job_with_name($cmd_file, 'production', $this_mem, $time, $slurm_out, $slurm_err, $job_name);
    $submitted_jobs{$job_id} = $cmd_file;
  }
}

WormSlurm::wait_for_jobs(keys %submitted_jobs);
$log->write_to("All GFF dump jobs have completed!\n");
my @problem_cmds;
for my $job_id (keys %submitted_jobs) {
    my $exit_code = WormSlurm::get_exit_code($job);
    if ($exit_codes == 0) {
	unlink $submitted_jobs{$job_id};
  } else {
    $log->write_to("Job $job_id (" . $submitted_jobs{$job_id} . ") exited non zero\n");
    push @problem_cmds, $submitted_jobs{$job_id};
  }
}


if (@problem_cmds and scalar(@problem_cmds) < 120 and $rerun_if_failed) { 
  ##################################################################
  # now try re-runnning any commands that failed
  if (@batch_chrs) {
    &stop_giface_server();
    &start_giface_server();
  }
  

  my $out = "$scratch_dir/wormpubGFFdump.rerun.slurmout";
  my $job_name = "worm_".$wormbase->species."_gffbatch";

  my $this_mem;
  if (defined $mem){$this_mem = $mem} else {$this_mem = '4500m'}
  
  my $rerun_count = 0;
  my @new_problem_cmds;
  my %resubmitted_jobs;
  foreach my $cmd_file (@problem_cmds) {
    $log->write_to("*** Attempting to re-run job: $cmd_file\n");
    my $slurm_out = sprintf("%s/wormpubGFFdump.rerun.slurmout", $scratch_dir, ++$rerun_count);
    my $slurm_err = sprintf("%s/wormpubGFFdump.rerun.slurmerr", $scratch_dir, ++$rerun_count);
    my $job_name = sprintf("worm_gff.%s.rerun_", $species, $rerun_count);
    
    my $job_id = WormSlurm::submit_job_with_name($cmd_file, 'production', $this_mem, $time, $slurm_out, $slurm_err, $job_name);
    $resubmitted_jobs{$job_id} = $cmd_file;
  }
  WormSlurm::wait_for_jobs(keys %resubmitted_jobs);
  
  my $failed = 0;
  for my $job_id ( keys %resubmitted_jobs ) {
      my $exit_code = WormSlurm::get_exit_code($job_id);
      if ($exit_code == 0) {
	  unlink $job->history->command;
      } else {
	  $failed++;
	  $log->error("\n\nERROR: Job $job_id (" . $resubmitted_jobs{$job_id} . ") exited non zero at the second attempt. See $slurm_out\n");
      push @new_problem_cmds, $resubmitted_jobs{$job_id};
    }
    $log->write_to("\n\nNumber of jobs that failed after the second attempt: $failed\n");
  }
  @problem_cmds = @new_problem_cmds;
  $lsf->clear;
}

if (@problem_cmds) {
  $log->error("ERROR: There are ". scalar @problem_cmds ." GFF_method_dump.pl Slurm jobs that have failed. See Slurm output files ($scratch_dir)\n");      
}


if (@batch_chrs and @individual_chrs) {
  my @to_delete;

  if (not $log->report_errors) {
    $log->write_to("Concatenating per-sequence dumps to all-sequence dump file(s)...\n");

    #certain steps require GFF files containing all sequences (for a particular method), so we aggergate the sequences here. 
    if (@methods) {
      foreach my $method (@methods) {
        my $target_gff = sprintf("%s/%s.gff%s", $dump_dir, $method, ($gff3) ? "3" : "");
        foreach my $seq (@individual_chrs) {
          my $source_gff = sprintf("%s/%s_%s.gff%s", $dump_dir, $seq, $method, ($gff3) ? "3" : "");
          $wormbase->run_command("cat $source_gff >> $target_gff", $log);
          push @to_delete, $source_gff;
        }
      } 
    } else {
      my $target_gff = sprintf("%s/%s.gff%s", $dump_dir, $species, ($gff3) ? "3" : "");
      foreach my $seq (@individual_chrs) {
        my $source_gff = sprintf("%s/%s.gff%s", $dump_dir, $seq, ($gff3) ? "3" : "");
        $wormbase->run_command("cat $source_gff >> $target_gff", $log);
        push @to_delete, $source_gff;
      }
    }
  }

  if (not $log->report_errors) {
    $log->write_to("Appended @to_delete so deleting\n");
    unlink @to_delete;
  } else {
    $log->write_to("Something went wrong with appending @to_delete, so will not delete\n");
  }
}

&stop_giface_server() if @batch_chrs;

$log->mail;
exit(0);


sub start_giface_server {
  $wormbase->run_command("$giface_server $database $port 1200:6000000:1000:600000000 >/dev/null 2>&1 &",$log);
  sleep 20;
}

sub stop_giface_server {
  open (my $write_fh,"| $giface_client $host -port $port -userid wormpub -pass blablub");
  print $write_fh "shutdown now\n";
  close $write_fh;

  # wait a couple of minutes and if the process is still hanging around, kill it
  sleep 180;
  my $ps_string=`ps waux|grep sgiface|grep \$USER|grep -v grep`;
  $ps_string=~/\w+\s+(\d+)/;
  my $server_pid=$1;
  $wormbase->run_command("kill $server_pid",$log) if $server_pid;
}


=pod

=head1 dump_gff_batch.pl

  Use this in conjunction with GFF_method_dump.pl to dump GFF files in parallel using a cluster eg (cbi1)

=head2 SYNOPSIS

  This script is used to create distributed batch jobs running GFF_method_dump.pl.  It builds up a command line including options for said script and submits them to the queueing system

=head2 ARGUMENTS

=over4

  -database:s    - which database to dump data from
  -dump_dir:s    - where to put the output gff files
  -method:s      - comma separated list of methods to dump (does all if not specified)
  -chromosomes:s - comma separated list of chromsosomes to dump (does all if not specified)

=back

=head1 EXAMPLES

=over4

=item perl dump_gff_batch.pl -database ~wormpub/DATABASES/current_DB -dump_dir ~wormpub/GFFdump -method curated,RNAi -chromosome I,II

  will create 4 jobs to dump the following files in ~wormpub/GFFdump
  
=over8

  CHROMOSOME_I.curated.gff
  CHROMOSOME_I.RNAi.gff
  CHROMOSOME_II.curated.gff
  CHROMOSOME_II.RNAi.gff

=back

=over4

=item perl dump_gff_batch.pl -database ~wormpub/DATABASES/current_DB -dump_dir ~wormpub/GFFdump 

  will create 6 jobs to dump everything foreach chromosome.

=back

=cut
