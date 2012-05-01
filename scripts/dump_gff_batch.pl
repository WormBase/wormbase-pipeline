#!/usr/local/bin/perl -w

use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use strict;
use Log_files;
use Storable;

my ($debug, $test, $database,$species, $verbose, $giface, $gff3);

my $dump_dir;
my $dumpGFFscript = "GFF_method_dump.pl";
my $methods;
my $chrom_choice;
my $store;
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;

GetOptions (
  "debug:s"       => \$debug,
  "test"          => \$test,
  "verbose"       => \$verbose,
  "database:s"    => \$database,
  "dump_dir:s"    => \$dump_dir,
  "methods:s"     => \$methods,
  "chromosomes:s" => \$chrom_choice,
  "store:s"       => \$store,
  "giface:s"      => \$giface,
  "gff3"          => \$gff3,
  "species:s"	    => \$species, # for debug purposes
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

my $log = Log_files->make_build_log($wormbase);
my $scratch_dir = $wormbase->build_lsfout;

my @chroms = $wormbase->get_chromosome_names(-prefix => 1,-mito => 1);
$wormbase->checkLSF;

my @methods     = split(/,/,join(',',$methods)) if $methods;

my @chromosomes = $chrom_choice ? split(/,/,join(',',$chrom_choice)):@chroms;

$database = $wormbase->autoace    unless $database;
$dump_dir = $wormbase->gff_splits unless $dump_dir;

$log->write_to("Dumping from DATABASE : $database\n\tto $dump_dir\n\n");
	      
$log->write_to("\t chromosomes ".@chromosomes."\n");
if( @methods ){
  $log->write_to("\tmethods ".@methods."\n\n");
}
else {
  $log->write_to("\tno method specified\n\n");
}

$log->write_to("bsub commands . . . . \n\n");
my $submitchunk=0;
my $lsf = LSF::JobManager->new();

my $host = qx('hostname');chomp $host;
my $port = 23100;
if (scalar(@chromosomes) > 50){
	$wormbase->run_command("(/nfs/users/nfs_w/wormpub/bin/sgifaceserver_phase_patch $database $port 600:6000000:1000:600000000>/dev/null)>&/dev/null &",$log);
	sleep 20;
}

my $cmd_dir = $wormbase->autoace . "/TMP";
mkdir $cmd_dir, 0777;
my $cmd_file_root = "${cmd_dir}/dump_gff_batch_$$"; # root of name of files to hold commands
my $cmd_number = 0;		# count for making name of next command file
my $store_file = $wormbase->build_store; # get the store file to use in all commands

CHROMLOOP: foreach my $chrom ( @chromosomes ) {
  if ( @methods ) {
    foreach my $method ( @methods ) {
      my $out = scalar(@chromosomes) < 50 ? 
          "$scratch_dir/wormpubGFFdump.$chrom.$method.lsfout" 
          : "$scratch_dir/wormpubGFFdump.$submitchunk.$method.lsfout";
      my $job_name = "worm_".$wormbase->species."_gffbatch";

      my @bsub_options = (-o => "$out",
			  -J => $job_name);

      if (scalar(@chromosomes) < 50) {
        push @bsub_options, (-M => "3500000", 
                             -R => "\"select[mem>3500] rusage[mem=3500]\"");
      }

      my $cmd = "$dumpGFFscript -database $database -dump_dir $dump_dir -method $method -species $species";
      $cmd .= " -host $host" if scalar(@chromosomes) > 50;
      $cmd .= " -debug $debug" if $debug;
      $cmd .= " -chromosome $chrom" if scalar(@chromosomes) < 50;
      $cmd .= " -giface $giface"  if defined $giface;
      $cmd .= " -gff3" if $gff3;
      $cmd = $wormbase->build_cmd_line($cmd, $store_file);      
      $log->write_to("Command: $cmd\n") if ($verbose);
      print "Command: $cmd\n" if ($verbose);

      # write out the command file to be executed
      my $cmd_file = "$cmd_file_root." . $cmd_number++;
      open (CMD, ">$cmd_file") || die "Can't open file $cmd_file: $!";
      print CMD "#!/bin/csh\n";
      print CMD "$cmd\n";
      close (CMD);
      chmod 0777, $cmd_file;

      $lsf->submit(@bsub_options, $cmd_file);
    }
    last CHROMLOOP if scalar(@chromosomes) > 50;
  }
  else {
    # for large chromosomes, ask for a memory limit of 3.5 Gb
    my $job_name = "worm_".$wormbase->species."_gffbatch";
    my @bsub_options = scalar(@chromosomes) < 50 ? (-M => "3500000", 
						    -R => "\"select[mem>3500] rusage[mem=3500]\"",
						   ) : ();
    my $out = scalar(@chromosomes) < 50 
        ? "$scratch_dir/wormpubGFFdump.$chrom.lsfout" 
        : "$scratch_dir/wormpubGFFdump.$submitchunk.lsfout";
    push @bsub_options, (-o => "$out",
			 -J => $job_name);

    my $cmd = "$dumpGFFscript -database $database -dump_dir $dump_dir -species $species";
    $cmd .=" -chromosome $chrom" if scalar(@chromosomes) < 50;
    $cmd .=" -host $host" if scalar(@chromosomes) > 50;
    $cmd .=" -debug $debug" if $debug;
    $cmd .= " -giface $giface"  if defined $giface;
    $cmd .= " -gff3" if $gff3;
    $cmd = $wormbase->build_cmd_line($cmd, $store_file);
    $log->write_to("Command: $cmd\n") if ($verbose);
    print "Command: $cmd\n" if ($verbose);

    # write out the command file to be executed
    my $cmd_file = "$cmd_file_root." . $cmd_number++;
    open (CMD, ">$cmd_file") || die "Can't open file $cmd_file: $!";
    print CMD "#!/bin/csh\n";
    print CMD "$cmd\n";
    close (CMD);
    chmod 0777, $cmd_file;

    $lsf->submit(@bsub_options, $cmd_file);
    last CHROMLOOP if scalar(@chromosomes) > 50;
  }
  $submitchunk++;
}

$lsf->wait_all_children( history => 1 );
$log->write_to("All GFF dump jobs have completed!\n");
my @problem_cmds;
for my $job ( $lsf->jobs ) {
  if ($job->history->exit_status == 0) {
    unlink $job->history->command;
  } else {
    $log->write_to("Job $job (" . $job->history->command . ") exited non zero\n");
    push @problem_cmds, $job->history->command;
  }
}
$lsf->clear;

##################################################################
# now try re-runnning any commands that failed
$lsf = LSF::JobManager->new();
my @bsub_options = scalar(@chromosomes) < 50 ? (-M => "3500000", 
						-R => "\"select[mem>3500] rusage[mem=3500]\""
					       ) : ();

my $out = "$scratch_dir/wormpubGFFdump.rerun.lsfout";
my $job_name = "worm_".$wormbase->species."_gffbatch";
push @bsub_options, (-o => "$out",
		     -J => $job_name);
if (scalar @problem_cmds < 120) { # we don't want to re-run too many jobs!
  foreach my $cmd_file (@problem_cmds) {
    $log->write_to("*** Attempting to re-run job: $cmd_file\n");
    $lsf->submit(@bsub_options, $cmd_file);
  }
  $lsf->wait_all_children( history => 1 );
  my $failed = 0;
  for my $job ( $lsf->jobs ) {
    if ($job->history->exit_status == 0) {
      unlink $job->history->command;
    } else {
      $failed++;
      $log->error("\n\nERROR: Job $job (" . $job->history->command . ") exited non zero at the second attempt. See $out\n");
      push @problem_cmds, $job->history->command;
    }
    $log->write_to("\n\nNumber of jobs that failed after the second attempt: $failed\n");
  }
} else {
  $log->error("ERROR: There are ". scalar @problem_cmds ." GFF_method_dump.pl LSF jobs that have failed. See $out\n");
  $log->write_to("Not attempting to re-run all of these jobs.\nPlease investigate what went wrong!\n");
  for my $job ( $lsf->jobs ) {
    unlink $job->history->command;
  }
}
$lsf->clear;
##################################################################


if (scalar(@chromosomes) > 50){
	open (WRITEDB,"| saceclient $host -port $port -userid wormpub -pass blablub");
	print WRITEDB "shutdown now\n";
	close WRITEDB;
	sleep 180;
	my $ps_string=`ps waux|grep sgiface|grep \$USER|grep -v grep`;
	$ps_string=~/\w+\s+(\d+)/;
	my $server_pid=$1;
	$wormbase->run_command("kill $server_pid",$log) if $server_pid;
}

$log->mail;
exit(0);

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
