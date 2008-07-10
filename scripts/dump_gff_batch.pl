#!/usr/local/bin/perl -w

use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use strict;
use Log_files;
use Storable;

my ($debug, $test, $database,$species);

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
	    "database:s"    => \$database,
	    "dump_dir:s"    => \$dump_dir,
	    "methods:s"     => \$methods,
	    "chromosomes:s" => \$chrom_choice,
	    "store:s"       => \$store,
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
my $scratch_dir = $wormbase->logs;

my @chroms = $wormbase->get_chromosome_names(-prefix => 1,-mito => 1);
$wormbase->checkLSF;

my @methods     = split(/,/,join(',',$methods)) if $methods;

if (scalar @chroms > 50){
	my @bins;
	my $i=0;
	while ($i<scalar @chroms){
		push (@{$bins[$i % 64]},$chroms[$i]);
		$i++;
	}
	map {$_=join(',',@$_)} @bins;
	@chroms = @bins;
}

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
	$wormbase->run_command("(/software/worm/bin/acedb/sgifaceserver $database $port 600:6000000:1000:600000000>/dev/null)>&/dev/null &",$log);
	sleep 5;
}

foreach my $chrom ( @chromosomes ) {
  if ( @methods ) {
    foreach my $method ( @methods ) {
      my $err = scalar(@chromosomes) < 50 ? "$scratch_dir/wormpubGFFdump.$chrom.$method.err" : "$scratch_dir/wormpubGFFdump.$submitchunk.$method.err";
      my $out = scalar(@chromosomes) < 50 ? "$scratch_dir/wormpubGFFdump.$chrom.$method.out" : "$scratch_dir/wormpubGFFdump.$submitchunk.$method.out";
      my @bsub_options = (-e => "$err", -o => "$out");

      my $cmd = "$dumpGFFscript -database $database -dump_dir $dump_dir -chromosome $chrom -method $method -species $species";
      $cmd.=" -host $host" if scalar(@chromosomes) > 50;
      $log->write_to("$cmd\n");
      $cmd = $wormbase->build_cmd($cmd);

      $lsf->submit(@bsub_options, $cmd);
      
    }
  }
  else {
    # for large chromosomes, ask for a file size limit of 2 Gb and a memory limit of 3.5 Gb
    # See: http://scratchy.internal.sanger.ac.uk/wiki/index.php/Submitting_large_memory_jobs
    my @bsub_options = scalar(@chromosomes) < 50 ? (-F => "2000000", -M => "3500000", -R => "\"select[mem>3500] rusage[mem=3500]\"") : ();
    my $err = scalar(@chromosomes) < 50 ? "$scratch_dir/wormpubGFFdump.$chrom.err" : "$scratch_dir/wormpubGFFdump.$submitchunk.err";
    my $out = scalar(@chromosomes) < 50 ? "$scratch_dir/wormpubGFFdump.$chrom.out" : "$scratch_dir/wormpubGFFdump.$submitchunk.out";
    push @bsub_options, (-e => "$err", -o => "$out");

    my $cmd = "$dumpGFFscript -database $database -dump_dir $dump_dir -chromosome $chrom -species $species";
    $cmd.=" -host $host" if scalar(@chromosomes) > 50;
    $log->write_to("$cmd\n");
    print "$cmd\n";
    $cmd = $wormbase->build_cmd($cmd);

    $lsf->submit(@bsub_options, $cmd);

  }
  $submitchunk++;
}

$lsf->wait_all_children( history => 1 );
$log->write_to("All GFF dump jobs have completed!\n");
for my $job ( $lsf->jobs ) {
  $log->error("Job $job (" . $job->history->command . ") exited non zero\n") if $job->history->exit_status != 0;
}
$lsf->clear;

if (scalar(@chromosomes) > 50){
	$wormbase->run_command("(/software/worm/bin/acedb/sgifaceserver $database $port 600:6000000:1000:600000000>/dev/null)>&/dev/null &",$log);
	open (WRITEDB,"| saceclient $host -port $port -userid wormpub -pass yslef4");
	print WRITEDB "shutdown now\n";
	close WRITEDB;
	sleep 120;
	my $ps_string=`ps waux|grep sgiface|grep -v grep`;
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