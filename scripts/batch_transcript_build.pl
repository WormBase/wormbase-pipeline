#!/usr/local/bin/perl5.8.0 -w
#
# batch_transcript_builder.pl
#
# by Anthony Rogers
#
# wrapper script for running transcript_builder.pl
#
# Last edited by: $Author: klh $
# Last edited on: $Date: 2014-06-17 15:27:11 $

use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use strict;
use Coords_converter;
use Storable;
use Modules::WormSlurm;

my $dump_dir;
my $database;
my $builder_script = "transcript_builder.pl";
my $scratch_dir = "/tmp";
my $gff_dir;
my ($store, $debug, $test, $use_previous_run, $no_load, $species, @outfile_names, $mem, $time);

GetOptions (
    "database:s"       => \$database,
    "dump_dir:s"       => \$dump_dir,
    "gff_dir:s"        => \$gff_dir,
    "store:s"          => \$store,
    "debug:s"          => \$debug,
    "test"             => \$test,
    "usepreviousrun"   => \$use_previous_run,
    "noload"           => \$no_load,
    "species:s"	       => \$species,
    "mem:s"            => \$mem,
    "time:s"           => \$time
    );

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species
			     );
}

my $log = Log_files->make_build_log($wormbase);

unless (defined $mem){
  $mem = '3500m';
}
unless ($mem =~ /^\d+[m|g|M|G]$/) {
    $log->log_and_die('Expecting memory parameter in format /^\d+[m|g|M|G]$/' . "\n");
}

unless (defined $time) {
    $time = '3:00:00';
}
unless ($time =~ /^(\d+\-[0-2]\d:[0-5]\d:[0-5]\d|[0-2]?\d:[0-5]\d:[0-5]\d)$/) {
    $log->log_and_die("Expecting time parameter in format D-HH:MM:SS\n");
}

$wormbase->check_slurm($log);

$database = $wormbase->autoace unless $database;

my @chrs = $wormbase->get_chromosome_names;
my $chunk_total = 24;
$chunk_total = scalar(@chrs) if $chunk_total > scalar(@chrs);

$gff_dir  = $wormbase->gff_splits unless $gff_dir;
$dump_dir = $wormbase->transcripts unless $dump_dir;

# make a Coords_converter to write the coords files. Otherwise all 6 processes try and do it.
my $coords = Coords_converter->invoke($database,1, $wormbase);
# this extract paired read info from the database and writes it to EST_pairs file
my $pairs = "$database/COMMON_DATA/EST_pairs.txt";
if (not -e $pairs) {
  my $cmd = "select cdna, pair from cdna in class cDNA_sequence where exists_tag cdna->paired_read, pair in cdna->paired_read";
  my $tace = $wormbase->tace;

  open (TACE, "echo '$cmd' | $tace $database |") or die "cant open tace to $database using $tace\n";
  open ( PAIRS, ">$pairs") or die "cant open $pairs :\t$!\n";
  while ( <TACE> ) {
    chomp;
    if (/^\s*$/) {next;}
    if (/^Format/) {next;}
    if (/^\/\//) {next;}
    if (/^acedb/) {next;}
    s/Sequence://g;
    my @data = split;
    print PAIRS "$data[0]\t$data[1]\n";
  }
  close PAIRS;
}
  
my $job_name = "worm_".$wormbase->species."_transcript";

# create and submit Slurm jobs.
$log->write_to("Slurm commands . . . . \n\n");
my %jobs;
foreach my $chunk_id (1..$chunk_total) {
  my $batchname = "batch_${chunk_id}";
  my $outfname = "transcripts_${batchname}.ace";
  my $foutname = sprintf("%s/%s", $wormbase->transcripts, $outfname);
  my $pfname = "problems_${batchname}.txt";
  my $fpfname = sprintf("%s/%s", $wormbase->transcripts, $pfname);

  push @outfile_names, sprintf("%s/%s", $wormbase->transcripts, $outfname);
  if (not $use_previous_run or not -e $foutname) {
    unlink $foutname if -e $foutname;
    unlink $fpfname if -e $fpfname;

    my $err = "$scratch_dir/transcript_builder.$batchname.err.$$";
    my $cmd = "$builder_script -database $database -chunkid $chunk_id -chunktotal $chunk_total -acefname $outfname -problemfname $pfname";
    
    $log->write_to("$cmd\n");
    print "$cmd\n";
    $cmd = $wormbase->build_cmd($cmd);

    my $job_id = WormSlurm::submit_job_with_name($cmd, 'production', $mem, $time, '/dev/null', $err, $job_name);
    $jobs{$job_id} = $cmd;
  } else {
    $log->write_to("For batch $chunk_id, using pre-generated acefile $foutname\n");
  }
}  
WormSlurm::wait_for_jobs(keys %jobs);

$log->write_to("All transcript builder jobs have completed!\n");
my $critical_error = 0;
for my $job ( keys %jobs ) {
    my $exit_code = WormSlurm::get_exit_code($job);
  $log->error("Job $job (" . $jobs{$job} . ") exited non zero\n") if $exit_code != 0;
  $critical_error++ if $exit_code != 0;
}

$log->log_and_die("There were $critical_error critical errors in the transcript_builder jobs, please check and re-run with -mem 6000 if you suspect memory issues.\n") if $critical_error;


if ($wormbase->species eq 'elegans') {
  my $source = sprintf("%s/MtDNA_Transcripts.ace", $wormbase->misc_dynamic);
  my $target = sprintf("%s/transcripts_MtDNA.ace", $wormbase->transcripts);
  $wormbase->run_command("cp $source $target", $log);
  push @outfile_names, $target;
}

$log->write_to("all batch jobs done - cating outputs to ".$wormbase->transcripts."/transcripts.ace\n");
my $allfile = sprintf("%s/%s_transcripts.ace",  $wormbase->transcripts, $wormbase->species );
$wormbase->run_command("cat @outfile_names > $allfile", $log);

if (not $no_load) {
  $log->write_to("loading file to ".$wormbase->autoace."\n");
  $wormbase->load_to_database($wormbase->autoace,$allfile,'transcript_builder', $log);
  
  $log->write_to("batch_dumping GFF files\n");
  $wormbase->run_script("dump_gff_batch.pl -method Coding_transcript -time $time", $log);
    
  $log->write_to("Updating common data\n");
  $wormbase->run_script("update_Common_data.pl -worm_gene2geneID", $log);

  ##################
  # Check the files
  ##################
  #CHROMOSOME_III  Coding_transcript       protein_coding_primary_transcript       3933602 3935885 .       -       .       Transcript "F37A8.4"
  
  if ($wormbase->assembly_type ne 'contig') { # elegans
    foreach my $sequence ( $wormbase->get_chromosome_names(-prefix => 1, -mito => 1) ) {
      if($wormbase->species eq 'elegans') {
        
        my %sizes = (
                     'CHROMOSOME_I'       => 4500000,
                     'CHROMOSOME_II'      => 4500000,
                     'CHROMOSOME_III'     => 4000000,
                     'CHROMOSOME_IV'      => 5000000,
                     'CHROMOSOME_V'       => 6000000,
                     'CHROMOSOME_X'       => 5000000,
                     'CHROMOSOME_MtDNA'   =>    2000,
		    );
        $wormbase->check_file("$gff_dir/${sequence}_Coding_transcript.gff", $log,
                              minsize => $sizes{$sequence},
                              lines => ['^##', 
                                        "^${sequence}\\s+Coding_transcript\\s+(protein_coding_primary_transcript|intron|exon)\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+\\s+Transcript\\s+\\S+",
                                        "^${sequence}\\s+Link\\s+region\\s+1\\s+\\d+\\s+\\.\\s+\\+\\s+\\.\\s+Sequence\\s+\\\"${sequence}\\\"",
				       ],
                              gff => 1,
			     );   
	
      }
    }
  } else { # briggsae, remanei, brenneri, japonica, prist, het
    
    $wormbase->check_file("$gff_dir/Coding_transcript.gff", $log,
                          lines => ['^##', 
                                    "^\\S+\\s+Coding_transcript\\s+(protein_coding_primary_transcript|intron|exon)\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+\\s+Transcript\\s+\\S+",
                                    "^\\S+\\s+Link\\s+region\\s+1\\s+\\d+\\s+\\.\\s+\\+\\s+\\.\\s+Sequence\\s+\\\"\\S+\\\"",
                                    "^\\S+\\s+Genomic_canonical\\s+region\\s+1\\s+\\d+\\s+\\.\\s+\\+\\s+\\.\\s+Sequence\\s+\\\"\\S+\\\"",
                                    ],
                          gff => 1,
                          );   
  }
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
