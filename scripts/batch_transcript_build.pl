#!/usr/local/bin/perl5.8.0 -w
#
# batch_transcript_builder.pl
#
# by Anthony Rogers
#
# wrapper script for running transcript_builder.pl
#
# Last edited by: $Author: gw3 $
# Last edited on: $Date: 2011-05-11 09:00:54 $

use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use strict;
use Coords_converter;
use Storable;
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;

my $dump_dir;
my $database;
my $builder_script = "transcript_builder.pl";
my $scratch_dir = "/tmp";
my $gff_dir;
my ($store, $debug, $test, @no_run, $no_load, $species, @outfile_names);

GetOptions (
	    "database:s"    => \$database,
	    "dump_dir:s"    => \$dump_dir,
	    "gff_dir:s"     => \$gff_dir,
	    "store:s"       => \$store,
	    "debug:s"       => \$debug,
	    "test"          => \$test,
	    "no_run:s"      => \@no_run,
            "noload"        => \$no_load,
	    "species:s"	    => \$species

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
$wormbase->checkLSF($log);

$database = $wormbase->autoace unless $database;

if (@no_run) {
  @outfile_names = @no_run;
} else {
  my @chromosomes = @{$wormbase->get_binned_chroms()}; # no MtDNA

  $gff_dir  = $wormbase->gff_splits unless $gff_dir;
  $dump_dir = $wormbase->transcripts unless $dump_dir;

  # make a Coords_converter to write the coords files. Otherwise all 6 processes try and do it.
  my $coords = Coords_converter->invoke($database,1, $wormbase);

  # clean up output files created in any possible previous runs of this script
  $wormbase->run_command("rm -f ".$wormbase->transcripts."/transcripts*.ace", $log);

  # this extract paired read info from the database and writes it to EST_pairs file
  my $cmd = "select cdna, pair from cdna in class cDNA_sequence where exists_tag cdna->paired_read, pair in cdna->paired_read";
  my $tace = $wormbase->tace;
  my $pairs = "$database/EST_pairs.txt";

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
  
  my $job_name = "worm_".$wormbase->species."_transcript";

  # create and submit LSF jobs.
  $log->write_to("bsub commands . . . . \n\n");
  my $lsf = LSF::JobManager->new();
  foreach my $chrom ( @chromosomes ) {
    my $batchname = $chrom;
    if ($chrom =~ /,/) {
      # if submitting a batch of chromsomes, use the name of the first
      # one - it will be sufficiently unique to identify the batch
      ($batchname) = split(/,/, $chrom);
    }
    my $outfname = "transcripts_${batchname}.ace";
    my $err = "$scratch_dir/transcript_builder.$batchname.err.$$";
    my $cmd = "$builder_script -database $database -chromosome $chrom -acefname $outfname ";
    $log->write_to("$cmd\n");
    print "$cmd\n";
    $cmd = $wormbase->build_cmd($cmd);
    my @bsub_options = (-e => "$err", 
                        -F => "2000000", 
                        -M => "3500000", 
                        -R => "\"select[mem>3500] rusage[mem=3500]\"",
                        -J => $job_name);
    $lsf->submit(@bsub_options, $cmd);
    push @outfile_names, sprintf("%s/%s", $wormbase->transcripts, $outfname);
  }  
  $lsf->wait_all_children( history => 1 );
  $log->write_to("All transcript builder jobs have completed!\n");
  for my $job ( $lsf->jobs ) {
    $log->error("Job $job (" . $job->history->command . ") exited non zero\n") if $job->history->exit_status != 0;
  }
  $lsf->clear;   
}



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
  $wormbase->run_script("dump_gff_batch.pl -method Coding_transcript");

  ##################
  # Check the files
  ##################
  #CHROMOSOME_III  Coding_transcript       protein_coding_primary_transcript       3933602 3935885 .       -       .       Transcript "F37A8.4"

  if ($wormbase->assembly_type ne 'contig') { # elegans, briggsae
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
      } elsif ($wormbase->species eq 'briggsae') {
        
        my %sizes = (
                     'chrI'          => 2000000,
                     'chrI_random'   =>    6000,
                     'chrII'         => 2500000,
                     'chrII_random'  =>  350000,
                     'chrIII'        => 2300000,
                     'chrIII_random' =>   20000,
                     'chrIV'         => 2500000,
                     'chrIV_random'  =>    2000,
                     'chrV'          => 3100000,
                     'chrV_random'   =>   30000,
                     'chrX'          => 3400000,
                     'chrUn'         =>  950000,
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
  } else { # remanei, brenneri, japonica, prist, het
    
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
