#!/usr/local/bin/perl5.8.0 -w
#
# batch_transcript_builder.pl
#
# by Anthony Rogers
#
# wrapper script for running transcript_builder.pl
#
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2009-03-05 11:31:49 $

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
my $chrom_choice;
my $gff_dir;
my ($store, $debug, $test, $no_run, $species);

GetOptions (
	    "database:s"    => \$database,
	    "dump_dir:s"    => \$dump_dir,
	    "gff_dir:s"     => \$gff_dir,
	    "chromosomes:s" => \$chrom_choice,
	    "chromosome:s"  => \$chrom_choice,
	    "store:s"       => \$store,
	    "debug:s"       => \$debug,
	    "test"          => \$test,
	    "no_run"        => \$no_run,
	    "species:s"		=> \$species

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
unless ( $no_run ){
  my @chromosomes = @{$wormbase->get_binned_chroms()}; # mo MtDNA
  if ($chrom_choice) {
  	my @chrs = split(/,/,join(',',$chrom_choice)) ;
  	$chromosomes[0] = @chrs;
  }  

  $gff_dir  = $wormbase->gff_splits unless $gff_dir;
  $dump_dir = $wormbase->transcripts unless $dump_dir;

  # make a Coords_converter to write the coords files. Otherwise all 6 processes try and do it.
  my $coords = Coords_converter->invoke($database,undef, $wormbase);

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

  # create and submit LSF jobs.
  $log->write_to("bsub commands . . . . \n\n");
  my $lsf = LSF::JobManager->new();
 	 foreach my $chrom ( @chromosomes ) {
 	 	my $name = substr($chrom,0, 15);
 	    my $err = "$scratch_dir/transcipt_builder.$name.err.$$";
 	    my $cmd = "$builder_script -database $database -chromosome $chrom";
	    $log->write_to("$cmd\n");
	    print "$cmd\n";
	    $cmd = $wormbase->build_cmd($cmd);
	  	my @bsub_options = (-e => "$err", -F => "2000000", -M => "3500000", -R => "\"select[mem>3500] rusage[mem=3500]\"");
	    $lsf->submit(@bsub_options, $cmd);
	  }  
  $lsf->wait_all_children( history => 1 );
  $log->write_to("All transcript builder jobs have completed!\n");
  for my $job ( $lsf->jobs ) {
    $log->error("Job $job (" . $job->history->command . ") exited non zero\n") if $job->history->exit_status != 0;
  }
  $lsf->clear;   
}


$log->write_to("all batch jobs done - cating outputs to ".$wormbase->transcripts."/transcripts.ace\n");
$wormbase->run_command("cp " .$wormbase->misc_dynamic."/MtDNA_Transcripts.ace ".$wormbase->transcripts."/chromosome_MtDNA.ace", $log) if ($wormbase->species eq 'elegans'); #hard coded cp of fixed MtDNA Transcripts.

my $allfile = $wormbase->transcripts."/".$wormbase->species."_transcripts.ace";
$wormbase->run_command("cat ".$wormbase->transcripts."/transcripts*.ace > $allfile", $log);

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
		   'CHROMOSOME_I'       => 4800000,
		   'CHROMOSOME_II'      => 5000000,
		   'CHROMOSOME_III'     => 4500000,
		   'CHROMOSOME_IV'      => 5000000,
		   'CHROMOSOME_V'       => 6000000,
		   'CHROMOSOME_X'       => 5000000,
		   'CHROMOSOME_MtDNA'   =>    2000,
		  );
      $wormbase->check_file("$gff_dir/${sequence}_Coding_transcript.gff", $log,
			    minsize => $sizes{$sequence},
			    lines => ['^##', 
				      "^${sequence}\\s+Coding_transcript\\s+(protein_coding_primary_transcript|intron|exon)\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+Transcript\\s+\\S+",
				      "^${sequence}\\s+Link\\s+region\\s+1\\s+\\d+\\s+\\.\\s+\\+\\s+\\.\\s+Sequence\\s+\\\"${sequence}\\\"",
				     ],
			    gff => 1,
			   );   
    } elsif ($wormbase->species eq 'briggsae') {

      my %sizes = (
		   'chr_I'          => 2000000,
		   'chr_I_random'   =>  500000,
		   'chr_II'         => 2500000,
		   'chr_II_random'  =>  350000,
		   'chr_III'        => 2300000,
		   'chr_III_random' =>  100000,
		   'chr_IV'         => 2500000,
		   'chr_IV_random'  =>   80000,
		   'chr_V'          => 3100000,
		   'chr_V_random'   =>  500000,
		   'chr_X'          => 3400000,
		   'chr_Un'         =>  950000,
		  );
      $wormbase->check_file("$gff_dir/${sequence}_Coding_transcript.gff", $log,
			    minsize => $sizes{$sequence},
			    lines => ['^##', 
				      "^${sequence}\\s+Coding_transcript\\s+(protein_coding_primary_transcript|intron|exon)\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+Transcript\\s+\\S+",
				      "^${sequence}\\s+Link\\s+region\\s+1\\s+\\d+\\s+\\.\\s+\\+\\s+\\.\\s+Sequence\\s+\\\"${sequence}\\\"",
				     ],
			    gff => 1,
			   );   
    }
  }
} else { # remanei, brenneri, japonica, prist, het
      
  $wormbase->check_file("$gff_dir/Coding_transcript.gff", $log,
			lines => ['^##', 
				  "^\\S+\\s+Coding_transcript\\s+(protein_coding_primary_transcript|intron|exon)\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\.\\s+Transcript\\s+\\S+",
				  "^\\S+\\s+Link\\s+region\\s+1\\s+\\d+\\s+\\.\\s+\\+\\s+\\.\\s+Sequence\\s+\\\"\\S+\\\"",
				  "^\\S+\\s+Genomic_canonical\\s+region\\s+1\\s+\\d+\\s+\\.\\s+\\+\\s+\\.\\s+Sequence\\s+\\\"\\S+\\\"",
				 ],
			gff => 1,
		       );   
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
