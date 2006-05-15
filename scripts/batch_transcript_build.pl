#!/usr/local/bin/perl5.8.0 -w
#
# batch_transcript_builder.pl
#
# by Anthony Rogers
#
# wrapper script for running transcript_builder.pl
#
# Last edited by: $Author: gw3 $
# Last edited on: $Date: 2006-05-15 13:04:20 $

use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use strict;
use Coords_converter;
use Storable;

my $dump_dir;
my $database;
my $builder_script = $ENV{CVS_DIR}."/transcript_builder.pl";
my $scratch_dir = "/tmp";
my $chrom_choice;
my $gff_dir;
my ($store, $debug, $test, $no_run);

GetOptions (
	    "database:s"    => \$database,
	    "dump_dir:s"    => \$dump_dir,
	    "gff_dir:s"     => \$gff_dir,
	    "chromosomes:s" => \$chrom_choice,
	    "chromosome:s"  => \$chrom_choice,
	    "store:s"       => \$store,
	    "debug:s"       => \$debug,
	    "test"          => \$test,
	    "no_run"        => \$no_run

	   );

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

my $log = Log_files->make_build_log($wormbase);
$wormbase->checkLSF($log);

$database = $wormbase->autoace unless $database;
unless ( $no_run ){
  my @chromosomes = split(/,/,join(',',$chrom_choice));

  $gff_dir  = $wormbase->gff unless $gff_dir;
  $dump_dir = $wormbase->transcripts unless $dump_dir;
  #MtDNA removed from @chromosomes as the MtDNA Transcripts have been hand curated.
  @chromosomes = qw(I II III IV V X) unless @chromosomes;

  # make a Coords_converter to write the coords files. Otherwise all 6 processes try and do it.
  my $coords = Coords_converter->invoke($database,1, $wormbase);

  # this extract paired read info from the database and writes it to EST_pairs file
  my $cmd = "select cdna, pair from cdna in class cDNA_sequence where exists_tag cdna->paired_read, pair in cdna->paired_read";
  my $tace = $wormbase->tace;
  my $pairs = "$database/EST_pairs.txt";

  open (TACE, "echo '$cmd' | $tace $database |") or die "cant open tace to $database using $tace\n";
  open ( PAIRS, ">$pairs") or die "cant open $pairs :\t$!\n";
  while ( <TACE> ) {
    chomp;
    s/\"//g;
    my @data = split;
    print PAIRS "$data[0]\t$data[1]\n";
  }
  close PAIRS;

  # create and submit LSF jobs.
  foreach my $chrom ( @chromosomes ) {
    my $err = "$scratch_dir/transcipt_builder.$chrom.err.$$";
    my $out = "$dump_dir/CHROMOSOME_${chrom}_transcript.ace";
    my $bsub = "bsub -e $err \"perl $builder_script -database $database -chromosome $chrom -store $store\"";
    print "$bsub\n";
    $wormbase->run_command("$bsub");
  }
}

$log->write_to("waiting for LSF jobs to finish\n");
$wormbase->wait_for_LSF;

$log->write_to("all batch jobs done - cating outputs to ".$wormbase->transcripts."/transcripts.ace\n");
$wormbase->run_command("cp " .$wormbase->misc_dynamic."/MtDNA_Transcripts.ace ".$wormbase->transcripts."/chromosome_MtDNA.ace"); #hard coded cp of fixed MtDNA Transcripts.
$wormbase->run_command("cat ".$wormbase->transcripts."/*.ace > ".$wormbase->transcripts."/transcripts.ace");

$log->write_to("loading file to ".$wormbase->autoace."\n");
$wormbase->load_to_database($wormbase->autoace,$wormbase->transcripts."/transcripts.ace",'transcript_builder', $log);

$log->write_to("batch_dumping GFF files\n");
$wormbase->run_script("dump_gff_batch.pl -method Coding_transcript");

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
