#!/usr/bin/perl -w

use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use strict;

my @chroms = qw( I II III IV V X );
my $dump_dir;
my $database;
my $builder_script = glob("~ar2/wormbase/scripts/transcript_builder.pl");
my $scratch_dir = "/tmp";
my $chrom_choice;
my $gff_dir;

&checkLSF;

GetOptions (
	    "database:s"    => \$database,
	    "dump_dir:s"    => \$dump_dir,
	    "gff_dir:s"     => \$gff_dir,
	    "chromosomes:s" => \$chrom_choice,
	    "chromosome:s"  => \$chrom_choice
	   );

my @chromosomes = split(/,/,join(',',$chrom_choice));

$database = glob("~wormpub/TRANSCRIPTS") unless $database;
$gff_dir  = "$database/GFF" unless $gff_dir;
$dump_dir = "$database/ACE" unless $dump_dir;
@chromosomes = qw(I II III IV V X) unless @chromosomes;


my $cmd = "select cdna, pair from cdna in class cDNA_sequence where exists_tag cdna->paired_read, pair in cdna->paired_read";
my $tace = &tace;
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
   

  

foreach my $chrom ( @chromosomes ) {
  my $err = "$scratch_dir/transcipt_builder.$chrom.err.$$";
  my $out = "$dump_dir/CHROMOSOME_${chrom}_transcript.ace";
  my $bsub = "bsub -e $err \"$builder_script -database $database -chromosome $chrom -gff_dir $gff_dir \"";
  print "$bsub\n";
  system("$bsub");
}


sub checkLSF 
  {
    my $lshosts = `lshosts`;
    die "You need to be on cbi1 or other LSF enabled server to run this" 
      unless $lshosts =~ /HOST_NAME/;
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
