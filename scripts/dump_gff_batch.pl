#!/usr/local/bin/perl -w

use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use strict;

my @chroms = qw( I II III IV V X );
my $dump_dir;
my $database;
my $dumpGFFscript = glob("~wormpub/TEST_BUILD/scripts/GFF_method_dump.pl");
my $scratch_dir = "/tmp";
my $methods;
my $chrom_choice;

&checkLSF;

GetOptions (
	    "database:s"    => \$database,
	    "dump_dir:s"    => \$dump_dir,
	    "method:s"      => \$methods,
	    "methods:s"     => \$methods,
	    "chromosomes:s" => \$chrom_choice,
	    "chromosome:s"  => \$chrom_choice
	   );

my @methods     = split(/,/,join(',',$methods));
my @chromosomes = split(/,/,join(',',$chrom_choice));

$dump_dir = glob("~wormpub/GFFdump") unless $dump_dir;
$database = glob("~wormpub/DATABASES/current_DB") unless $database;
@chromosomes = @chroms unless @chromosomes;

foreach my $chrom ( @chromosomes ) {
  if ( @methods ) {
    foreach my $method ( @methods ) {
      my $err = "$scratch_dir/wormpubGFFdump.$chrom.$method.err";
      my $out = "$scratch_dir/wormpubGFFdump.$chrom.$method.out";
      my $bsub = "bsub -e $err -o $out \"perl5.6.1 $dumpGFFscript -database $database -dump_dir $dump_dir -chromosome $chrom -method $method\"";
      system("$bsub");
    }
  }
  else {
    my $err = "$scratch_dir/wormpubGFFdump.$chrom.err";
    my $out = "$scratch_dir/wormpubGFFdump.$chrom.out";
    my $bsub = "bsub -e $err -o $out \"perl5.6.1 $dumpGFFscript -database $database -dump_dir $dump_dir -chromosome $chrom";
    system("$bsub");
  }
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
