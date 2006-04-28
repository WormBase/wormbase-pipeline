#!/usr/local/bin/perl -w

use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use strict;
use Log_files;
use Storable;

my ($debug, $test, $database);

my @chroms = qw( I II III IV V X MtDNA);
my $dump_dir;
my $dumpGFFscript = $ENV{'CVS_DIR'}."/GFF_method_dump.pl";
my $scratch_dir = "/tmp";
my $methods;
my $chrom_choice;
my $store;


GetOptions (
	    "debug:s"       => \$debug,
	    "test"          => \$test,
	    "database:s"    => \$database,
	    "dump_dir:s"    => \$dump_dir,
	    "method:s"      => \$methods,
	    "methods:s"     => \$methods,
	    "chromosomes:s" => \$chrom_choice,
	    "chromosome:s"  => \$chrom_choice,
	    "store:s"       => \$store
	   );
my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}

my $log = Log_files->make_build_log($wormbase);

$wormbase->checkLSF;

my @methods     = split(/,/,join(',',$methods)) if $methods;
my @chromosomes = split(/,/,join(',',$chrom_choice));



$database = $wormbase->autoace    unless $database;
$dump_dir = $wormbase->gff_splits unless $dump_dir;

@chromosomes = @chroms unless @chromosomes;

$log->write_to("Dumping from DATABASE : $database\n\tto $dump_dir\n\n");
	      
$log->write_to("\t chromosomes".@chromosomes."\n");
if( @methods ){
  $log->write_to("\tmethods ".@methods."\n\n");
}
else {
  $log->write_to("\tno method specified\n\n");
}

$log->write_to("bsub commands . . . . \n\n");
foreach my $chrom ( @chromosomes ) {
  if ( @methods ) {
    foreach my $method ( @methods ) {
      my $err = "$scratch_dir/wormpubGFFdump.$chrom.$method.err";
      my $out = "$scratch_dir/wormpubGFFdump.$chrom.$method.out";
      my $bsub = "bsub -e $err -o $out \"perl $dumpGFFscript -store $store -database $database -dump_dir $dump_dir -chromosome $chrom -method $method\"";
      $log->write_to("$bsub\n");
      $wormbase->run_command("$bsub", $log);
    }
  }
  else {
    my $err = "$scratch_dir/wormpubGFFdump.$chrom.err";
    my $out = "$scratch_dir/wormpubGFFdump.$chrom.out";
    my $bsub = "bsub -e $err -o $out \"perl $dumpGFFscript -store $store -database $database -dump_dir $dump_dir -chromosome $chrom";
    $log->write_to("$bsub\n");
    $wormbase->run_command("$bsub", $log);
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
