#!/usr/local/bin/perl -w
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2006-02-17 11:20:55 $

use strict;
use lib  $ENV{'CVS_DIR'};
use Wormbase;
use Storable;
use Getopt::Long;
use Log_files;

our ($help, $debug, $test, $stage);
my $store;
my @chromosomes;

GetOptions ("help"         => \$help,
            "debug=s"      => \$debug,
	    "test"         => \$test,
	    "store:s"      => \$store,
	    "stage:s"      => \$stage,
	    "chromosomes:s"=> \@chromosomes,
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
$log->log_and_die("stage not specified\n") unless defined $stage;

my $cmd;
if ( $stage eq 'final' ){
  $cmd = "dump_gff_batch.pl -database ".$wormbase->autoace." -dump_dir ".$wormbase->chromosomes;
  $cmd .= " -chromosomes ". join(",",@chromosomes) if @chromosomes;
}
else {
  my $methods;
  my @methods;
 READARRAY: while (<DATA>) {
    chomp;
    my ($type,$method) = split;
    push(@methods,"$method") if ( $type eq $stage) ;
  }
  $methods = join(',',@methods);

  $log->write_to("Dumping methods $methods from ".$wormbase->autoace."\n");

  $cmd = "dump_gff_batch.pl -database ".$wormbase->autoace." -methods $methods -dump_dir ".$wormbase->autoace."/GFF_SPLITS";
  $cmd .= " -chromosomes ". join(",",@chromosomes) if @chromosomes;
}

$wormbase->run_script($cmd,$log);

$log->mail();
exit(0);

__DATA__
init Genomic_canonical
init Link
init Pseudogene
init Transposon
init Transposon_CDS
init curated
init history
init miRNA
init tRNAscan-SE-1.23
init Non_coding_transcript
init SAGE_transcript
init snRNA
init miRNA
init rRNA
init scRNA
init snoRNA
init tRNA
init stRNA
init snRNA
init GenePairs
init Oligo_set
init SAGE_transcript
init RNAi_primary
init RNAi_secondary
init Expr_profile
blat BLAT_EMBL_BEST
blat BLAT_EMBL_OTHER
blat BLAT_EST_BEST
blat BLAT_EST_OTHER
blat BLAT_NEMATODE
blat BLAT_OST_BEST
blat BLAT_OST_OTHER
blat BLAT_TC1_BEST
blat BLAT_TC1_OTHER
blat BLAT_mRNA_BEST
blat BLAT_mRNA_OTHER
blat BLAT_ncRNA_BEST
blat BLAT_ncRNA_OTHER
blat BLAT_WASHU
blat BLAT_NEMBASE
blat SL1
blat SL2
blat polyA_signal_sequence
blat polyA_site
homol waba_coding
homol waba_strong
homol waba_weak
homol tandem
homol RepeatMasker
homol wublastx
homol inverted
__END__
