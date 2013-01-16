#!/usr/local/bin/perl -w
# Last edited by: $Author: klh $
# Last edited on: $Date: 2013-01-16 15:07:46 $

use strict;
use lib  $ENV{'CVS_DIR'};
use Wormbase;
use Storable;
use Getopt::Long;
use Log_files;

our ($help, $debug, $test, $stage, $gff3, $giface, $cmd, $wormbase, $dump_dir);
my $store;
my @chromosomes;

GetOptions ("help"         => \$help,
            "debug=s"      => \$debug,
	    "test"         => \$test,
	    "store:s"      => \$store,
	    "stage:s"      => \$stage,
	    "chromosomes:s"=> \@chromosomes,
            "gff3"         => \$gff3,
            "dumpdir:s"    => \$dump_dir,
            "giface:s"     => \$giface,
	   );


if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}

$dump_dir = $wormbase->gff if not defined $dump_dir;

my $log = Log_files->make_build_log($wormbase);
$log->log_and_die("stage not specified\n") unless defined $stage;


if ( $stage eq 'final' ){
  $cmd = "dump_gff_batch.pl -database ".$wormbase->autoace." -dump_dir $dump_dir";
  $cmd .= " -gff3" if $gff3;
  $cmd .= " -giface $giface" if $giface;
  $cmd .= " -chromosomes ". join(",",@chromosomes) if @chromosomes;
}
else {
  my $methods;
  my @methods;
 READARRAY: while (<DATA>) {
    chomp;
    my ($type,$method,@speciestodo) = split;
    push(@methods,"$method") if ( $type eq $stage and (grep($wormbase->species eq $_,@speciestodo)) ) ;
  }
  $methods = join(',',@methods);

  $log->write_to("Dumping methods $methods from ".$wormbase->autoace."\n");

  $cmd = "dump_gff_batch.pl -database ".$wormbase->autoace." -methods $methods -dump_dir ".$wormbase->gff_splits;
  $cmd .= " -chromosomes ". join(",",@chromosomes) if @chromosomes;
}

$wormbase->run_script($cmd,$log);

$log->mail();
exit(0);

__DATA__
init Genomic_canonical elegans briggsae remanei pristionchus japonica brenneri
init Link elegans briggsae remanei pristionchus japonica brenneri
init Pseudogene elegans briggsae remanei pristionchus japonica brenneri
init Transposon elegans pristionchus japonica brenneri briggsae remanei
init Transposon_CDS elegans pristionchus japonica brenneri briggsae remanei
init Transposon_Pseudogene elegans pristionchus japonica brenneri briggsae remanei
init curated  elegans briggsae remanei pristionchus japonica brenneri
init history elegans briggsae pristionchus japonica brenneri
init history_pseudogene elegans briggsae pristionchus japonica brenneri
init history_transcript elegans briggsae pristionchus japonica brenneri
init tRNAscan-SE-1.23 elegans
init tRNAscan-SE-1.3 briggsae pristionchus japonica brenneri remanei
init Non_coding_transcript elegans briggsae pristionchus japonica brenneri remanei
init SAGE_tag elegans
init snlRNA elegans
init snRNA elegans
init rRNA elegans
init scRNA elegans
init snoRNA elegans
init tRNA elegans
init stRNA elegans
init ncRNA elegans briggsae
init miRNA elegans
init miRNA_primary_transcript elegans
init Oligo_set elegans
init SAGE_transcript elegans
blat BLAT_EST_BEST elegans briggsae remanei pristionchus japonica brenneri
blat BLAT_EST_OTHER elegans briggsae remanei pristionchus japonica brenneri
blat BLAT_NEMATODE elegans pristionchus remanei briggsae japonica brenneri
blat BLAT_OST_BEST elegans
blat BLAT_OST_OTHER elegans
blat BLAT_RST_BEST elegans
blat BLAT_RST_OTHER elegans
blat BLAT_TC1_BEST elegans
blat BLAT_TC1_OTHER elegans
blat BLAT_mRNA_BEST elegans briggsae remanei pristionchus japonica brenneri
blat BLAT_mRNA_OTHER elegans briggsae remanei pristionchus japonica brenneri
blat BLAT_ncRNA_BEST elegans
blat BLAT_ncRNA_OTHER elegans
blat BLAT_WASHU elegans briggsae remanei pristionchus japonica brenneri
blat BLAT_NEMBASE elegans briggsae remanei pristionchus japonica brenneri
blat SL1 elegans briggsae remanei pristionchus japonica brenneri
blat SL2 elegans briggsae remanei pristionchus japonica brenneri
blat polyA_signal_sequence elegans briggsae remanei pristionchus japonica brenneri
blat polyA_site elegans briggsae remanei pristionchus japonica brenneri
blat GenePairs elegans
blat Orfeome elegans
blat RNAi_primary elegans
blat RNAi_secondary elegans
blat Expr_profile elegans
homol waba_coding elegans
homol waba_strong elegans
homol waba_weak elegans
homol tandem elegans briggsae remanei pristionchus japonica brenneri
homol RepeatMasker elegans briggsae remanei pristionchus japonica brenneri
homol wublastx elegans briggsae remanei pristionchus japonica brenneri
homol inverted elegans briggsae remanei pristionchus japonica brenneri
variation Allele elegans
variation CGH_allele elegans
variation Deletion_allele elegans
variation Deletion_and_insertion_allele elegans
variation Insertion_allele elegans
variation KO_consortium_allele elegans
variation Million_mutation elegans
variation Mos_insertion elegans
variation NBP_knockout_allele elegans
variation NemaGENETAG_consortium_allele elegans
variation SNP elegans
variation SNP_Swan elegans
variation SNP_Wicks elegans
variation Substitution_allele elegans
variation Transposon_insertion elegans
variation WGS_Andersen elegans
variation WGS_De_Bono elegans
variation WGS_Hawaiian_Waterston elegans
variation WGS_Hobert elegans
variation WGS_Jarriault elegans
variation WGS_McGrath elegans
variation WGS_Pasadena_Quinlan elegans
variation WGS_Stein elegans
variation WGS_Yanai elegans

__END__
