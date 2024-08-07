#!/usr/local/bin/perl -w
# Last edited by: $Author: klh $
# Last edited on: $Date: 2015-06-01 11:27:08 $

use strict;
use lib  $ENV{'CVS_DIR'};
use Wormbase;
use Storable;
use Getopt::Long;
use Log_files;

our ($help, $debug, $test, $stage, $gff3, $giface, $giface_server, $giface_client, $cmd, $wormbase, $dump_dir, $species);
my $store;
my @chromosomes;

GetOptions ("help"           => \$help,
            "debug=s"        => \$debug,
	    "test"           => \$test,
	    "store:s"        => \$store,
	    "stage:s"        => \$stage,
	    "chromosomes:s"  => \@chromosomes,
            "gff3"           => \$gff3,
            "giface:s"       => \$giface,
            "gifaceserver:s" => \$giface_server,
            "gifaceclient:s" => \$giface_client,
            "dumpdir:s"      => \$dump_dir,
	    "species:s"      => \$species,
	   );


if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			     -organism => $species,
			   );
}

my $log = Log_files->make_build_log($wormbase);
$log->log_and_die("stage not specified\n") unless defined $stage;

$dump_dir = $wormbase->gff_splits if not defined $dump_dir;

my (@methods);

READARRAY: while (<DATA>) {
  chomp;
  my ($type,$method,@speciestodo) = split;
  push(@methods,"$method") if ( $type eq $stage and (grep($wormbase->species eq $_,@speciestodo)) ) ;
}
my $methods = join(',',@methods);

$log->write_to("Dumping methods $methods from ".$wormbase->autoace."\n");

$cmd = "dump_gff_batch.pl -database ".$wormbase->autoace." -methods $methods -dump_dir $dump_dir -jobtime 00:10:00";
$cmd .= " -chromosomes ". join(",",@chromosomes) if @chromosomes;
$cmd .= " -giface $giface" if defined $giface;
$cmd .= " -gifaceserver $giface_server" if defined $giface_server;
$cmd .= " -gifaceclient $giface_client" if defined $giface_client;

$wormbase->run_script($cmd,$log);

$log->mail();
exit(0);

__DATA__
init Genomic_canonical elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
init Link elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
init Pseudogene elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
init Transposon elegans pristionchus japonica brenneri briggsae remanei brugia ovolvulus sratti tmuris
init Transposon_CDS elegans pristionchus japonica brenneri briggsae remanei brugia ovolvulus sratti tmuris
init Transposon_Pseudogene elegans pristionchus japonica brenneri briggsae remanei brugia ovolvulus sratti tmuris
init Transposon_ncRNA elegans pristionchus japonica brenneri briggsae remanei brugia ovolvulus sratti tmuris
init rRNA_Pseudogene elegans pristionchus japonica brenneri briggsae remanei brugia ovolvulus sratti tmuris
init tRNA_Pseudogene elegans pristionchus japonica brenneri briggsae remanei brugia ovolvulus sratti tmuris
init curated  elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
init history elegans briggsae pristionchus japonica brenneri brugia ovolvulus sratti tmuris
init history_pseudogene elegans briggsae pristionchus japonica brenneri brugia ovolvulus sratti tmuris
init history_transcript elegans briggsae pristionchus japonica brenneri brugia ovolvulus sratti tmuris
init Non_coding_transcript elegans briggsae pristionchus japonica brenneri remanei brugia ovolvulus sratti tmuris
init snlRNA elegans briggsae pristionchus japonica brenneri remanei brugia ovolvulus sratti tmuris
init snRNA elegans briggsae pristionchus japonica brenneri remanei brugia ovolvulus sratti tmuris
init rRNA elegans briggsae pristionchus japonica brenneri remanei brugia ovolvulus sratti tmuris
init scRNA elegans briggsae pristionchus japonica brenneri remanei brugia ovolvulus sratti tmuris
init snoRNA elegans briggsae pristionchus japonica brenneri remanei brugia ovolvulus sratti tmuris
init tRNA elegans briggsae pristionchus japonica brenneri remanei brugia ovolvulus sratti tmuris
init stRNA elegans briggsae pristionchus japonica brenneri remanei brugia ovolvulus sratti tmuris
init ncRNA elegans briggsae brugia japonica brenneri remanei pristionchus ovolvulus sratti tmuris
init miRNA elegans briggsae pristionchus japonica brenneri remanei brugia ovolvulus sratti tmuris
init pre_miRNA elegans briggsae pristionchus japonica brenneri remanei brugia ovolvulus sratti tmuris 
init miRNA_primary_transcript elegans briggsae pristionchus japonica brenneri remanei brugia ovolvulus sratti tmuris
init asRNA elegans briggsae pristionchus japonica brenneri remanei brugia ovolvulus sratti tmuris
init lincRNA elegans briggsae pristionchus japonica brenneri remanei brugia ovolvulus sratti tmuris
init piRNA elegans briggsae pristionchus japonica brenneri remanei brugia ovolvulus sratti tmuris
init circRNA elegans briggsae pristionchus japonica brenneri remanei brugia ovolvulus sratti tmuris
init 7kncRNA elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
blat BLAT_EST_BEST elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
blat BLAT_EST_OTHER elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
blat BLAT_NEMATODE elegans pristionchus remanei briggsae japonica brenneri brugia ovolvulus sratti tmuris
blat BLAT_OST_BEST elegans
blat BLAT_OST_OTHER elegans
blat BLAT_RST_BEST elegans
blat BLAT_RST_OTHER elegans
blat BLAT_TC1_BEST elegans
blat BLAT_TC1_OTHER elegans
blat BLAT_mRNA_BEST elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
blat BLAT_mRNA_OTHER elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
blat BLAT_ncRNA_BEST elegans
blat BLAT_ncRNA_OTHER elegans
blat BLAT_WASHU elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
blat BLAT_NEMBASE elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
blat BLAT_Trinity_BEST elegans briggsae brugia ovolvulus japonica tmuris
blat BLAT_Trinity_OTHER elegans briggsae brugia ovolvulus japonica tmuris
blat BLAT_IsoSeq_BEST tmuris brugia
blat BLAT_IsoSeq_OTHER tmuris brugia
blat BLAT_Nanopore_BEST elegans
blat BLAT_Nanopore_OTHER elegans
blat SL1 elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
blat SL2 elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
blat polyA_signal_sequence elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
blat polyA_site elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
blat GenePairs elegans
blat Orfeome elegans
blat RNAi_primary elegans
blat RNAi_secondary elegans
blat Expr_profile elegans
blat Oligo_set_mapping elegans briggsae remanei brenneri japonica pristionchus
blat Oligo_set elegans
homol waba_coding elegans
homol waba_strong elegans
homol waba_weak elegans
homol tandem elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
homol RepeatMasker elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
homol wublastx elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
homol inverted elegans briggsae remanei pristionchus japonica brenneri brugia ovolvulus sratti tmuris
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
