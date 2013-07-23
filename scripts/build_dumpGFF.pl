#!/usr/local/bin/perl -w
# Last edited by: $Author: klh $
# Last edited on: $Date: 2013-07-23 10:40:47 $

use strict;
use lib  $ENV{'CVS_DIR'};
use Wormbase;
use Storable;
use Getopt::Long;
use Log_files;

our ($help, $debug, $test, $stage, $gff3, $giface, $giface_server, $giface_client, $cmd, $wormbase);
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
	   );


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
$cmd .= " -giface $giface" if defined $giface;
$cmd .= " -gifaceserver $giface_server" if defined $giface_server;
$cmd .= " -gifaceclient $giface_client" if defined $giface_client;

$wormbase->run_script($cmd,$log);

$log->mail();
exit(0);

__DATA__
init Genomic_canonical elegans briggsae remanei pristionchus japonica brenneri brugia
init Link elegans briggsae remanei pristionchus japonica brenneri brugia
init Pseudogene elegans briggsae remanei pristionchus japonica brenneri brugia
init Transposon elegans pristionchus japonica brenneri briggsae remanei brugia
init Transposon_CDS elegans pristionchus japonica brenneri briggsae remanei brugia
init Transposon_Pseudogene elegans pristionchus japonica brenneri briggsae remanei brugia
init curated  elegans briggsae remanei pristionchus japonica brenneri brugia
init history elegans briggsae pristionchus japonica brenneri brugia
init history_pseudogene elegans briggsae pristionchus japonica brenneri brugia
init history_transcript elegans briggsae pristionchus japonica brenneri brugia
init tRNAscan-SE-1.23 elegans
init tRNAscan-SE-1.3 briggsae pristionchus japonica brenneri remanei brugia
init Non_coding_transcript elegans briggsae pristionchus japonica brenneri remanei brugia
init snlRNA elegans
init snRNA elegans
init rRNA elegans
init scRNA elegans
init snoRNA elegans
init tRNA elegans
init stRNA elegans
init ncRNA elegans briggsae brugia
init miRNA elegans
init miRNA_primary_transcript elegans
blat BLAT_EST_BEST elegans briggsae remanei pristionchus japonica brenneri brugia
blat BLAT_EST_OTHER elegans briggsae remanei pristionchus japonica brenneri brugia
blat BLAT_NEMATODE elegans pristionchus remanei briggsae japonica brenneri brugia
blat BLAT_OST_BEST elegans
blat BLAT_OST_OTHER elegans
blat BLAT_RST_BEST elegans
blat BLAT_RST_OTHER elegans
blat BLAT_TC1_BEST elegans
blat BLAT_TC1_OTHER elegans
blat BLAT_mRNA_BEST elegans briggsae remanei pristionchus japonica brenneri brugia
blat BLAT_mRNA_OTHER elegans briggsae remanei pristionchus japonica brenneri brugia
blat BLAT_ncRNA_BEST elegans
blat BLAT_ncRNA_OTHER elegans
blat BLAT_WASHU elegans briggsae remanei pristionchus japonica brenneri brugia
blat BLAT_NEMBASE elegans briggsae remanei pristionchus japonica brenneri brugia
blat SL1 elegans briggsae remanei pristionchus japonica brenneri brugia
blat SL2 elegans briggsae remanei pristionchus japonica brenneri brugia
blat polyA_signal_sequence elegans briggsae remanei pristionchus japonica brenneri brugia
blat polyA_site elegans briggsae remanei pristionchus japonica brenneri brugia
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
homol tandem elegans briggsae remanei pristionchus japonica brenneri brugia
homol RepeatMasker elegans briggsae remanei pristionchus japonica brenneri brugia
homol wublastx elegans briggsae remanei pristionchus japonica brenneri brugia
homol inverted elegans briggsae remanei pristionchus japonica brenneri brugia
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
