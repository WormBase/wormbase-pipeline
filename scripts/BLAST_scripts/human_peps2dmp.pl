#!/usr/local/bin/perl5.6.1 -w
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2003-01-13 16:29:06 $

use strict;
use Getopt::Long;
use DB_File;

my ($help, $debug, $file);
GetOptions ("help"      => \$help,
	    "file:s"    => \$file,
            "debug=s"   => \$debug);

my $DB_dir = "/wormsrv2/wormbase/ensembl_dumps/DBMs";

unless (defined $file) {
  $file = "$DB_dir/9606.SPEns.FASTAC";
}
my $fasta = "$DB_dir/human_pep.fasta";

#wormbase/ensembl_dumps/DBMs
open (DATA, "<$file") or die "cant open $file\n";
#open (FASTA, ">$fasta") or die "cant write $fasta\n";

# DBM files - all keyed off the primary database id (ENS  > SW  > TR)

my $ace_info_dbm = "$DB_dir/ace_info.dbm";
my $sequence_dbm = "$DB_dir/sequence.dbm";
my %ACE_INFO;
my %SEQUENCE;
tie  (%ACE_INFO , 'DB_File',"$ace_info_dbm") or die "cant open $ace_info_dbm :$!\n";

$ACE_INFO{'ant'} = "ROGERS";
#tie  %SEQUENCE , "$sequence_dbm", 0666 or die "cant open sequence_dbm$ :$!\n";

my %swiss_id2gene;
my %acc2id;
&getSwissGeneName;

my %ENSpep_gene;
&makeENSgenes;

#>SWISS-PROT:O43933|TREMBL:Q96S69;Q96S70|ENSEMBL:ENSP00000248633 Tax_Id=9606 Peroxisome biogenesis factor 1
#MWGSDRLAGAGGGGAAVTVAFTNARDCFLHLPRRLVAQLHLLQNQAIEVVWSHQPAFLSW
#VEGRHFSDQGENVAEINRQVGQKLGLSNGGQVFLKPCSHVVSCQQVEVEPLSADDWEILE
#LHAVSLEQHLLDQIRIVFPKAIFPVWVDQQTYIFIQIVALIPAASYGRLETDTKLLIQPK
#TRRAKENTFSKADAEYKKLHSYGRDQKGMMKELQTKQLQSNTVGITESNENESEIPVDSS
#SVASLWTMIGSIFSFQSEKKQETSWGLTEINAFKNMQSKVVPLDNIFRVCKSQPPSIYNA

my $prim_DB_id;
while (<DATA>) {
  if( /^>/ ) {
    chomp;
    my @data = split(/\s+/,$'); # everything after the >
    my @databases = split(/\|/,$data[0]);
    my %databases;
    foreach ( @databases ) {
      my ($db,$acc) = split(/:/,$_);
      if( $acc =~ /(\w+);\S+/ )
	{ $acc = $1;	}
      $databases{$db} = $acc;
    }


    # select primary database
    my $prim_DB;

    my %swissprot;
    my %ensembl;
    my $trembl_id;

    if( $databases{'ENSEMBL'} ) {
      $prim_DB = "ENSEMBL";
      $prim_DB_id = $databases{'ENSEMBL'};
      $ensembl{id} = $databases{'ENSEMBL'};
      $ensembl{gene} = $ENSpep_gene{"$ensembl{id}"};
    }
    if( $databases{'SWISS-PROT'} ) {

      $swissprot{acc} = $databases{'SWISS-PROT'};
      $swissprot{id} = $acc2id{ "$swissprot{acc}" };
      $swissprot{gene} = $swiss_id2gene{ "$swissprot{id}" };

      unless ($prim_DB) {
	$prim_DB = "SW";
	$prim_DB_id = $databases{'SWISS-PROT'};
	$prim_DB_id = $acc2id{"$prim_DB_id"} if (defined $acc2id{"$prim_DB_id"}); # use ID if available
      }
    }

    if ($databases{'TREMBL'}){
      $trembl_id = $databases{'TREMBL'};

      unless ($prim_DB) {
	$prim_DB = "TR";
	$prim_DB_id = $trembl_id;
      }
    }

    my $i = $#data;
    my $description = join(" ",@data[2 .. $i]);
#    $ACE_INFO{"$prim_DB_id"} = "
#Protein : \"$prim_DB:$prim_DB_id\"\n
#Peptide \"$prim_DB:$prim_DB_id\"\n";

    my $database_IDS;
    $database_IDS .= "Database ENSEMBL ENSEMBL:".$ensembl{id}." " if $ensembl{id};
    $database_IDS .= "Database ENSEMBL_GENE ".$ensembl{gene}." " if $ensembl{gene};

    $database_IDS .= "Database SWISSPROT SW:".$swissprot{id}." " if $swissprot{id};
    $database_IDS .= "Database SWISSPROT SW:".$swissprot{acc}." " if $swissprot{acc};
    $database_IDS .= "Database SWISSPROT_GENE ".$swissprot{gene}." " if $swissprot{gene};

    $database_IDS .= "Database TREMBL TR".$trembl_id." " if $trembl_id;
    #print $database_IDS;

    #$ACE_INFO{"$prim_DB_id"} = "$database_IDS";
    #print $ACE_INFO{"$prim_DB_id"};

  }
  else {
   $SEQUENCE{"$prim_DB_id"} .= "$_\n";
  }
}


sub getSwissGeneName
  {
    open (GETZ, "getz -f \"ID PrimAccNumber DBxref\" \"[swissprot-NCBI_TaxId#9606:9606]\" | ");
    my ($id, $acn, $gene);
    my %counts;
    while (<GETZ>) {
      #print $_;
      chomp;
      if( /^ID\s+(\S+)/ ) {
	$id = $1;
	$counts{ids}++;
      }
      elsif( /^AC\s+(\S+) /) {
	$acn = $1;
	$acn =~ s/;//g;
	$acc2id{"$acn"} = $id; 
	$counts{acn}++;
      }
      elsif( /DR\s+Genew;\s+\w+:\d+;\s+(\w+)/ ) {
	    # DR   Genew; HGNC:989; BCL10
	$gene = $1;
	$swiss_id2gene{$id} = $gene;
	$counts{genes}++;
      }
    }
    foreach (keys %swiss_id2gene ) {
      print "ERROR \t\t$_\n" unless $swiss_id2gene{$_};
    }
    foreach (keys %counts) {
      print "$_ $counts{$_}\n";
    }
  }

sub makeENSgenes 
  {
    my $file = glob("~ar2/ensgene_pep");
    open (ENS,"$file") or die "cant open ENSEMBL genes file\n";
    while (<ENS>) {
      #ENSG00000173097	ENSP00000311007
      my @data = split;
      $ENSpep_gene{$data[1]} = $data[0];
    }
  }
