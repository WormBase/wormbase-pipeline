#!/usr/local/bin/perl5.6.1 -w                  
#
# Last updated by: $Author: wormpipe $     
# Last updated on: $Date: 2003-03-07 15:27:12 $      

use strict;
use Getopt::Long;
use DB_File;


my ($help, $debug, $file);
GetOptions ("help"      => \$help,
	    "file:s"    => \$file,
            "debug=s"   => \$debug);

unless (defined $file) {
  $file = glob("~wormpipe/ipi_human/ipi.HUMAN.fasta");
}
my $output_dir = "/acari/work2a/wormpipe/dumps/";
my $fasta = glob("~wormpipe/BlastDB/ipi_human.pep");
open (DATA, "<$file") or die "cant open $file\n";
open (OUT, ">$output_dir/ipi_human.ace") or die "cant write $output_dir\n";
open (FASTA, ">$fasta") or die "cant write $fasta\n";



#>SWISS-PROT:O43933|TREMBL:Q96S69;Q96S70|ENSEMBL:ENSP00000248633 Tax_Id=9606 Peroxisome biogenes
#is factor 1
#MWGSDRLAGAGGGGAAVTVAFTNARDCFLHLPRRLVAQLHLLQNQAIEVVWSHQPAFLSW
#VEGRHFSDQGENVAEINRQVGQKLGLSNGGQVFLKPCSHVVSCQQVEVEPLSADDWEILE
#LHAVSLEQHLLDQIRIVFPKAIFPVWVDQQTYIFIQIVALIPAASYGRLETDTKLLIQPK
#TRRAKENTFSKADAEYKKLHSYGRDQKGMMKELQTKQLQSNTVGITESNENESEIPVDSS
#SVASLWTMIGSIFSFQSEKKQETSWGLTEINAFKNMQSKVVPLDNIFRVCKSQPPSIYNA

# remove previous versions of the databases
`rm -f $output_dir/acc2db.dbm` if ( -e "$output_dir/acc2db.dbm" );
`rm -f $output_dir/desc.dbm` if ( -e "$output_dir/desc.dbm" );
`rm -f $output_dir/peptide.dbm` if ( -e "$output_dir/peptide.dbm" );
`rm -f $output_dir/databases.dbm` if ( -e "$output_dir/databases.dbm" );

dbmopen my %ACC2DB, "$output_dir/acc2db.dbm", 0777 or die "cannot open acc2db in $output_dir \n";
dbmopen my %DESC, "$output_dir/desc.dbm", 0777 or die "cannot open desc in $output_dir \n";
dbmopen my %PEPTIDE, "$output_dir/peptide.dbm", 0777 or die "cannot open peptide in $output_dir \n";
dbmopen my %DATABASES, "$output_dir/databases.dbm", 0777 or die "cannot open databases in $output_dir \n";

my $prim_DB_id;
my $prim_DB;
my $dont_read_pep;

while (<DATA>) {
  if( /^>/ ) {
    my @data = split(/\s+/,$'); # everything after the >
    my @databases = split(/\|/,$data[0]);
    my %databases;
    
    # if the entry doesnt have a SWALL or ENSEMBL acc dont include it.
    # use $dont_read_pep to prevent sequence being added to previous protein
    $dont_read_pep = 0;

    foreach ( @databases ) {
      my ($db,$acc) = split(/:/,$_);
      next if( ("$db" eq "IPI" ) or ($db =~ /REFSEQ/) );
      if( $acc =~ /(\w+);\S+/ )
	{ $acc = $1;	}
      $databases{$db} = $acc;
    }
    # select primary database
    if( $databases{'ENSEMBL'} ) {
      $prim_DB = "ENSEMBL";
      $prim_DB_id = $databases{'ENSEMBL'};
    }
    elsif( $databases{'SWISS-PROT'} ) {
      $prim_DB = "SW";
      $prim_DB_id = $databases{'SWISS-PROT'};
    }
    elsif( $databases{'TREMBL'} ){
      $prim_DB = "TR";
      $prim_DB_id = $databases{'TREMBL'};
    }
    else {
      $dont_read_pep = 1;
      next;
    }
    
    $ACC2DB{$prim_DB_id} = $prim_DB;   # record database for each accession

    print OUT "\nProtein : \"$prim_DB:$prim_DB_id\"\n";
    print OUT "Peptide \"$prim_DB:$prim_DB_id\"\n";
    foreach (keys %databases) {
      print OUT "Database $_ $databases{$_} $databases{$_}\n";
      $DATABASES{$prim_DB_id} .= "$_:$databases{$_} ";
    }
    my $i = 2;
    if( $data[$i] ) {
      print OUT "Title \"";
      while( $data[$i] ){ 
	print OUT "$data[$i] ";
	$DESC{$prim_DB_id} .= "$data[$i] ";
	$i++;
      }
    }
    print OUT "\"\n";
    print OUT "\nPeptide : \"$prim_DB:$prim_DB_id\"\n";
    print FASTA ">$prim_DB_id\n";
      
  }
  else {
    if( $dont_read_pep == 0 ){
      print OUT ;
      print FASTA;
      $PEPTIDE{$prim_DB_id} .= "$_";
    }
  }
}



dbmclose %ACC2DB;
dbmclose %DESC;
dbmclose %PEPTIDE;
dbmclose %DATABASES;

exit(0);

sub check 
  {
    dbmopen my %DATA, "/acari/work2a/wormpipe/dumps/acc2db.dbm",0777 or die "fail";
    foreach (keys %DATA){
      print "$_ $DATA{$_}\n";
    }
  }
