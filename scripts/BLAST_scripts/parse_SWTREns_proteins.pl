#!/usr/local/bin/perl5.6.1 -w                  
#
# Last updated by: $Author: ar2 $     
# Last updated on: $Date: 2003-01-07 14:38:18 $      

use strict;
use Getopt::Long;

my ($help, $debug, $file);
GetOptions ("help"      => \$help,
	    "file:s"    => \$file,
            "debug=s"   => \$debug);

unless (defined $file) {
  $file = glob("~ar2/9606.SPEns.FASTAC");
}
my $output = glob("~ar2/SWRTEns.ace");
my $fasta = glob("~ar2/SWTREns.fasta");
open (DATA, "<$file") or die "cant open $file\n";
open (OUT, ">$output") or die "cant write $output\n";
open (FASTA, ">$fasta") or die "cant write $fasta\n";


#>SWISS-PROT:O43933|TREMBL:Q96S69;Q96S70|ENSEMBL:ENSP00000248633 Tax_Id=9606 Peroxisome biogenes
#is factor 1
#MWGSDRLAGAGGGGAAVTVAFTNARDCFLHLPRRLVAQLHLLQNQAIEVVWSHQPAFLSW
#VEGRHFSDQGENVAEINRQVGQKLGLSNGGQVFLKPCSHVVSCQQVEVEPLSADDWEILE
#LHAVSLEQHLLDQIRIVFPKAIFPVWVDQQTYIFIQIVALIPAASYGRLETDTKLLIQPK
#TRRAKENTFSKADAEYKKLHSYGRDQKGMMKELQTKQLQSNTVGITESNENESEIPVDSS
#SVASLWTMIGSIFSFQSEKKQETSWGLTEINAFKNMQSKVVPLDNIFRVCKSQPPSIYNA


while (<DATA>) {
  if( /^>/ ) {
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
    my $prim_DB_id;
    my $prim_DB;
    if( $databases{'ENSEMBL'} ) {
      $prim_DB = "ENSEMBL";
      $prim_DB_id = $databases{'ENSEMBL'};
    }
    elsif( $databases{'SWISS-PROT'} ) {
      $prim_DB = "SW";
      $prim_DB_id = $databases{'SWISS-PROT'};
    }
    else {
      $prim_DB = "TR";
      $prim_DB_id = $databases{'TREMBL'};
    }
      
    print OUT "\nProtein : \"$prim_DB:$prim_DB_id\"\n";
    print OUT "Peptide \"$prim_DB:$prim_DB_id\"\n";
    foreach (keys %databases) {
      print OUT "Database $_ $databases{$_} $databases{$_}\n";
    }
    print OUT "Title \"";
    my $i = 2;
    while( $data[$i] ){ 
      print OUT "$data[$i] ";
      $i++;
    }
    print OUT "\"\n";
    print OUT "\nPeptide : \"$prim_DB:$prim_DB_id\"\n";
    print FASTA ">$prim_DB_id\n";
      
  }
  else {
    print OUT ;
    print FASTA;
  }
}


