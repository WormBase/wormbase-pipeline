#!/usr/local/ensembl/bin/perl -w                  
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2006-11-13 14:27:25 $      

use strict;
use Getopt::Long;
use GDBM_File;


my ($help, $debug, $file);
GetOptions ("help"      => \$help,
            "debug=s"   => \$debug,
	    "file=s"    => \$file,
    );

die "\nYou need to give me a fasta file of the ipi_human proteins with the date appended eg ipi_human_03_05 (5th March)\n
Get this from ftp://ftp.ebi.ac.uk/pub/databases/IPI/current/ipi.HUMAN.fasta.gz\n\n" unless (defined $file);

#make sure on farm-login  -> otherwise GDBM_File fails.
my $host = `hostname`;
chomp $host;
unless ( $host =~ /^bc/ ) {
  die "You must run this from an farm-login machine, otherwise GDBM_File causes a segmentation fault\n\n";
}

#extract date from file name
$file =~ /ipi_human(_\d+_\d+)/;

my $datestamp = $1;
my $output_dir = "/lustre/work1/ensembl/wormpipe/dumps/";
print "processing $file\n" if $debug;
my $fasta = glob("~wormpipe/BlastDB/ipi_human$datestamp.pep");

if( $debug ) {
  $output_dir .= "_test";
}
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

# new UniProt style 
# UniProt/Swiss-Prot:O95180-2|UniProt/TrEMBL:Q9NYY7;Q9NYY6

# remove previous versions of the databases
`rm -f $output_dir/acc2db.dbm` if ( -e "$output_dir/acc2db.dbm" );
`rm -f $output_dir/desc.dbm` if ( -e "$output_dir/desc.dbm" );
`rm -f $output_dir/peptide.dbm` if ( -e "$output_dir/peptide.dbm" );
`rm -f $output_dir/databases.dbm` if ( -e "$output_dir/databases.dbm" );

tie my %ACC2DB, 'GDBM_File',"$output_dir/acc2db.dbm", &GDBM_WRCREAT,0777 or die "cannot open acc2db in $output_dir \n";
tie my %DESC, 'GDBM_File',"$output_dir/desc.dbm",&GDBM_WRCREAT, 0777 or die "cannot open desc in $output_dir \n";
tie my %PEPTIDE,'GDBM_File', "$output_dir/peptide.dbm",&GDBM_WRCREAT, 0777 or die "cannot open peptide in $output_dir \n";
tie my %DATABASES,'GDBM_File', "$output_dir/databases.dbm",&GDBM_WRCREAT, 0777 or die "cannot open databases in $output_dir \n";

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
      if( $acc =~ /^(\w+(-\d+)?);/)
	{ $acc = $1; }
 
      # IPI indicate splice variants with -2 style nomenclature so we strip this and only use the primary version ( -1 )
      if( $acc =~ /(\w+)-(\d+)/ ) {
	if( $2 == 1 ){
	  $acc = $1;
	}
	else {
	  next;
	}	      
      }
      $db =~ s/UniProt\///;
      $databases{$db} = $acc;
    }
    # select primary database
    if( $databases{'ENSEMBL'} ) {
      $prim_DB = "ENSEMBL";
      $prim_DB_id = $databases{'ENSEMBL'};
    }
    elsif( $databases{'Swiss-Prot'} ) {
      $prim_DB = "SW";
      $prim_DB_id = $databases{'Swiss-Prot'};
    }
    elsif( $databases{'TrEMBL'} ){
      $prim_DB = "TR";
      $prim_DB_id = $databases{'TrEMBL'};
    }
    else {
      $dont_read_pep = 1;
      next;
    }

    if( defined $ACC2DB{$prim_DB_id}) {
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
      print OUT "\"\n";
    }
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



untie %ACC2DB;
untie %DESC;
untie %PEPTIDE;
untie %DATABASES;

exit(0);

sub check 
  {
    tie my %DATA, 'GDBM_File',"/lustre/work1/ensembl/wormpipe/dumps/acc2db.dbm",&GDBM_WRCREAT,0777 or die "fail";
    foreach (keys %DATA){
      print "$_ $DATA{$_}\n";
    }
    untie %DATA;
  }
