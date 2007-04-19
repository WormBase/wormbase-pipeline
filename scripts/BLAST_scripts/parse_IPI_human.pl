#!/usr/local/ensembl/bin/perl -w                  
#
# Last updated by: $Author: ar2 $     
# Last updated on: $Date: 2007-04-19 10:41:40 $      

use strict;
use Getopt::Long;
use File::Copy;
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


#>IPI:IPI00000006.1|SWISS-PROT:P01112|ENSEMBL:ENSP00000309845|REFSEQ:NP_005334|H-INV:HIT000080764|
#VEGA:OTTHUMP00000041426;OTTHUMP00000064451;OTTHUMP00000064452;OTTHUMP00000162769;OTTHUMP00000166055
#Tax_Id=9606 Gene_Symbol=HRAS GTPase HRas precursor

# new UniProt style 
# UniProt/Swiss-Prot:O95180-2|UniProt/TrEMBL:Q9NYY7;Q9NYY6

# delete any lefovers 
map {unlink "/tmp/$_"} qw(acc2db.dbm desc.dbm peptide.dbm detebase.dm);

# remove previous versions of the databases
tie my %ACC2DB, 'GDBM_File',"/tmp/acc2db.dbm", &GDBM_WRCREAT,0777 or die "cannot open /tmp/acc2db.dbm\n";
tie my %DESC, 'GDBM_File',"/tmp/desc.dbm",&GDBM_WRCREAT, 0777 or die "cannot open /tmp/desc.dbm \n";
tie my %PEPTIDE,'GDBM_File', "/tmp/peptide.dbm",&GDBM_WRCREAT, 0777 or die "cannot open /tmp/peptide.dbm\n";
tie my %DATABASES,'GDBM_File', "/tmp/databases.dbm",&GDBM_WRCREAT, 0777 or die "cannot open /tmp/databases.dbm\n";

my %database_field = ( 	'VEGA' 		=> 'OtterID',
						'ENSEMBL' 	=> 'ENSEMBL_PEP_ID',
						'SWISS-PROT'=> 'SwissProt_AC',
						'TREMBL' 	=> 'TrEMBL_ID',
						'H-INV' 	=> 'H-INV_ID'
					);

my $prim_DB_id;
my $prim_DB;
my $dont_read_pep;
my $prim_DB_field;

while (<DATA>) {

  if( /^>/ ) {
    my @data = split(/\s+/,$'); # '
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
    elsif( $databases{'VEGA'} ){
      $prim_DB = "VG";
      $prim_DB_id = $databases{'VEGA'};
    }
    elsif( $databases{'TrEMBL'} ){
      $prim_DB = "TR";
      $prim_DB_id = $databases{'TrEMBL'};
    }
    elsif( $databases{'H-INV'} ){
      $prim_DB = "HI";
      $prim_DB_id = $databases{'H-INV'};
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
      print OUT "Database $_ $database_field{$_} $databases{$_}\n";
      $DATABASES{$prim_DB_id} .= "$_:$databases{$_} ";
    }
    if( $data[2] =~ /Gene_Symbol=(\w+)/) {
    	print OUT "Gene_name \"$1\"\n";
    }
    splice @data,0,3;
    print OUT "Title \"".join(" ",@data)."\"\n";#"
	$DESC{$prim_DB_id} .= join(" ",@data);
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

&cleanup($output_dir);

exit(0);

sub check {
    tie my %DATA, 'GDBM_File',"/lustre/work1/ensembl/wormpipe/dumps/acc2db.dbm",&GDBM_WRCREAT,0777 or die "fail";
    foreach (keys %DATA){
      print "$_ $DATA{$_}\n";
    }
    untie %DATA;
}


# moves the databases from tmp to $outputdir
sub cleanup {
	my $to_dir = shift;
	my @dbfiles=qw(acc2db.dbm desc.dbm peptide.dbm databases.dbm);
	foreach my $f (@dbfiles) {
	  system("mv /tmp/$f $to_dir/$f");
	}
}


