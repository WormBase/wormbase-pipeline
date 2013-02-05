#!/usr/bin/env perl 
#

use strict;



use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my (
  $dbname,
  $dbhost,
  $dbuser,
  $dbport,
  $dbpass,
  @hit_names,
    );


&GetOptions(
	'dbname=s' => \$dbname,
	'dbuser=s' => \$dbuser,
	'dbhost=s' => \$dbhost,
	'dbport=s' => \$dbport,
	'dbpass=s' => \$dbpass
);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	'-dbname' => $dbname,
	'-host' => $dbhost,
	'-user' => $dbuser,
	'-port' => $dbport,
	'-pass' => $dbpass
);


my $ipi_flatfile = shift;

my $ipi_ensembl_map = &get_ipi_ensembl_map($ipi_flatfile);

my $ana = $db->get_AnalysisAdaptor->fetch_by_logic_name('ipi_humanx');
my $ana_db_id = $ana->dbID;

my $sth = $db->dbc->prepare("SELECT distinct(hit_name) from protein_align_feature " .
                            "WHERE analysis_id = $ana_db_id");


$sth->execute;
while(my ($val) = $sth->fetchrow_array) {
  push @hit_names, $val;
  if (not exists $ipi_ensembl_map->{$val}) {
    die "Could not find IPI id for $val\n";
  }
}
$sth->finish;

$sth = $db->dbc->prepare("UPDATE protein_align_feature set hit_name = ? ".
                         "WHERE analysis_id = $ana_db_id AND hit_name = ?");
foreach my $id (@hit_names) {
  my $change_to = $ipi_ensembl_map->{$id};
  print "Will update $id to $change_to\n";
  $sth->execute($change_to, $id);
}
$sth->finish;



sub get_ipi_ensembl_map {
  my $file = shift;
  
  my ($fh, %acc_map);

  if ($file =~ /\.gz$/) {
    open($fh, "gunzip -c $file |") or die "Could not open gzip stream to $file\n";
  } else {
    open($fh, $file) or die "Could not open $file for reading\n";
  }

  while(<$fh>) {  
    if( /^>/ ) {
      
      my ($dbs) = ($_ =~ /^\>(\S+)/);
      
      my @databases = split(/\|/,$dbs);
      my %databases;
      
      my ($ipi_acc, $prim_DB, $prim_DB_id);
      
      
      foreach ( @databases ) {
        my ($db,$acc) = split(/:/,$_);
        
        if ($db eq 'IPI') {
          $ipi_acc = $acc;
          $ipi_acc =~ s/\.\d+//;
          next;
        }
        
        next if( ($db =~ /REFSEQ/) );
        
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
      
      $acc_map{$prim_DB_id} = $ipi_acc;
    }
    
  }

  return \%acc_map;
}

