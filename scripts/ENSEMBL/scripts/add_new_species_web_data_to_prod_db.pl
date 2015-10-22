#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBConnection;

my (
  $verbose,
  $host, $port, $user, $pass, $database,
  $mhost, $muser, $mport, $mpass, $mdatabase, 
  $clear_existing,
    );
    


GetOptions( 
  "mhost=s"     => \$mhost,
  "mport=s"     => \$mport,
  "muser=s"     => \$muser,
  "mpass=s"     => \$mpass,
  "mdbname=s"   => \$mdatabase,
  "host=s"      => \$host,
  "user=s"      => \$user,
  "pass=s"      => \$pass,
  "port=i"      => \$port,
  "dbname=s"    => \$database,
  "verbose"     => \$verbose,
  "clearexisting" =>  \$clear_existing,
    );


my $core = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host   => $host, 
    -port   => $port, 
    -user   => $user, 
    -pass   => $pass, 
    -dbname => $database,
    );
    
if (not defined $core) {
  die "Could not connect to core database '$database' using supplied connection details\n";
}


# Get the production database

my $prod = Bio::EnsEMBL::DBSQL::DBConnection->new(
    -host   => $mhost, 
    -port   => $mport, 
    -user   => $muser, 
    -pass   => $mpass, 
    -dbname => $mdatabase,
);

if (not defined $prod) {
  die "Could not connect to production database '$mdatabase' using supplied connection details\n";
}


my $species = $core->get_MetaContainer->get_production_name();

my @analyses = @{$core->get_AnalysisAdaptor->fetch_all};

my $species_prod_sql  = "SELECT species_id FROM species WHERE db_name = ?";

my $res_aref = $prod->sql_helper()->execute(-SQL=>$species_prod_sql, 
                                            -PARAMS => [$species]);
my ($species_id) = @{$res_aref->[0]};

if (not $species_id) {
  die "Could not find species_id for $species\n";
}

if ($clear_existing) {
  print STDERR "DELETING existing web_data for $species\n";

  my $delete_existing_sql = "DELETE FROM analysis_web_data WHERE species_id = ?";

  $res_aref = $prod->sql_helper()->execute_update(-SQL => $delete_existing_sql,
                                                  -PARAMS => [$species_id]);
}  


my $analysis_prod_sql = "SELECT analysis_description_id, default_web_data_id "
    . "FROM analysis_description "
    . "WHERE is_current=1 AND logic_name = ?";

my $species_analysis_sql = "SELECT analysis_web_data_id "
    . "FROM analysis_web_data awd, analysis_description ad "
    . "WHERE awd.analysis_description_id = ad.analysis_description_id AND logic_name = ? AND species_id = ?";
my $analysis_species_insert_sql = "INSERT INTO analysis_web_data (analysis_description_id,web_data_id,species_id,db_type,displayable,created_at,modified_at) "
    . "VALUES (?,?,?,?,?,NOW(),NOW())";


foreach my $analysis (@analyses) {
  
  my $logic_name = $analysis->logic_name;
  my $analysis_id = $analysis->dbID;
	
  print STDERR "logic_name, $logic_name\n";
	
  # check the row is not in there yet
	
  $res_aref = $prod->sql_helper()->execute(-SQL    => $species_analysis_sql, 
                                           -PARAMS => [$logic_name, $species_id]);
  if (@$res_aref > 0) {
    # print STDERR "row already in prod db for logic_name, $logic_name, so no need to add it\n";
  }
  else {
    my $analysis_res_aref = $prod->sql_helper()->execute(-SQL=>$analysis_prod_sql, 
                                                         -PARAMS => [$logic_name]);
	    
    if (@$analysis_res_aref == 0) {
      warn "logic_name, $logic_name, can not be found in the production db, make sure it's been added first, before running this script\n";
      #print "$species\t$analysis_id\t$logic_name\n";
    }
    else {
      my ($analysis_description_id, $default_web_data_id) = @{$analysis_res_aref->[0]};
      
      my $now = time();            
      print STDERR "Inserting row for species_id, analysis_description_id, $species_id, $analysis_description_id\n";
      
      $prod->sql_helper()->execute_update(-SQL=>$analysis_species_insert_sql, -PARAMS => [$analysis_description_id, $default_web_data_id, $species_id, "core", 1]);
    }
  }
}
