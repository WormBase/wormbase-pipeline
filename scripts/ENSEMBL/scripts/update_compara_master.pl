#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::Method;
use Bio::EnsEMBL::Compara::GenomeDB;
use Bio::EnsEMBL::Compara::SpeciesSet;
use Bio::EnsEMBL::Compara::Utils::MasterDatabase;

use Bio::EnsEMBL::ApiVersion;

my (
  @species,
  %species_data,
  $compara_code,
  $collection_name,
  $createmlss,
  $recreatedb,
  $reg_conf,
  $create_tree_mlss,
  $tax_host,
  $tax_port,
  $tax_user,
  $tax_dbname,
  $master_host,
  $master_port, 
  $master_user,
  $master_pass, 
  $master_dbname,
  $sfile,
  $no_dnafrags,
  $locators,
  $keep_old_species_sets,
    ); 
    

GetOptions(
  "reg_conf=s"         => \$reg_conf,
  "masterdbname=s"     => \$master_dbname,
  "taxdbname=s"        => \$tax_dbname,
  "collectionname=s"   => \$collection_name,
  "comparacode=s"      => \$compara_code,
  "recreate"           => \$recreatedb,
  "treemlss"           => \$create_tree_mlss,
  "species=s@"         => \@species,
  "sfile=s"            => \$sfile,
  "nodnafrags"         => \$no_dnafrags,
  "locators"           => \$locators,
  "keepoldspeciessets" => \$keep_old_species_sets,
    );

die("must specify registry conf file on commandline\n") unless($reg_conf);
die("Must specify -reg_conf, -masterdbname") unless $reg_conf and $master_dbname;

Bio::EnsEMBL::Registry->load_all($reg_conf);


$collection_name = "wormbase" if not defined $collection_name;

#
# 0. Get species data
#

if (defined $sfile) {
  open(my $fh, $sfile) or die "Could not open species file for reading\n";
  while(<$fh>) {
    next if /^\#/;
    /^(\S+)/ and push @species, $1;
  }
}
@species = map { split(/,/, $_) } @species;

my @core_dbs;
if (@species) {
  @core_dbs = map { Bio::EnsEMBL::Registry->get_DBAdaptor($_, 'core') } @species; 
} else {
  @core_dbs = @{Bio::EnsEMBL::Registry->get_all_DBAdaptors(-GROUP => 'core')};
}
die "No core databases found! " unless @core_dbs;

if ($recreatedb) { 
  die("When creating the db from scratch, you must give -comparacode") unless $compara_code ;

  $tax_dbname = "ncbi_taxonomy" if not defined $tax_dbname;

  print STDERR "Re-creating database from scratch\n";

  {
    my $mdbh = Bio::EnsEMBL::Registry->get_DBAdaptor($master_dbname, 'compara');
    my $taxdbh =  Bio::EnsEMBL::Registry->get_DBAdaptor($tax_dbname, 'compara');

    $master_host = $mdbh->dbc->host;
    $master_port = $mdbh->dbc->port;
    $master_user = $mdbh->dbc->user;
    $master_pass = $mdbh->dbc->password;
    $master_dbname = $mdbh->dbc->dbname;

    $tax_host = $taxdbh->dbc->host;
    $tax_port = $taxdbh->dbc->port;
    $tax_user = $taxdbh->dbc->user;
    $tax_dbname = $taxdbh->dbc->dbname;
  }

  my $compara_connect = "-u $master_user -p${master_pass} -h $master_host -P $master_port";
  my $cmd = "mysql $compara_connect -e 'DROP DATABASE IF EXISTS $master_dbname; CREATE DATABASE $master_dbname'";
  system($cmd) and die "Could not drop and re-create existing database\n";
  
  $cmd = "cat $compara_code/sql/table.sql | mysql $compara_connect $master_dbname";
  system($cmd) and die "Could not populate new database with compara schema\n"; 
  
  $cmd = "mysqlimport --local $compara_connect $master_dbname  $compara_code/sql/method_link.txt";
  system($cmd) and die "Could not populate new database with method_link entries\n"; 
  
  #
  # Populate ncbi taxonomy tables
  #

  open(my $tgt_fh, "| mysql -u $master_user -p$master_pass -h $master_host -P $master_port -D $master_dbname")
      or die "Could not open pipe to target db $master_dbname\n";
  foreach my $table ('ncbi_taxa_node', 'ncbi_taxa_name') {
    open(my $src_fh, "mysqldump -u $tax_user -h $tax_host -P $tax_port ncbi_taxonomy $table |")
        or die "Could not open mysqldump stream from ncbi_taxonomy\n";
    while(<$src_fh>) {
      print $tgt_fh $_;
    }
  }
  close($tgt_fh) or die "Could not successfully close pipe to target db $master_dbname\n";
}


my $compara_dbh = Bio::EnsEMBL::Registry->get_DBAdaptor($master_dbname, 'compara');


#
# Populate genome_dbs
#
my @genome_dbs;
foreach my $core_db (@core_dbs) {

  my $prod_name = $core_db->get_MetaContainer->get_production_name();

  my $gdb;
  eval {
    $gdb = $compara_dbh->get_GenomeDBAdaptor->fetch_by_name_assembly($prod_name);
  };
  if ($@ or not defined $gdb) {
    # not present; create it
    print STDERR sprintf("[%d/%d] Creating new GenomeDB for %s \n", scalar(@genome_dbs)+1, scalar(@core_dbs), $prod_name);
    $gdb = Bio::EnsEMBL::Compara::GenomeDB->new_from_DBAdaptor($core_db);
    $gdb->last_release(undef);
    $gdb->first_release(software_version());
    if ($locators) {
      my $loc = sprintf("Bio::EnsEMBL::DBSQL::DBAdaptor/host=%s;port=%s;user=%s;pass=%s;dbname=%s;disconnect_when_inactive=1", 
                        $core_db->dbc->host, 
                        $core_db->dbc->port, 
                        $core_db->dbc->user, 
                        $core_db->dbc->password, 
                        $core_db->dbc->dbname);
      $gdb->locator($loc);
    }
    $compara_dbh->get_GenomeDBAdaptor->store($gdb);
    Bio::EnsEMBL::Compara::Utils::MasterDatabase::update_dnafrags($compara_dbh, $gdb, $core_db) unless $no_dnafrags;
    #&update_dnafrags($compara_dbh, $gdb, $core_db) unless $no_dnafrags;
  } else {
    print STDERR "Updating existing GenomeDB for $prod_name\n";
    # update the assembly and genebuild data
    my $gdb_tmp =  Bio::EnsEMBL::Compara::GenomeDB->new_from_DBAdaptor( $core_db);
    my $genebuild = $gdb_tmp->genebuild;
    my $assembly = $gdb_tmp->assembly;
    my $old_assembly = $gdb->assembly;
    
    $gdb->genebuild($genebuild);
    $gdb->assembly($assembly);
    $gdb->last_release(undef);

    if ($assembly ne $old_assembly) {
      $gdb->first_release(software_version());
    }
    if ($locators) {
      my $loc = sprintf("Bio::EnsEMBL::DBSQL::DBAdaptor/host=%s;port=%s;user=%s;pass=%s;dbname=%s;disconnect_when_inactive=1", 
                        $core_db->dbc->host, 
                        $core_db->dbc->port, 
                        $core_db->dbc->user, 
                        $core_db->dbc->password, 
                        $core_db->dbc->dbname);
      $gdb->locator($loc);      
    }

    $compara_dbh->get_GenomeDBAdaptor->update($gdb);

    if ($assembly ne $old_assembly and not $no_dnafrags) {
      Bio::EnsEMBL::Compara::Utils::MasterDatabase::update_dnafrags($compara_dbh, $gdb, $core_db) unless $no_dnafrags;
      #&update_dnafrags($compara_dbh, $gdb, $core_db) unless $no_dnafrags;
    }
  }

  push @genome_dbs, $gdb;
}


#
# Make the species set
#

if ($keep_old_species_sets) {
  my $collection_ss = $compara_dbh->get_SpeciesSetAdaptor->fetch_collection_by_name($collection_name);

  # Already have a species set for this collection. Need to check that it is up-to-date
  if ($collection_ss) {
    
    my %collection_gdb_ids = (map {$_->dbID => $_} @{$collection_ss->genome_dbs});
    my %new_gdb_ids = (map {$_->dbID => $_} @genome_dbs);
    
    if (grep { not exists $collection_gdb_ids{$_} } keys %new_gdb_ids or
        grep { not exists $new_gdb_ids{$_} } keys %collection_gdb_ids) {
      # mismatch between collection species set and new species set. 
      print STDERR "Collection composition has changed - updating\n";
      $compara_dbh->get_SpeciesSetAdaptor->update_collection($collection_name, $collection_ss, \@genome_dbs);
      exit(0);
    } else {
      # hunky dory, no update necessary
      print STDERR "Collection composition unchanged - not updating\n";
    }
    
  } else {
    $compara_dbh->get_SpeciesSetAdaptor->update_collection($collection_name, undef, \@genome_dbs);
  }
} else {
  &clean_out_mlss_data( $compara_dbh );
  $compara_dbh->get_SpeciesSetAdaptor->update_collection($collection_name, undef, \@genome_dbs);
}
unless($compara_dbh->get_SpeciesSetAdaptor->fetch_collection_by_name($collection_name)){
  die "Could not add collection for $collection_name";
}

#
# Finally, add the protein-tree MLSSs
#
if ($create_tree_mlss) {
  # For orthologs
  system("perl $compara_code/scripts/pipeline/create_mlss.pl --compara $master_dbname --reg_conf $reg_conf --collection $collection_name --source wormbase --method_link_type ENSEMBL_ORTHOLOGUES --f --pw  -use_genomedb_ids --release") 
      and die "Could not create MLSS for orthologs\n";
  
  # For same-species paralogues
  system("perl $compara_code/scripts/pipeline/create_mlss.pl --compara $master_dbname --reg_conf $reg_conf --collection $collection_name --source wormbase --method_link_type ENSEMBL_PARALOGUES --f --sg --use_genomedb_ids --release") 
      and die "Could not create MLSS for within-species paralogs\n"; 
  
# For protein trees
  system("perl $compara_code/scripts/pipeline/create_mlss.pl --compara $master_dbname --reg_conf $reg_conf --collection $collection_name --source wormbase --method_link_type PROTEIN_TREES --f --name protein_trees_${collection_name} --release ") 
      and die "Could not create MLSS for protein trees\n";

}

print STDERR "Updated database\n";
exit(0);



#############################

sub update_dnafrags {
  my ($compara_dba, $genome_db, $species_dba) = @_;
  die "update_dnafrags got deprecated!";
  my $dnafrag_adaptor = $compara_dba->get_adaptor("DnaFrag");
  my $old_dnafrags = $dnafrag_adaptor->fetch_all_by_GenomeDB_region($genome_db);
  my $old_dnafrags_by_id;
  foreach my $old_dnafrag (@$old_dnafrags) {
    $old_dnafrags_by_id->{$old_dnafrag->dbID} = $old_dnafrag;
  }

  my $sql1 = qq{
      SELECT
        cs.name,
        sr.name,
        sr.length
      FROM
        coord_system cs,
        seq_region sr,
        seq_region_attrib sra,
        attrib_type at
      WHERE
        sra.attrib_type_id = at.attrib_type_id
        AND at.code = 'toplevel'
        AND sr.seq_region_id = sra.seq_region_id
        AND sr.coord_system_id = cs.coord_system_id
    };
  my $sth1 = $species_dba->dbc->prepare($sql1);
  $sth1->execute();

  while (my ($coordinate_system_name, $name, $length) = $sth1->fetchrow_array) {

    #Find out if region is_reference or not
    my $slice = $species_dba->get_SliceAdaptor->fetch_by_region($coordinate_system_name,$name);
    my $is_reference = $slice->is_reference;

    my $new_dnafrag = new Bio::EnsEMBL::Compara::DnaFrag(
            -genome_db => $genome_db,
            -coord_system_name => $coordinate_system_name,
            -name => $name,
            -length => $length,
            -is_reference => $is_reference
        );
    my $dnafrag_id = $dnafrag_adaptor->update($new_dnafrag);
    delete $old_dnafrags_by_id->{$dnafrag_id};
  }

  foreach my $deprecated_dnafrag_id (keys %$old_dnafrags_by_id) {
    $compara_dba->dbc->do("DELETE FROM dnafrag WHERE dnafrag_id = ".$deprecated_dnafrag_id) ;
  }
}

#########################################

sub clean_out_mlss_data {
  my ($compara_dba) = @_;

  foreach my $table ('method_link_species_set',
                     'method_link_species_set_tag',
                     'species_set',
                     'species_set_tag',
                     'species_set_header') {

    my $sth = $compara_dba->dbc->prepare("TRUNCATE TABLE $table");
    $sth->execute();
  }
}
