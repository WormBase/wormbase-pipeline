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
  $division
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
  "division"	       => \$division
    );

die("must specify registry conf file on commandline\n") unless($reg_conf);
# die("Must specify -reg_conf, -masterdbname") unless $reg_conf and $master_dbname;

eval { require($reg_conf) };
$master_dbname = $Parasite::Compara::Registry::MASTER_DB_NAME if not defined $master_dbname;
Bio::EnsEMBL::Registry->load_all($reg_conf);

$collection_name = "wormbase" if not defined $collection_name;
$division = "parasite" if not defined $division;

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
foreach my $s (@species) {
}
  @core_dbs = map { Bio::EnsEMBL::Registry->get_DBAdaptor($_, 'core') } @species;
} else {
  @core_dbs = @{Bio::EnsEMBL::Registry->get_all_DBAdaptors(-GROUP => 'core')};
}
die "No core databases found! " unless @core_dbs;

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
    my $gdb_tmp =  Bio::EnsEMBL::Compara::GenomeDB->new_from_DBAdaptor($core_db);
    my $genebuild = $gdb_tmp->genebuild;
    my $assembly = $gdb_tmp->assembly;
    my $old_assembly = $gdb->assembly;

    $compara_dbh->get_GenomeDBAdaptor->store($gdb);
    Bio::EnsEMBL::Compara::Utils::MasterDatabase::update_dnafrags($compara_dbh, $gdb, $core_db)
    # $gdb->genebuild($genebuild);
    # $gdb->assembly($assembly);
    # $gdb->last_release(undef);
    #
    # if ($assembly ne $old_assembly) {
    #   $gdb->first_release(software_version());
    # }
    # if ($locators) {
    #   my $loc = sprintf("Bio::EnsEMBL::DBSQL::DBAdaptor/host=%s;port=%s;user=%s;pass=%s;dbname=%s;disconnect_when_inactive=1",
    #                     $core_db->dbc->host,
    #                     $core_db->dbc->port,
    #                     $core_db->dbc->user,
    #                     $core_db->dbc->password,
    #                     $core_db->dbc->dbname);
    #   $gdb->locator($loc);
    # }
    #
    # $compara_dbh->get_GenomeDBAdaptor->update($gdb);
    #
    # if ($assembly ne $old_assembly and not $no_dnafrags) {
    #   Bio::EnsEMBL::Compara::Utils::MasterDatabase::update_dnafrags($compara_dbh, $gdb, $core_db) unless $no_dnafrags;
    #   #&update_dnafrags($compara_dbh, $gdb, $core_db) unless $no_dnafrags;
    # }
  }

  push @genome_dbs, $gdb;
}

#
# Make the species set
#
# clean_out_mlss_data( $compara_dbh );
$compara_dbh->get_SpeciesSetAdaptor->update_collection($collection_name, \@genome_dbs, 1);
unless($compara_dbh->get_SpeciesSetAdaptor->fetch_collection_by_name($collection_name)){
  die "Could not add collection for $collection_name";
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
