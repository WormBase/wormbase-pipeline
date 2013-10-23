
use strict;
use Getopt::Long;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::Method;
use Bio::EnsEMBL::Compara::GenomeDB;
use Bio::EnsEMBL::Compara::SpeciesSet;


my (
  @species,
  %species_data,
  $master_dbname, 
  $host,
  $port,
  $user, 
  $pass,
  $compara_code,
  $collection_name,
  $createmlss,
  $recreatedb,
  $reg_conf,
  $create_tree_mlss,
    ); 
    

GetOptions(
  "host=s" => \$host,
  "port=s" => \$port,
  "user=s" => \$user,
  "pass=s" => \$pass,
  "masterdbname=s"  => \$master_dbname,
  "species=s@"      => \@species,
  "collectionname"  => \$collection_name,
  "comparacode=s"   => \$compara_code,
  "recreate"        => \$recreatedb,
  "reg_conf=s"      => \$reg_conf,
  "treemlss"        => \$create_tree_mlss,
    );


$collection_name = "worms" if not defined $collection_name;

#
# 0. Get species data
#
if (@species) {
  @species = map { split(/,/, $_) } @species;
} else {
  if (not @species) { 
    while(<DATA>) {
      chomp;
      push @species, $_;
    }
  }
}

foreach my $species (sort @species) {
  my $db_name = "worm_ensembl_${species}";
  my $dbh = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host => $host,
    -user => $user,
    -port => $port,
    -pass => $pass,
    -dbname => $db_name);

  if (not defined $dbh) {
    die "Could not connect to database $db_name\n";
  }

  # properties we need for each species:
  # - taxon_id
  # - production_name
  # - assembly
  # - genebuild
  my $taxon_id = $dbh->get_MetaContainer->get_taxonomy_id;
  my $prod_name = $dbh->get_MetaContainer->get_production_name;
  my $genebuild = $dbh->get_MetaContainer->get_genebuild;

  my ($cs) = sort { $a->rank <=> $b->rank } @{$dbh->get_CoordSystemAdaptor->fetch_all};

  $species_data{$species} = {
    taxon_id => $taxon_id,
    production_name => $prod_name,
    genebuild => $genebuild,
    assembly => $cs->version,
  };
}

if ($recreatedb) { 
  #
  # Create database shell
  #
  print STDERR "Re-creating database from scratch\n";

  my $compara_connect = "mysql -u $user -p${pass} -h $host -P $port";
  my $cmd = "$compara_connect -e 'DROP DATABASE IF EXISTS $master_dbname; CREATE DATABASE $master_dbname'";
  system($cmd) and die "Could not re-create existing database\n";
  
  $cmd = "cat $compara_code/sql/table.sql | $compara_connect $master_dbname";
  system($cmd) and die "Could not populate new database with compara schema\n"; 
  
  my $compara_dbh = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    -dbname => $master_dbname,
    -user => $user,
    -pass => $pass,
    -port => $port, 
    -host => $host);
  
  #
  # Populate method_link
  #
  my $ml_adp = $compara_dbh->get_MethodAdaptor();
  my @ml_list = (Bio::EnsEMBL::Compara::Method->new(-type => "ENSEMBL_ORTHOLOGUES", -class => "Homology.homology"), 
                 Bio::EnsEMBL::Compara::Method->new(-type => "ENSEMBL_PARALOGUES", -class => "Homology.homology"), 
                 Bio::EnsEMBL::Compara::Method->new(-type => "PROTEIN_TREES", -class => "ProteinTree.protein_tree_node"));
  foreach my $ml (@ml_list) {
    $ml_adp->store($ml);
  }
}


my $compara_dbh = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
  -dbname => $master_dbname,
  -user => $user,
  -pass => $pass,
  -port => $port, 
  -host => $host);


#
# Populate genome_dbs
#
my @genome_dbs;
foreach my $species (sort { $species_data{$a}->{taxon_id} <=> $species_data{$b}->{taxon_id} } keys %species_data) {
  my $sdata = $species_data{$species};
  my $prod_name = $sdata->{production_name};
  my $assembly = $sdata->{assembly};
  my $genebuild = $sdata->{genebuild};
  my $taxid = $sdata->{taxon_id};

  my $gdb;
  eval {
    $gdb = $compara_dbh->get_GenomeDBAdaptor->fetch_by_name_assembly($prod_name);
  };
  if ($@) {
    # not present; create it
    print STDERR "Creating new GenomeDB for $species\n";
    $gdb = Bio::EnsEMBL::Compara::GenomeDB->new(
      -name      => $prod_name,
      -assembly  => $assembly, 
      -taxon_id  => $taxid,
      -genebuild => $genebuild,
        );

    $compara_dbh->get_GenomeDBAdaptor->store($gdb);
  } else {
    print STDERR "Updating existing GenomeDB for $species\n";
    # update the assembly and genebuild data
    $gdb->assembly($assembly);
    $gdb->genebuild($genebuild);
    $compara_dbh->get_GenomeDBAdaptor->update($gdb);
  }
  push @genome_dbs, $gdb;
}


#
# Make the species set
#
print STDERR "Storing Species set for collection\n";
my $ss = Bio::EnsEMBL::Compara::SpeciesSet->new( -genome_dbs => \@genome_dbs );
$ss->add_tag("name", "collection-${collection_name}");
$compara_dbh->get_SpeciesSetAdaptor->store($ss);

#
# Finally, add the protein-tree MLSSs
#
if ($create_tree_mlss) {
  # For orthologs
  system("perl $compara_code/scripts/pipeline/create_mlss.pl --compara worm_compara_master --reg_conf $reg_conf --collection $collection_name --source wormbase --method_link_type ENSEMBL_ORTHOLOGUES --f --pw") 
      and die "Could not create MLSS for orthologs\n";
  
# For between-species paralogues  
  system("perl $compara_code/scripts/pipeline/create_mlss.pl --compara worm_compara_master --reg_conf $reg_conf --collection $collection_name --source wormbase --method_link_type ENSEMBL_PARALOGUES --f --pw") 
      and die "Could not create MLSS for between-species paralogs\n"; 
  
# For same-species paralogues
  system("perl $compara_code/scripts/pipeline/create_mlss.pl --compara worm_compara_master --reg_conf $reg_conf --collection $collection_name --source wormbase --method_link_type ENSEMBL_PARALOGUES --f --sg") 
      and die "Could not create MLSS for within-species paralogs\n"; 
  
# For protein trees
  system("perl $compara_code/scripts/pipeline/create_mlss.pl --compara worm_compara_master --reg_conf $reg_conf --collection $collection_name --source wormbase --method_link_type PROTEIN_TREES --f") 
      and die "Could not create MLSS for protein trees\n";
}

print STDERR "Updated database\n";
exit(0);

##################################################
__DATA__
elegans
briggsae
brenneri
remanei
japonica
cangaria
csp5
csp11
pristionchus
hcontortus
brugia
asuum
mhapla
heterorhabditis
bxylophilus
sratti
tspiralis
loaloa
panagrellus
dimmitis
