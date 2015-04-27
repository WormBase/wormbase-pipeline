#!/usr/bin/env perl

#
# script to add external_db_ids to the hit_names
# in protein- & dna-align_feature tables
#
#


use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my %mapping = (
  
  brepepx        => 'wormbase_id',
  brigpepx       => 'wormbase_id',
  ppapepx        => 'wormbase_id',  
  wormpepx       => 'wormbase_id', 
  jappepx        => 'wormbase_id',  
  remapepx       => 'wormbase_id',  
  brugpepx       => 'wormbase_id',  
  ovolpepx       => 'wormbase_id',  
  ipi_humanx     => 'IPI',
  gadflyx        => 'gadfly_translation_cgid',
  yeastx         => 'SGD',  
  slimswissprotx => 'Uniprot/SWISSPROT',
  slimtremblx    => 'Uniprot/SPTREMBL',

  celegans_est   => 'wormbase_id',
  celegans_mrna  => 'wormbase_id',
  cbriggsae_est   => 'wormbase_id',
  cbriggsae_mrna  => 'wormbase_id',
  cjaponica_est   => 'wormbase_id',
  cjaponica_mrna  => 'wormbase_id',
  cremanei_est   => 'wormbase_id',
  cremanei_mrna  => 'wormbase_id',
  cbrenneri_est   => 'wormbase_id',
  cbrenneri_mrna  => 'wormbase_id',
  ppacificus_est   => 'wormbase_id',
  ppacificus_mrna  => 'wormbase_id',
  ovolvulus_est   => 'wormbase_id',
  ovolvulus_mrna  => 'wormbase_id',
  sratti_est   => 'wormbase_id',
  sratti_mrna  => 'wormbase_id',
  celegans_ost     => 'wormbase_id',
  celegans_rst     => 'wormbase_id',

);


my ($host, $user, $pass, $dbname, $port, $conf_file);

&GetOptions(
  'dbhost|host:s'  => \$host,
  'dbuser|user:s'  => \$user,
  'dbpass|pass:s'  => \$pass,
  'port|dbport:s'  => \$port,
  'dbname:s'       => \$dbname,
  'conf:s'         => \$conf_file,
);

if (defined $conf_file) {
  %mapping = ();

  open(my $conffh, $conf_file) or die "Could not open $conf_file for reading\n";
  while(<$conffh>) {
    /^(\S+)\s+(\S+)/ and do {
      $mapping{$1} = $2;
    };
  }
}


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $host,
  -user   => $user,
  -port   => $port,
  -dbname => $dbname,
  -pass   => $pass,
) or die "can not connect to $dbname@$host\n";


#Get list of external dbs
my $extdb_hash = &get_extdb_names($db);

my $prot_sql = "UPDATE protein_align_feature SET external_db_id = ? WHERE analysis_id = ?";
my $dna_sql = "UPDATE dna_align_feature SET external_db_id = ? WHERE analysis_id = ?";

my $prot_sth = $db->dbc->prepare($prot_sql);
my $dna_sth = $db->dbc->prepare($dna_sql);

foreach my $logic_name (keys %mapping) {
  my $ext_db = $mapping{$logic_name};
  my $ana = $db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);

  if (not defined $ana) {
    warn("Could not find analysis with logic_name '$logic_name' so skipping it");
    next;
  }
  
  if (not exists $extdb_hash->{$ext_db}) {
    warn("Could not find external_db entry with name '$ext_db' - skipping");
    next;
  }

  my $ana_id = $ana->dbID;
  my $ext_db_id = $extdb_hash->{$ext_db};

  $prot_sth->execute($ext_db_id, $ana_id);
  $dna_sth->execute($ext_db_id, $ana_id);
}

$prot_sth->finish;
$dna_sth->finish;

exit(0);


###################################
sub get_extdb_names {
  my ($db) = @_;

  my %namehash;

  my $query = "select db_name, external_db_id from external_db " ;
  my $sth = $db->dbc->prepare($query) ;
  $sth->execute() ;

  while ( my $arrayref = $sth->fetchrow_arrayref() ) {
    my ( $name, $dbid ) = @$arrayref ;
    $namehash{$name} = $dbid ;
  }

  $sth->finish ;

  return \%namehash;
}

sub get_logic_names_hash {
  my ($db) = @_;

  my $aa = $db->get_AnalysisAdaptor;

  my @anals = @{$aa->fetch_all};

  my %analhash;
  foreach my $anal (@anals) {
    $analhash{$anal->logic_name} = $anal;
  }

  return \%analhash;
}


1;
