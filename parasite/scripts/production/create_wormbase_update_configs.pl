
use strict;
use YAML;


foreach my $var (qw(WORMBASEVERSION PARASITEVERSION ENSEMBLVERSION ENSEMBL_CVS_ROOT_DIR PARASITE_CONF)) {
  if (not exists $ENV{$var}) {
    die "$var environment variable not set. You must set up your environment before running this script\n";
  }
}

my $WORMBASEVERSION = $ENV{WORMBASEVERSION};
my $PARASITEVERSION = $ENV{PARASITEVERSION};
my $ENSEMBLVERSION  = $ENV{ENSEMBLVERSION};
my $ENSEMBL_CODE    = $ENV{ENSEMBL_CVS_ROOT_DIR};
my $PARASITE_CONF   = $ENV{PARASITE_CONF};

my $templates = &read_templates();

&write_worm_lite_config($templates->{WORM_LITE}, "$PARASITE_CONF/ensembl_lite.wb_update.conf");
&write_registry($templates->{STAGING_REGISTRY}, "$PARASITE_CONF/staging_pipelines.registry.pm" );
&write_registry($templates->{COMPARA_REGISTRY}, "$PARASITE_CONF/compara.registry.pm" );

######################################################
sub read_templates {
  my (%templates, $current);
  while(<DATA>) {
    if (/^BEGIN_(\S+)_TEMPLATE/) {
      $current = $1;
      next;
    }

    if (/^END_\S+\_TEMPLATE/) {
      $current = undef;
      next;
    }

    $templates{$current} .= $_ if defined $current;
  }

  return \%templates;
}


######################################################
sub write_worm_lite_config {
  my ($template, $yfile) = @_;

  my $conf = YAML::Load($template);
  
  my $pan_con = &get_connection_details("mysql-pan-1");
  my $staging_con = ($PARASITEVERSION % 2) 
      ? &get_connection_details("mysql-ps-staging-1-ensrw") 
      : &get_connection_details("mysql-ps-staging-2-ensrw");
  
  #
  # Generics
  #

  my $generics = $conf->{generics};
  delete $conf->{generics};
  
  $generics->{cvsdir} = $ENSEMBL_CODE;
  
  foreach my $k (keys %$pan_con) {
    $generics->{taxonomy_database}->{$k} = $pan_con->{$k};
    $generics->{production_database}->{$k} = $pan_con->{$k};
  }
  $generics->{"meta.genebuild.version"} = "WS".$WORMBASEVERSION;
  
  
  #
  # Per species
  #
  foreach my $species (sort keys %$conf) {
    foreach my $k (keys %$staging_con) {
      $conf->{$species}->{core_database}->{$k} = $staging_con->{$k};
    }
    my $dbname = $conf->{$species}->{core_database}->{dbname};
    $dbname =~ s/WORMBASEVERSION/$WORMBASEVERSION/g;
    $dbname =~ s/PARASITEVERSION/$PARASITEVERSION/g;
    $dbname =~ s/ENSEMBLVERSION/$ENSEMBLVERSION/g;
    $conf->{$species}->{core_database}->{dbname} = $dbname;
    
    $conf->{$species}->{gff3} =~ s/WORMBASEVERSION/WS$WORMBASEVERSION/g; 
  }
  
  $conf->{generics} = $generics;

  open(my $fh, ">$yfile") or die "Could not open $yfile for reading\n";
  
  print $fh YAML::Dump($conf);
  close($fh);
}

######################################################
sub write_registry {
  my ($template, $reg_file) = @_;

  my %vhosts = ( STAGING => ($PARASITEVERSION % 2) ? "mysql-ps-staging-1" : "mysql-ps-staging-2",
                 PROD    => "mysql-ps-prod",
                 PAN     => "mysql-pan-1" );

  my @replacements = ( { before => "ENSEMBLVERSION", 
                         after => $ENSEMBLVERSION } );

  foreach my $vhost_key (keys %vhosts) {
    my $vhost_name = $vhosts{$vhost_key};

    my $con_ro = &get_connection_details($vhost_name);
    my $con_rw = &get_connection_details("${vhost_name}-ensrw");

    my $url_ro = &get_connection_url($vhost_name);
    my $url_rw = &get_connection_url("${vhost_name}-ensrw");

    push @replacements, { before => "${vhost_key}HOST", 
                          after  => $con_ro->{host} };
    push @replacements, { before => "${vhost_key}PORT", 
                          after  => $con_ro->{port} };
    push @replacements, { before => "${vhost_key}USERRO", 
                          after  => $con_ro->{user} };
    push @replacements, { before => "${vhost_key}USERRW", 
                          after  => $con_rw->{user} };
    push @replacements, { before => "${vhost_key}PASSRW", 
                          after  => $con_rw->{password} };
    push @replacements, { before => "${vhost_key}URLRO", 
                          after  => $url_ro };
    push @replacements, { before => "${vhost_key}URLRW", 
                          after  => $url_rw };
  }

  foreach my $repl (@replacements) {
    my ($bef, $aft) = ($repl->{before}, $repl->{after});

    $template =~ s/$bef/$aft/g; 
  }

  open(my $rfh, ">$reg_file") or die "Could not open $reg_file for writing\n";
  print $rfh $template;
  close($rfh);

}


#########################################
sub get_connection_details {
  my ($db_cmd) = @_;

  open(my $fh, "$db_cmd details script |") or die "Could not run $db_cmd command\n";

  my %conn;
  while(<$fh>) {
    /host\s+(\S+)/ and $conn{host}     = $1;
    /port\s+(\d+)/ and $conn{port}     = $1;
    /user\s+(\S+)/ and $conn{user}     = $1;
    /pass\s+(\S+)/ and $conn{password} = $1;
  }

  return \%conn;
}

##########################################
sub get_connection_url {
  my ($db_cmd) = @_;

  my $url;

  open(my $fh, "$db_cmd details url |") or die "Could not run $db_cmd command\n";
  while(<$fh>) {
    /^(mysql:\S+)/ and $url = $1;
  }

  $url =~ s/\/$//; 

  return $url;
}


__DATA__

BEGIN_WORM_LITE_TEMPLATE
generics: 
  cvsdir: 
  taxonomy_database:
    host: 
    port: 
    user: 
    dbname: ncbi_taxonomy
  production_database:
    host: 
    user: 
    port: 
    dbname: ensembl_production_parasite
  meta.genebuild.version: 
caenorhabditis_elegans:
  seleno: WBGene00015553
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WORMBASEVERSION.annotations.gff3.gz
  core_database:
    dbname: caenorhabditis_elegans_core_PARASITEVERSION_ENSEMBLVERSION_WORMBASEVERSION
caenorhabditis_briggsae:
  seleno: WBGene00028139
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/c_briggsae/PRJNA10731/c_briggsae.PRJNA10731.WORMBASEVERSION.annotations.gff3.gz  
  core_database:
    dbname: caenorhabditis_briggsae_core_PARASITEVERSION_ENSEMBLVERSION_WORMBASEVERSION
caenorhabditis_brenneri:
  seleno: WBGene00158831
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/c_brenneri/PRJNA20035/c_brenneri.PRJNA20035.WORMBASEVERSION.annotations.gff3.gz
  core_database:
    dbname: caenorhabditis_brenneri_core_PARASITEVERSION_ENSEMBLVERSION_WORMBASEVERSION
caenorhabditis_remanei:
  seleno: WBGene00068657
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/c_remanei/PRJNA53967/c_remanei.PRJNA53967.WORMBASEVERSION.annotations.gff3.gz
  core_database:
    dbname: caenorhabditis_remanei_core_PARASITEVERSION_ENSEMBLVERSION_WORMBASEVERSION
caenorhabditis_japonica:
  seleno: WBGene00122465
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/c_japonica/PRJNA12591/c_japonica.PRJNA12591.WORMBASEVERSION.annotations.gff3.gz 
  core_database:
    dbname: caenorhabditis_japonica_core_PARASITEVERSION_ENSEMBLVERSION_WORMBASEVERSION
pristionchus_pacificus:  
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/p_pacificus/PRJNA12644/p_pacificus.PRJNA12644.WORMBASEVERSION.annotations.gff3.gz 
  core_database:
    dbname: pristionchus_pacificus_prjna12644_core_PARASITEVERSION_ENSEMBLVERSION_WORMBASEVERSION
brugia_malayi:
  seleno: WBGene00222286
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/b_malayi/PRJNA10729/b_malayi.PRJNA10729.WORMBASEVERSION.annotations.gff3.gz
  core_database:
    dbname: brugia_malayi_prjna10729_core_PARASITEVERSION_ENSEMBLVERSION_WORMBASEVERSION
onchocerca_volvulus:
  seleno: WBGene00241445
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/o_volvulus/PRJEB513/o_volvulus.PRJEB513.WORMBASEVERSION.annotations.gff3.gz 
  core_database:
    dbname: onchocerca_volvulus_prjeb513_core_PARASITEVERSION_ENSEMBLVERSION_WORMBASEVERSION
strongyloides_ratti:
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/s_ratti/PRJEB125/s_ratti.PRJEB125.WORMBASEVERSION.annotations.gff3.gz
  core_database:
    dbname: strongyloides_ratti_prjeb125_core_PARASITEVERSION_ENSEMBLVERSION_WORMBASEVERSION
END_WORM_LITE_TEMPLATE

BEGIN_STAGING_REGISTRY_TEMPLATE
use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->no_version_check(1);
Bio::EnsEMBL::Registry->no_cache_warnings(1);
{
  Bio::EnsEMBL::Registry->load_registry_from_url('STAGINGURLRW/ENSEMBLVERSION');

  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host    => 'PANHOST',
    -port    => 'PANPORT',
    -user    => 'PANUSERRO',
    -dbname  => 'ensembl_production_parasite',
    -species => 'multi',
    -group   => 'production'
  );
}
1;
END_STAGING_REGISTRY_TEMPLATE

BEGIN_COMPARA_REGISTRY_TEMPLATE
use strict;
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

Bio::EnsEMBL::Registry->load_registry_from_url('STAGINGURLRO/ENSEMBLVERSION');
Bio::EnsEMBL::Registry->load_registry_from_url('PRODURLRW/ENSEMBLVERSION');

Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    -host => 'PANHOST'
    -port => 'PANPORT',
    -user => 'PANUSERRW',
    -pass => 'PANPASSRW',
    -species => 'ensembl_compara_master_parasite',
    -dbname  => 'ensembl_compara_master_parasite',
);

Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    -host => 'PANHOST',
    -user => 'PANUSERRO',
    -port => 'PANPORT',             
    -species => 'ncbi_taxonomy',
    -dbname => 'ncbi_taxonomy',
);
END_COMPARA_REGISTRY_TEMPLATE

__END__










