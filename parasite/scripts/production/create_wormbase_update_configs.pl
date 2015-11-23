
use strict;
use YAML;


foreach my $var (qw(WORMBASEVERSION PARASITEVERSION ENSEMBLVERSION ENSEMBL_CVS_ROOT_DIR PARASITE_CONF)) {
  if (not exists $ENV{$var}) {
    die "$var environment variable not set. You must set up your environment before running this script\n";
  }
}

my $WORMBASEVERSION      =  $ENV{WORMBASEVERSION};
my $PARASITEVERSION      = $ENV{PARASITEVERSION};
my $ENSEMBLVERSION       = $ENV{ENSEMBLVERSION};
my $ENSEMBL_CVS_ROOT_DIR = $ENV{ENSEMBL_CVS_ROOT_DIR};
my $PARASITE_CONF        = $ENV{PARASITE_CONF};

my $templates = &read_templates();

&write_config($templates->{WORM_LITE}, "$PARASITE_CONF/ensembl_lite.wb_update.conf");
&write_config($templates->{STAGING_REGISTRY}, "$PARASITE_CONF/staging_pipelines.registry.pm" );
&write_config($templates->{COMPARA_REGISTRY}, "$PARASITE_CONF/compara.registry.pm" );


&write_xref_schema("$ENSEMBL_CVS_ROOT_DIR/ensembl/misc-scripts/xref_mapping");

my $ylite_conf = YAML::LoadFile("$PARASITE_CONF/ensembl_lite.wb_update.conf");
foreach my $spe (keys %$ylite_conf) {
  next if $spe eq 'generics';
  my $template = $templates->{XREF_INPUT};
  &write_config($template, 
                "$PARASITE_CONF/xref_mapping.$spe.input", 
                { before => "COREDBNAME", after => $ylite_conf->{$spe}->{core_database}->{dbname} });
}

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
sub write_config {
  my ($template, $file, @replacements) = @_;
  
  my %vhosts = ( STAGING => ($PARASITEVERSION % 2) ? "mysql-ps-staging-1" : "mysql-ps-staging-2",
                 PROD    => "mysql-ps-prod",
                 PAN     => "mysql-pan-1" );
  
  push @replacements,  ( { before => "ENSEMBLVERSION", 
                           after => $ENSEMBLVERSION },
                         { before => "PARASITEVERSION", 
                           after => $PARASITEVERSION },
                         { before => "WORMBASEVERSION", 
                           after => $WORMBASEVERSION },
                         { before => "ENSEMBL_CVS_ROOT_DIR",
                           after => $ENSEMBL_CVS_ROOT_DIR }
  );
  
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

  open(my $rfh, ">$file") or die "Could not open $file for writing\n";
  print $rfh $template;
  close($rfh);

}

######################################################
sub write_xref_schema {
  my ($xref_code_root) = @_;

  my $conf_ini = "$xref_code_root/xref_config.ini";
  my $conf_script = "$xref_code_root/xref_config2sql.pl";
  my $conf_sql = "$xref_code_root/sql/populate_metadata.sql";

  my $wsrel = "WS".$WORMBASEVERSION;
  
  my $new_conf = "${conf_ini}.${wsrel}";
  
  open(my $infh, $conf_ini) or die "Could not open $conf_ini for reading\n";
  open(my $outfh, ">$new_conf") or die "Could not open $new_conf for writing\n";
  
  my $source;
  while(<$infh>) {
    my $line = $_;
    if (/^\[source (\S+)\]/) {
      $source = $1;
    }
    
    if ($source =~ /^wormbase/) {
      $line =~ s/RELEASE/$wsrel/g; 
    }
    print $outfh $line;
  }
  close($outfh) or die "Could not close $new_conf after writing\n";
  
  open(my $sqlout, ">$conf_sql") or die "Could not open $conf_sql for writing\n";
  open(my $cmd, "perl $conf_script $new_conf |") or die "Could not run $conf_script\n";
  while(<$cmd>) {
    print $sqlout $_;
  }
  close($sqlout) or die "Could not close $conf_sql after writing\n";
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
  cvsdir: ENSEMBL_CVS_ROOT_DIR
  taxonomy_database:
    host: PANHOST
    port: PANPORT
    user: PANUSERRO
    dbname: ncbi_taxonomy
  production_database:
    host: PANHOST
    port: PANPORT
    user: PANUSERRO
    dbname: ensembl_production_parasite
  meta.genebuild.version: WSWORMBASEVERSION
caenorhabditis_elegans:
  seleno: WBGene00015553
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WORMBASEVERSION.annotations.gff3.gz
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    pass: STAGINGPASSRW
    dbname: caenorhabditis_elegans_core_ENSEMBLVERSION_WORMBASEVERSION
caenorhabditis_briggsae:
  seleno: WBGene00028139
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/c_briggsae/PRJNA10731/c_briggsae.PRJNA10731.WORMBASEVERSION.annotations.gff3.gz  
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    pass: STAGINGPASSRW
    dbname: caenorhabditis_briggsae_core_ENSEMBLVERSION_WORMBASEVERSION
caenorhabditis_brenneri:
  seleno: WBGene00158831
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/c_brenneri/PRJNA20035/c_brenneri.PRJNA20035.WORMBASEVERSION.annotations.gff3.gz
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    pass: STAGINGPASSRW
    dbname: caenorhabditis_brenneri_core_ENSEMBLVERSION_WORMBASEVERSION
caenorhabditis_remanei:
  seleno: WBGene00068657
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/c_remanei/PRJNA53967/c_remanei.PRJNA53967.WORMBASEVERSION.annotations.gff3.gz
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    pass: STAGINGPASSRW
    dbname: caenorhabditis_remanei_core_ENSEMBLVERSION_WORMBASEVERSION
caenorhabditis_japonica:
  seleno: WBGene00122465
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/c_japonica/PRJNA12591/c_japonica.PRJNA12591.WORMBASEVERSION.annotations.gff3.gz 
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    pass: STAGINGPASSRW
    dbname: caenorhabditis_japonica_core_ENSEMBLVERSION_WORMBASEVERSION
pristionchus_pacificus:  
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/p_pacificus/PRJNA12644/p_pacificus.PRJNA12644.WORMBASEVERSION.annotations.gff3.gz 
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    pass: STAGINGPASSRW
    dbname: pristionchus_pacificus_prjna12644_core_PARASITEVERSION_ENSEMBLVERSION_WORMBASEVERSION
brugia_malayi:
  seleno: WBGene00222286
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/b_malayi/PRJNA10729/b_malayi.PRJNA10729.WORMBASEVERSION.annotations.gff3.gz
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    pass: STAGINGPASSRW
    dbname: brugia_malayi_prjna10729_core_PARASITEVERSION_ENSEMBLVERSION_WORMBASEVERSION
onchocerca_volvulus:
  seleno: WBGene00241445
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/o_volvulus/PRJEB513/o_volvulus.PRJEB513.WORMBASEVERSION.annotations.gff3.gz 
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    pass: STAGINGPASSRW
    dbname: onchocerca_volvulus_prjeb513_core_PARASITEVERSION_ENSEMBLVERSION_WORMBASEVERSION
strongyloides_ratti:
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WORMBASEVERSION/species/s_ratti/PRJEB125/s_ratti.PRJEB125.WORMBASEVERSION.annotations.gff3.gz
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    pass: STAGINGPASSRW
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

BEGIN_XREF_INPUT_TEMPLATE
xref
host=PRODHOST
port=PRODPORT
user=PRODUSERRW
password=PRODPASSRW
dbname=xref_parasite_SPECIES_PARASITEVERSION
dir=/nfs/nobackup/ensemblgenomes/wormbase/parasite/production/xrefs/mapping/SPECIES

species=caenorhabditis_elegans
taxon=wormbase
host=STAGINGHOST
port=STAGINGPORT
user=STAGINGUSERRW
password=STAGINGPASSRW
dbname=COREDBNAME
dir=/nfs/nobackup/ensemblgenomes/wormbase/parasite/production/xrefs/mapping/SPECIES

farm
queue=production-rh6
exonerate=/nfs/panda/ensemblgenomes/external/exonerate-2/bin/exonerate

END_XREF_INPUT_TEMPLATE

__END__










