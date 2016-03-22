
use strict;
use YAML;
use Getopt::Long;
use File::Path qw(mkpath);

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my ($write_configs, 
    $update_cores,
    $xref_parsing,
    $xref_mapping,
    $xref_loading,
    $xref_cleanup,
    $all_stages,
    $verbose,
    @species);

foreach my $var (qw(WORMBASE_VERSION PARASITE_VERSION ENSEMBL_VERSION ENSEMBL_CVS_ROOT_DIR PARASITE_CONF PARASITE_SCRATCH WORM_CODE)) {
  if (not exists $ENV{$var}) {
    die "$var environment variable not set. You must set up your environment before running this script\n";
  }
}

my $WORMBASE_VERSION     = $ENV{WORMBASE_VERSION};
my $PARASITE_VERSION     = $ENV{PARASITE_VERSION};
my $ENSEMBL_VERSION      = $ENV{ENSEMBL_VERSION};
my $ENSEMBL_CVS_ROOT_DIR = $ENV{ENSEMBL_CVS_ROOT_DIR};
my $PARASITE_CONF        = $ENV{PARASITE_CONF};
my $WORM_CODE            = $ENV{WORM_CODE};
my $PARASITE_SCRATCH     = $ENV{PARASITE_SCRATCH};

my $WORK_DIR   = "$PARASITE_SCRATCH/core_update";
my $LOG_DIR    = "$WORK_DIR/logs";
my $BACKUP_DIR = "$WORK_DIR/dbbackups";

&GetOptions("species=s" => \@species,
            "writeconfig" => \$write_configs,
            "updatecores" => \$update_cores,
            "xrefparsing" => \$xref_parsing,
            "xrefmapping" => \$xref_mapping,
            "xrefloading" => \$xref_loading,
            "xrefcleanup" => \$xref_cleanup,
            "allstages"   => \$all_stages,
            "verbose"     => \$verbose,
    );


&write_configs(@species) if $write_configs or $all_stages;
&update_cores(@species) if $update_cores or $all_stages;
&xref_parsing(@species) if $xref_parsing or $all_stages;
&xref_mapping(@species) if $xref_mapping or $all_stages;
&xref_loading(@species) if $xref_loading or $all_stages;
&xref_cleanup(@species) if $xref_cleanup or $all_stages;


#####################################################
sub write_configs {
  my $templates = &read_templates();
  
  mkpath($PARASITE_CONF, { verbose => 0, mode => 0775 }) if not -e $PARASITE_CONF; 

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
                  { before => "COREDBNAME", 
                    after => $ylite_conf->{$spe}->{core_database}->{dbname} },
                  { before => "SPECIES",
                    after => $spe });
  }
}

#####################################################
sub update_cores {
  my @species = @_;

  foreach my $spe (@species) {
    my $cmd = "perl $WORM_CODE/scripts/ENSEMBL/scripts/worm_lite.pl -yfile $PARASITE_CONF/ensembl_lite.wb_update.conf -load_genes -species $spe";
    &write_log("Running: $cmd\n");
    system($cmd) and die "Command failure: $cmd\n";

    $cmd = "perl $WORM_CODE/scripts/ENSEMBL/scripts/worm_lite.pl -yfile $PARASITE_CONF/ensembl_lite.wb_update.conf -load_meta -species $spe";

    &write_log("Running: $cmd\n");
    system($cmd) and die "Command failure: $cmd\n";

    my $xref_attr = &parse_xref_inputconf("$PARASITE_CONF/xref_mapping.$spe.input");

    my $backup = sprintf("%s/%s.post_parsing.sql.gz", 
                         $BACKUP_DIR,
                         $xref_attr->{xref}->{dbname});
    my $backup_cmd = sprintf("mysqldump --host=%s --port=%s --user=%s --password=%s %s | gzip > %s", 
                             $xref_attr->{xref}->{host},
                             $xref_attr->{xref}->{port},
                             $xref_attr->{xref}->{user},
                             $xref_attr->{xref}->{password},
                             $xref_attr->{xref}->{dbname},
                             $backup);

  }
}

#####################################################
sub xref_parsing {
  my @species = @_;

  foreach my $spe (@species) {

    my $xref_attr = &parse_xref_inputconf("$PARASITE_CONF/xref_mapping.$spe.input");

    my $download_dir = "$WORK_DIR/xrefs/sources/$spe";
    system("rm -fr $download_dir") and die "Could not remove $download_dir for re-parsing\n";
    mkpath($download_dir, { verbose => 0, mode => 0775 });

    my $cmd = sprintf("cd %s/ensembl/misc-scripts/xref_mapping && perl xref_parser.pl --host %s --port %s --user %s --pass %s --dbname %s --download_dir %s -species %s -drop_db -delete_downloaded -create -stats > %s 2>&1",
                      $ENSEMBL_CVS_ROOT_DIR,
                      $xref_attr->{xref}->{host},
                      $xref_attr->{xref}->{port},
                      $xref_attr->{xref}->{user},
                      $xref_attr->{xref}->{password},
                      $xref_attr->{xref}->{dbname},
                      $download_dir,
                      $spe,
                      "$LOG_DIR/parsing.$spe.WS${WORMBASE_VERSION}.out");
    
    &write_log("Running: $cmd\n");
    system($cmd) and die "Command failure: $cmd\n";
    
    my $backup = sprintf("%s/%s.pre_xref.sql.gz", 
                         $BACKUP_DIR,
                         $xref_attr->{species}->{dbname});
    my $backup_cmd = sprintf("mysqldump --host=%s --port=%s --user=%s --password=%s %s | gzip > %s", 
                             $xref_attr->{species}->{host},
                             $xref_attr->{species}->{port},
                             $xref_attr->{species}->{user},
                             $xref_attr->{species}->{password},
                             $xref_attr->{species}->{dbname},
                             $backup);
    &write_log("Running: $backup_cmd\n");
    system($backup_cmd) and die "Could not backup core database after parsing\n";
  }
}

#####################################################
sub xref_mapping {
  my @species = @_;

  foreach my $spe (@species) {
    my $xref_conf = "$PARASITE_CONF/xref_mapping.$spe.input";
    my $xref_attr = &parse_xref_inputconf($xref_conf);
    my $mapping_dir = $xref_attr->{xref}->{dir};

    if (not -d $mapping_dir) {
      mkpath($mapping_dir, { verbose => 0, mode => 0775 } );
    } else {
      foreach my $file (glob("$mapping_dir/*.*")) {
        unlink $file;
      }
    }

    my $cmd = sprintf("cd %s/ensembl/misc-scripts/xref_mapping && perl xref_mapper.pl -file %s > %s 2>&1",
                      $ENSEMBL_CVS_ROOT_DIR,
                      $xref_conf,
                      "$LOG_DIR/mapping.$spe.WS${WORMBASE_VERSION}.out");

    &write_log("Running: $cmd\n");
    system($cmd) and die "Command failure: $cmd\n";

    
    my $backup = sprintf("%s/%s.post_mapping.sql.gz", 
                         $BACKUP_DIR,
                         $xref_attr->{xref}->{dbname});
    my $backup_cmd = sprintf("mysqldump --host=%s --port=%s --user=%s --password=%s %s | gzip > %s", 
                             $xref_attr->{xref}->{host},
                             $xref_attr->{xref}->{port},
                             $xref_attr->{xref}->{user},
                             $xref_attr->{xref}->{password},
                             $xref_attr->{xref}->{dbname},
                             $backup);
    &write_log("Running: $backup_cmd\n");
    system($backup_cmd) and die "Could not backup xref_database after parsing\n";
  }
}

#####################################################
sub xref_loading {
  my @species = @_;

  foreach my $spe (@species) {
    my $xref_conf = "$PARASITE_CONF/xref_mapping.$spe.input";

    my $cmd = sprintf("cd %s/ensembl/misc-scripts/xref_mapping && perl xref_mapper.pl -file %s -upload > %s 2>&1",
                      $ENSEMBL_CVS_ROOT_DIR,
                      $xref_conf,
                      "$LOG_DIR/loading.$spe.WS${WORMBASE_VERSION}.out");
    system($cmd) and die "Failed to run xref_mapping for $spe\n";    
  }
}

#####################################################
sub xref_cleanup {

  foreach my $spe (@species) {
    my $xref_conf = "$PARASITE_CONF/xref_mapping.$spe.input";
    my $xref_attr = &parse_xref_inputconf($xref_conf);

    my $summary_sql = "SELECT db_name, count(*) " 
        . "FROM object_xref, xref, external_db "
        . "WHERE object_xref.xref_id = xref.xref_id "
        . "AND xref.external_db_id = external_db.external_db_id "
        . "GROUP BY db_name";
    
    my $log_file = "$LOG_DIR/xrefcheckandclean.$spe.WS${WORMBASE_VERSION}.out";
    open(my $log_fh, ">$log_file") or die "Could not open $log_file for writing\n";

    # check and cleanup
    my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
      -host => $xref_attr->{species}->{host},
      -port => $xref_attr->{species}->{port},
      -user => $xref_attr->{species}->{user},
      -pass => $xref_attr->{species}->{password},
      -dbname => $xref_attr->{species}->{dbname});

    #
    # Summary before cleanup 
    #
    print $log_fh "Object Xref summary before cleanup:\n\n";

    printf $log_fh "%-20s\t%10s\n", "DBNAME", "COUNT";
    my $sth = $dba->dbc->prepare($summary_sql);
    $sth->execute();
    while( my $row = $sth->fetchrow_arrayref ) {
      printf $log_fh "%-20s\t%10s\n", @$row;
    }
    $sth->finish;


    #
    # Cleanup 1 -  delete UniProt_gn object_xref
    #
    $dba->dbc->do("DELETE object_xref.*, dependent_xref.* " .
                  "FROM  object_xref, dependent_xref, xref, external_db " .
                  "WHERE external_db.external_db_id = xref.external_db_id " .
                  "AND xref.xref_id = object_xref.xref_id " .
                  "AND object_xref.object_xref_id = dependent_xref.object_xref_id " .
                  "AND db_name = 'Uniprot_gn'");


    #
    # Summary after cleanup 
    #
    print $log_fh "\n\nObject Xref summary after cleanup:\n\n";

    printf $log_fh "%-20s\t%10s\n", "DBNAME", "COUNT";
    my $sth = $dba->dbc->prepare($summary_sql);
    $sth->execute();
    while( my $row = $sth->fetchrow_arrayref ) {
      printf $log_fh "%-20s\t%10s\n", @$row;
    }
    $sth->finish;
  }

}

######################################################
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
  
  my %vhosts = ( STAGING => ($PARASITE_VERSION % 2) ? "mysql-ps-staging-1" : "mysql-ps-staging-2",
                 PROD    => "mysql-ps-prod",
                 PAN     => "mysql-pan-1" );
  
  push @replacements,  ( { before => "ENSEMBL_VERSION", 
                           after => $ENSEMBL_VERSION },
                         { before => "PARASITE_VERSION", 
                           after => $PARASITE_VERSION },
                         { before => "WORK_DIR", 
                           after => $WORK_DIR },
                         { before => "WORMBASE_VERSION", 
                           after => $WORMBASE_VERSION },
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

  my $wsrel = "WS".$WORMBASE_VERSION;
  
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

###########################################
sub parse_xref_inputconf {
  my ($input_conf) = @_;
    
  my ($header, %xref_attr);
  open(my $fh, $input_conf) or die "Could not open $input_conf for reading\n";

  while(<$fh>) {
    /^(\S+)$/ and do {
      my ($k, $v) = split(/=/, $1);
      
      if (not $header) {
        $header = $k;
        if (defined $v) {
          $xref_attr{$header}->{$k} = $v;
        }
        next;
      }
      
      
      if (defined $v) {
        $xref_attr{$header}->{$k} = $v;
      } else {
        $xref_attr{$header} = 1;
      }
    };
      
    /^\s*$/ and do {
      $header = "";
      next;
    };
  }

  return \%xref_attr;
}

###########################################
sub write_log {
  my ($msg) = @_;

  $verbose and print STDERR "$msg\n";
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
  meta.genebuild.version: WSWORMBASE_VERSION
caenorhabditis_elegans:
  seleno: WBGene00015553
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WSWORMBASE_VERSION/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WSWORMBASE_VERSION.annotations.gff3.gz
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    password: STAGINGPASSRW
    dbname: caenorhabditis_elegans_core_ENSEMBL_VERSION_WORMBASE_VERSION
caenorhabditis_briggsae:
  seleno: WBGene00028139
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WSWORMBASE_VERSION/species/c_briggsae/PRJNA10731/c_briggsae.PRJNA10731.WSWORMBASE_VERSION.annotations.gff3.gz  
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    password: STAGINGPASSRW
    dbname: caenorhabditis_briggsae_core_ENSEMBL_VERSION_WORMBASE_VERSION
caenorhabditis_brenneri:
  seleno: WBGene00158831
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WSWORMBASE_VERSION/species/c_brenneri/PRJNA20035/c_brenneri.PRJNA20035.WSWORMBASE_VERSION.annotations.gff3.gz
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    password: STAGINGPASSRW
    dbname: caenorhabditis_brenneri_core_ENSEMBL_VERSION_WORMBASE_VERSION
caenorhabditis_remanei:
  seleno: WBGene00068657
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WSWORMBASE_VERSION/species/c_remanei/PRJNA53967/c_remanei.PRJNA53967.WSWORMBASE_VERSION.annotations.gff3.gz
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    password: STAGINGPASSRW
    dbname: caenorhabditis_remanei_core_ENSEMBL_VERSION_WORMBASE_VERSION
caenorhabditis_japonica:
  seleno: WBGene00122465
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WSWORMBASE_VERSION/species/c_japonica/PRJNA12591/c_japonica.PRJNA12591.WSWORMBASE_VERSION.annotations.gff3.gz 
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    password: STAGINGPASSRW
    dbname: caenorhabditis_japonica_core_ENSEMBL_VERSION_WORMBASE_VERSION
pristionchus_pacificus:  
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WSWORMBASE_VERSION/species/p_pacificus/PRJNA12644/p_pacificus.PRJNA12644.WSWORMBASE_VERSION.annotations.gff3.gz 
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    password: STAGINGPASSRW
    dbname: pristionchus_pacificus_prjna12644_core_PARASITE_VERSION_ENSEMBL_VERSION_WORMBASE_VERSION
brugia_malayi:
  seleno: WBGene00222286
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WSWORMBASE_VERSION/species/b_malayi/PRJNA10729/b_malayi.PRJNA10729.WSWORMBASE_VERSION.annotations.gff3.gz
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    password: STAGINGPASSRW
    dbname: brugia_malayi_prjna10729_core_PARASITE_VERSION_ENSEMBL_VERSION_WORMBASE_VERSION
onchocerca_volvulus:
  seleno: WBGene00241445
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WSWORMBASE_VERSION/species/o_volvulus/PRJEB513/o_volvulus.PRJEB513.WSWORMBASE_VERSION.annotations.gff3.gz 
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    password: STAGINGPASSRW
    dbname: onchocerca_volvulus_prjeb513_core_PARASITE_VERSION_ENSEMBL_VERSION_WORMBASE_VERSION
strongyloides_ratti:
  gff3: /nfs/ftp/pub/databases/wormbase/releases/WSWORMBASE_VERSION/species/s_ratti/PRJEB125/s_ratti.PRJEB125.WSWORMBASE_VERSION.annotations.gff3.gz
  core_database:
    host: STAGINGHOST
    port: STAGINGPORT
    user: STAGINGUSERRW
    password: STAGINGPASSRW
    dbname: strongyloides_ratti_prjeb125_core_PARASITE_VERSION_ENSEMBL_VERSION_WORMBASE_VERSION
END_WORM_LITE_TEMPLATE

BEGIN_STAGING_REGISTRY_TEMPLATE
use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->no_version_check(1);
Bio::EnsEMBL::Registry->no_cache_warnings(1);
{
  Bio::EnsEMBL::Registry->load_registry_from_url('STAGINGURLRW/ENSEMBL_VERSION');

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

Bio::EnsEMBL::Registry->load_registry_from_url('STAGINGURLRO/ENSEMBL_VERSION');
Bio::EnsEMBL::Registry->load_registry_from_url('PRODURLRW/ENSEMBL_VERSION');

Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    -host => 'PANHOST',
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
dbname=xref_parasite_SPECIES_WORMBASE_VERSION
dir=WORK_DIR/xrefs/mapping/SPECIES/WBPSPARASITE_VERSION

species=SPECIES
taxon=wormbase
host=STAGINGHOST
port=STAGINGPORT
user=STAGINGUSERRW
password=STAGINGPASSRW
dbname=COREDBNAME
dir=WORK_DIR/xrefs/mapping/SPECIES/WBPSPARASITE_VERSION

farm
queue=production-rh6
exonerate=/nfs/panda/ensemblgenomes/external/exonerate-2/bin/exonerate

END_XREF_INPUT_TEMPLATE

__END__










