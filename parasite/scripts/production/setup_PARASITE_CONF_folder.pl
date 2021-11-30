#!/usr/bin/env perl
use strict;
use File::Path qw(mkpath);
use ProductionMysql;
use Try::Tiny;

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
my $PARASITE_STAGING     = $ENV{PARASITE_STAGING_MYSQL};

my $templates = &read_templates();
  
mkpath($PARASITE_CONF, { verbose => 0, mode => 0775 }) if not -e $PARASITE_CONF; 

&write_config($templates->{STAGING_REGISTRY}, "$PARASITE_CONF/staging_pipelines.registry.pm" );
&write_config($templates->{COMPARA_REGISTRY}, "$PARASITE_CONF/compara.registry.pm" );

#####################################################
sub find_ensrw_db {
  my ($db) = @_;
  try {
    ProductionMysql->new("${db}-w")->conn;
    return "${db}-w";
  }
  catch {
    return "${db}-ensrw";
  }
}

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
  
  my %vhosts = ( STAGING => $PARASITE_STAGING, 
                 PROD    => "mysql-ps-prod-1",
                 PAN     => "mysql-eg-pan-prod" );
  
  push @replacements,  ( { before => "ENSEMBL_VERSION", 
                           after => $ENSEMBL_VERSION },
                         { before => "PARASITE_VERSION", 
                           after => $PARASITE_VERSION },
  );
  
  foreach my $vhost_key (keys %vhosts) {
    my $vhost_name = $vhosts{$vhost_key};

    my $con_ro = ProductionMysql->new($vhost_name)->conn;
    my $con_rw = ProductionMysql->new(find_ensrw_db(${vhost_name}))->conn;

    my $url_ro = ProductionMysql->new($vhost_name)->url;
    my $url_rw = ProductionMysql->new(find_ensrw_db(${vhost_name}))->url;

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

__DATA__

BEGIN_STAGING_REGISTRY_TEMPLATE
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Production::DBSQL::DBAdaptor;

Bio::EnsEMBL::Registry->no_version_check(1);
Bio::EnsEMBL::Registry->no_cache_warnings(1);
{
  Bio::EnsEMBL::Registry->load_registry_from_url('STAGINGURLRW/ENSEMBL_VERSION');

  Bio::EnsEMBL::Production::DBSQL::DBAdaptor->new(
    -host    => 'PRODHOST',
    -port    => 'PRODPORT',
    -user    => 'PRODUSERRO',
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
    -host => 'PRODHOST',
    -port => 'PRODPORT',
    -user => 'PRODUSERRW',
    -pass => 'PRODPASSRW',
    -species => 'ensembl_compara_master_parasite',
    -dbname  => 'ensembl_compara_master_parasite',
);

Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    -host => 'PRODHOST',
    -user => 'PRODUSERRO',
    -port => 'PRODPORT',             
    -species => 'ncbi_taxonomy_parasite',
    -dbname => 'ncbi_taxonomy_parasite',
);
END_COMPARA_REGISTRY_TEMPLATE

__END__










