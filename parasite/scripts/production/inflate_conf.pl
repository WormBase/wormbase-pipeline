#!/usr/bin/env perl

use YAML;
use ProductionMysql;

my $conf = Load(do { local $/; <> }) or die "Usage: $0 <conf file>";

$conf->{generics} //= {
  taxonomy_database => ProductionMysql->new("mysql-ps-prod-1")->conn("ncbi_taxonomy_parasite"),
  production_database => ProductionMysql->new("mysql-ps-prod-1")->conn("ensembl_production_parasite"),
  cvsdir => $ENV{ENSEMBL_CVS_ROOT_DIR} 
};

for my $k (keys %{$conf}){
  next if $k eq "generics";
  my $species_conf = $conf->{$k};
  my ( $spe, $cies, $bioproject ) = split "_", $k, 3;
  my $species = "${spe}_${cies}";
  my $staging_conn = ProductionMysql->staging_writable->conn;
  for my $conn_detail (keys %{$staging_conn}){
    $species_conf->{core_database}->{$conn_detail} //= $staging_conn->{$conn_detail};
  }  
  $species_conf->{core_database}->{dbname} //= 
   "${species}_${bioproject}_core_$ENV{PARASITE_VERSION}_$ENV{ENSEMBL_VERSION}_1";

  next unless $species_conf->{meta};
  
  my $taxon_id = $species_conf->{taxon_id};  
  my $assembly_version = $species_conf->{assembly_version};
  my $assembly_accession = $species_conf->{meta}->{"assembly.accession"};
  my $provider_name = $species_conf->{meta}->{"provider.name"};
  my $provider_url = $species_conf->{meta}->{"provider.url"};
  my $biosample = $species_conf->{meta}->{"species.biosample"};
  
  next unless $taxon_id and $provider_name and $provider_url and $biosample;
  
  $species_conf->{meta}->{"assembly.coverage_depth"} //= "medium";
  $species_conf->{meta}->{"species.db_name"} //= "${species}_${bioproject}";
  $species_conf->{meta}->{"species.display_name"} //=
    sprintf( "%s %s (%s)", ucfirst($spe), $cies, uc($bioproject) );
  $species_conf->{meta}->{"species.division"} //= "EnsemblParasite";
  $species_conf->{meta}->{"species.alias"} //=
    sprintf( "%s %s", ucfirst($spe), $cies );
  $species_conf->{meta}->{"species.bioproject_id"} //= uc($bioproject);
  $species_conf->{meta}->{"species.ftp_genome_id"} //= uc($bioproject);
  $species_conf->{meta}->{"species.production_name"} //=
    "${species}_${bioproject}";
  $species_conf->{meta}->{"species.scientific_name"} //=
    sprintf( "%s %s", ucfirst($spe), $cies );
  $species_conf->{meta}->{"species.species_taxonomy_id"} //= $taxon_id;
  $species_conf->{meta}->{"species.taxonomy_id"} //= $taxon_id;
  $species_conf->{meta}->{"species.url"} //= join( "_", ucfirst($spe), $cies, $bioproject);

  next unless $assembly_version and $assembly_accession;
  $species_conf->{meta}->{"assembly.name"} //= $assembly_version;
  $species_conf->{meta}->{"assembly.default"} //= $assembly_version;
  $species_conf->{meta}->{"assembly.default"} =~ s/\s+/_/g; # Ensembl convention too
}

print Dump($conf);
