#!/usr/bin/bash

for VAR in  HANDOVER_CONF  \
            HANDOVER_STAGING_MYSQL  \
            ENSEMBL_VERSION         \
            EG_VERSION
do
   if [ -z $(eval "echo \$$VAR") ]; then
      echo "$VAR must be set in your environment before running $0"
      exit 254
   fi
done

printf \
"use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Production::DBSQL::DBAdaptor;

Bio::EnsEMBL::Registry->no_version_check(1);
Bio::EnsEMBL::Registry->no_cache_warnings(1);
{
  Bio::EnsEMBL::Registry->load_registry_from_url('$($HANDOVER_STAGING_MYSQL-w details url)${ENSEMBL_VERSION}');

  Bio::EnsEMBL::MetaData::DBSQL::MetaDataDBAdaptor->new(
    -host    => '$($HANDOVER_STAGING_MYSQL host)',
    -port    => '$($HANDOVER_STAGING_MYSQL port)',
    -user    => '$($HANDOVER_STAGING_MYSQL user)',
    -dbname  => 'ensembl_metadata_${ENSEMBL_VERSION}',
    -species => 'multi',
    -group   => 'metadata',
  );

  Bio::EnsEMBL::Production::DBSQL::DBAdaptor->new(
    -host    => '$(meta1 host)',
    -port    => '$(meta1 port)',
    -user    => '$(meta1 user)',
    -dbname  => 'ensembl_production',
    -species => 'multi',
    -group   => 'production'
  );
}
1;
"

