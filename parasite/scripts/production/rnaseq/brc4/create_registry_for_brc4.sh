#!/usr/bin/bash

for VAR in  RNASEQ_CONF  \
            RNASEQ_STAGING_MYSQL  \
            PARASITE_ENSEMBL_VERSION \
            PARASITE_STAGING_MYSQL
do
   if [ -z $(eval "echo \$$VAR") ]; then
      echo "$VAR must be set in your environment before running $0"
      exit 254
   fi
done

printf \
"use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;

Bio::EnsEMBL::Registry->load_registry_from_db(
  -host       => '$($PARASITE_STAGING_MYSQL host)',
  -port       => '$($PARASITE_STAGING_MYSQL port)',
  -user       => '$($PARASITE_STAGING_MYSQL user)',
  -db_version => '${PARASITE_ENSEMBL_VERSION}',
);

Bio::EnsEMBL::DBSQL::OntologyDBAdaptor->new(
  -species => 'multi',
  -dbname => 'ensembl_ontology_${PARASITE_ENSEMBL_VERSION}',
  -group => 'ontology',
  -host => '$($PARASITE_STAGING_MYSQL host)',
  -port => '$($PARASITE_STAGING_MYSQL port)',
  -user => '$($PARASITE_STAGING_MYSQL user)',
);

1;
"

