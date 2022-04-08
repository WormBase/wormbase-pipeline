#!/usr/bin/bash

for VAR in  PARASITE_STAGING_MYSQL \
            PARASITE_VERSION \
            ENSEMBL_VERSION
do
   if [ -z $(eval "echo \$$VAR") ]; then
      echo "$VAR must be set in your environment before running $0"
      exit 254
   fi
done

printf \
"
#!/usr/bin/env perl

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;


"

VARIATION_GENOMES=( "caenorhabditis_elegans_prjna13758" "schistosoma_mansoni_prjea36577" )
for GENOME in ${VARIATION_GENOMES[@]}; do
  CORE_DB=$($PARASITE_STAGING_MYSQL -Ne "SHOW databases like \"%${GENOME}_core_${PARASITE_VERSION}_${ENSEMBL_VERSION}%\";" | head -n1)
  VARIATION_DB=$(echo $CORE_DB | sed 's/core/variation/g')
  printf \
  "#${GENOME}
Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -host    => '$($PARASITE_STAGING_MYSQL host)',
  -port    => '$($PARASITE_STAGING_MYSQL port)',
  -user    => '$($PARASITE_STAGING_MYSQL user)',
  -dbname  => '${CORE_DB}',
  -species => '${GENOME}'
);

Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(
  -host    => '$($PARASITE_STAGING_MYSQL-w host)',
  -port    => '$($PARASITE_STAGING_MYSQL-w port)',
  -user    => '$($PARASITE_STAGING_MYSQL-w user)',
  -pass    => '$($PARASITE_STAGING_MYSQL-w pass)',
  -dbname  => '${VARIATION_DB}',
  -species => '${GENOME}'
);

";
  done;

printf "1;"

