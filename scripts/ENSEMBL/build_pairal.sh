cd /nfs/acari/wormpipe/ensembl/ensembl-compara/scripts/pipeline
mysql -h ecs1f -u wormadmin -pworms -e 'DROP DATABASE IF EXISTS worm_compara_blastz;CREATE DATABASE worm_compara_blastz;'
cat ~wormpipe/ensembl/ensembl-hive/sql/tables.sql ~wormpipe/ensembl/ensembl-compara/sql/table.sql ~wormpipe/ensembl/ensembl-compara/sql/pipeline-tables.sql |mysql -h ecs1f -u wormadmin -pworms worm_compara_blastz
~wormpipe/ensembl/ensembl-compara/scripts/pipeline/comparaLoadGenomes.pl -conf compara-hive-pairaligner.conf
~wormpipe/ensembl/ensembl-compara/scripts/pipeline/loadPairAlignerSystem.pl -conf compara-hive-pairaligner.conf
