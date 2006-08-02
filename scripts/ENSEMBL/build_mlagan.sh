cd /nfs/acari/wormpipe/ensembl/ensembl-compara/scripts/pipeline
mysql -h ecs1f -u wormadmin -pworms -e 'DROP DATABASE IF EXISTS worm_compara_lagan;CREATE DATABASE worm_compara_lagan;'
cat ~wormpipe/ensembl/ensembl-hive/sql/tables.sql ~wormpipe/ensembl/ensembl-compara/sql/table.sql ~wormpipe/ensembl/ensembl-compara/sql/pipeline-tables.sql |mysql -h ecs1f -u wormadmin -pworms worm_compara_lagan
~wormpipe/ensembl/ensembl-compara/scripts/pipeline/comparaLoadGenomes.pl -conf compara-hive-multiplealigner.conf
~wormpipe/ensembl/ensembl-compara/scripts/pipeline/loadMultipleAlignerSystem.pl -conf compara-hive-multiplealigner.conf
