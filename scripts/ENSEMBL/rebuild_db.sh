cd /nfs/acari/wormpipe/ensembl/ensembl-compara/scripts/pipeline
mysql -h ecs1f -u wormadmin -pworms -e 'DROP DATABASE worm_compara;CREATE DATABASE worm_compara;'
cat ~wormpipe/ensembl/ensembl-hive/sql/tables.sql ~wormpipe/ensembl/ensembl-compara/sql/table.sql ~wormpipe/ensembl/ensembl-compara/sql/pipeline-tables.sql |mysql -h ecs1f -u wormadmin -pworms worm_compara
~wormpipe/ensembl/ensembl-compara/scripts/pipeline/comparaLoadGenomes.pl -conf my_compara.conf
~wormpipe/ensembl/ensembl-compara/scripts/pipeline/loadHomologySystem.pl -conf my_compara.conf
