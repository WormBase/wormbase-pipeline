rm -f  /nfs/wormpub/BUILD/autoace/wormpep_clw.sql.gz
/software/worm/postgres/bin/pg_dump --inserts -x -O -t clustal -h pgsrv1 -p 5436 -U worm_clw -t clustal worm_clw |perl ~wormpub/wormbase/scripts/p2msql.pl | gzip -9 -c > /nfs/wormpub/BUILD/autoace/MISC_OUTPUT/wormpep_clw.sql.gz
