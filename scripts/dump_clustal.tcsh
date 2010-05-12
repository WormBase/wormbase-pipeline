/software/worm/postgres/bin/pg_dump -d -x -O -t clustal -h pgsrv1 -p 5436 -U worm_clw -t clustal worm_clw |perl ~/wormbase/scripts/p2msql.pl |bzip2 -9 -c >! ~wormpub/BUILD/autoace/wormpep_clw.sql.bz2
