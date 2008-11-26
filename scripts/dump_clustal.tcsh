/software/worm/postgres/bin/pg_dump -d -x -O -t clustal -h deskpro16391.dynamic.sanger.ac.uk -U clx clx |perl ~/wormbase/scripts/p2msql.pl |bzip2 -9 -c >! ~wormpub/BUILD/autoace/wormpep_clw.sql.bz2
