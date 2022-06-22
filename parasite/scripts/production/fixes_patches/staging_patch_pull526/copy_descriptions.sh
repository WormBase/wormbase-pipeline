host1=$1
dbname1=$2
host2=$3
dbname2=$4

#This script prints the SQL commands needing to fix gene/transcript descriptions
#dbname2 by copying the gene/transcript descriptions from dbname1.
#Usage sh.copy_descriptions host1 port1 user1 dbname1 host2 port2 user2 dbname2

if [[ $($host1 $dbname1 -Ne "SELECT description FROM gene WHERE description IS NOT NULL;") ]];
  then
    $host1 $dbname1 -Ne "SELECT stable_id,description FROM gene WHERE description IS NOT NULL;" | while read gid desc;
    do
      printf "UPDATE gene SET description=\"${desc}\" WHERE stable_id=\"${gid}\";\n";
    done;
fi

if [[ $($host1 $dbname1 -Ne "SELECT description FROM transcript WHERE description IS NOT NULL;") ]];
  then
    $host1 $dbname1 -Ne "SELECT stable_id,description FROM transcript WHERE description IS NOT NULL;" | while read gid desc;
    do
      printf "UPDATE transcript SET description=\"${desc}\" WHERE stable_id=\"${gid}\";\n";
    done;
fi