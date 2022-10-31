#!/bin/bash
usage()
{
    echo "usage: <command> options:<c|r|d>"
    printf "This script produces sql patches files for renaming gene/transcripts stable_id in the specified core db\n
            based on an input TSV file. It can optionally do tha same with gene descriptions\n"
    printf "c: REQUIRED - Core DB which you want to rename stable_ids in.\n\n"
    printf "r: REQUIRED - Path to a TSV file with two columns without headers.\ncolumn1=old_gene_stable_id, column2=new_gene_id.\n\n"
    printf "d: OPTIONAL - Path to a TSV file with two columns without headers.\ncolumn1=old_gene_stable_id, column2=description.\n\n"
}
while getopts h:c:r:d: flag
do
    case "${flag}" in
        c) core_db=${OPTARG};;
        r) rename_file=${OPTARG};;
        d) descriptions_file=${OPTARG};;
    esac
done


if [ -z "$core_db" ]; then usage; exit; fi
if [ -z "$rename_file" ]; then usage; exit; fi

cdb=${core_db};
rename_file=${rename_file};

$PARASITE_STAGING_MYSQL $cdb -Ne "SELECT stable_id FROM gene" > gene.stable_id
$PARASITE_STAGING_MYSQL $cdb -Ne "SELECT stable_id FROM gene WHERE display_xref_id IS NOT NULL" > gene.stable_id_xref_not_null
$PARASITE_STAGING_MYSQL $cdb -Ne "SELECT stable_id FROM transcript" > tr.stable_id

cut -f1 ${rename_file} | while read gene; do grep -w $gene gene.stable_id; done > gene.to_rename
if [[ $(cat gene.stable_id | sort | uniq -d | wc -l) -ne 0 ]];
  then echo "ERROR: ${cdb}'s gene table have duplicates in stable ids."; exit;
  else echo "INFO: No stable id clashes.. Good!";
fi;

cut -f1 ${rename_file} | while read gene; do grep -w $gene tr.stable_id; done > tr.to_rename
if [[ $(cat tr.stable_id | sort | uniq -d | wc -l) -ne 0 ]]; then echo "ERROR: ${cdb}'s transcript table have duplicates in stable ids."; exit; fi

#Datacheck
cat gene.to_rename | while read gstid;
  do $PARASITE_STAGING_MYSQL $cdb -Ne \
  "SELECT tr.stable_id FROM gene JOIN transcript tr USING (gene_id) WHERE gene.stable_id='${gstid}';";
done > tr.to_rename_from_db;
if [[ $(cat tr.to_rename | sort | uniq | wc -l) -ne $(cat tr.to_rename_from_db | sort | uniq | wc -l) ]];
  then echo "ERROR: You're missing transcripts (or having excess) in tr.to_rename"; exit;
fi

# Genes SQL
printf "INFO: Creating gene.sql.patch\n"

rm gene.sql.patch;
printf "INFO: Creating gene.sql.patch\n"
echo "USE ${cdb};" > gene.sql.patch;
cat ${rename_file} | while read -r old_gene new_gene;
  do
  if grep -q {new_gene} gene.stable_id;
    then printf "ERROR: new gene id ${new_gene} already exists in the DB. Exiting\n"; exit;
  fi;
  if grep -q ${old_gene} gene.stable_id;
    then printf "UPDATE gene SET stable_id='${new_gene}', display_xref_id=NULL WHERE stable_id='${old_gene}';\n" >> gene.sql.patch;
    else echo "ERROR: gene ${old_gene} is not in the gene table of ${cdb}"
  fi;
done;
## Datacheck
if [[ $(cat gene.sql.patch | tail -n +2 | sort | uniq | wc -l) -ne $(cat gene.to_rename | sort | uniq | wc -l) ]];
  then echo "ERROR: You're missing gene (or having excess) in gene.sql.patch"; exit;
fi
printf "  Done\n\n"

# Transcripts SQL
printf "INFO: Creating tr.sql.patch\n"
rm tr.sql.patch;
echo "USE ${cdb};" > tr.sql.patch;
cat ${rename_file} | while read -r old_gene new_gene;
  do awk -v pat="$old_gene" '$1 ~ pat' tr.stable_id | \
  while read old_tr;
    do new_tr=$(echo $old_tr | sed -s s,"${old_gene}","${new_gene}",g);
    if grep -q ${new_tr} tr.stable_id;
      then printf "ERROR: new gene id ${new_tr} already exists in the DB. Exiting\n"; exit;
    fi;
    printf "UPDATE transcript SET stable_id='${new_tr}', display_xref_id=NULL WHERE stable_id='${old_tr}';\n";
    done;
  done >> tr.sql.patch
## Datacheck
if [[ $(cat tr.sql.patch | tail -n +2 | sort | uniq | wc -l) -ne $(cat tr.to_rename_from_db | sort | uniq | wc -l) ]];
  then echo "ERROR: You're missing transcripts (or having excess) in tr.sql.patch"; exit;
fi
printf "  Done\n\n"

if [ -z ${descriptions_file+x} ];
  then echo "INFO: descriptions_file is unset. Not going to create a descriptions_xref_script.sql"; exit;
  else echo "INFO: descriptions_file is set. Going to create a descriptions_xref_script.sql...";
fi;

if [[ $(cat ${descriptions_file} | sort | uniq -d | wc -l) -ne 0 ]]; then echo "ERROR: ${description_file} have duplicates."; fi
if [[ $(cat ${descriptions_file} | cut -f1 | sort | uniq -d | wc -l) -ne 0 ]]; then echo "ERROR: ${description_file} genes have duplicates."; fi
if [[ $(cat ${descriptions_file} | sort | uniq -d | wc -l) -ne $(cat gene.to_rename | sort | uniq -d | wc -l) ]]; then echo "ERROR: ${description_file} have different number of lines compared to gene.to_rename."; fi

# Descriptions SQL
printf "INFO: Creating descriptions.sql.patch\n"
rm descriptions.sql.patch
cat ${rename_file} | while read -r old_gene new_gene;
  do awk -v pat="$old_gene" '$1 ~ pat' ${descriptions_file} | \
  while read -r og desc;
    do printf "UPDATE gene SET description='${desc}', display_xref_id=NULL WHERE stable_id='${new_gene}';\n";
    done;
  done > descriptions.sql.patch

if [[ $(cat descriptions.sql.patch | sort | uniq | wc -l) -ne $(cat gene.to_rename | sort | uniq | wc -l) ]];
  then echo "ERROR: You're missing descriptions (or having excess) in descriptions.sql.patch";
fi;
printf "  Done\n"
printf "INFO: The descriptions.sql.patch should be applied after the Xref pipeline.\n\n"
printf "Done. Exiting\n"
exit