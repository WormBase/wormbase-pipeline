# This script was created to perform this patch: https://github.com/Ensembl/staging-patches/pull/492
# While the xrefs for the WormBase species were correct the display xref ids / description
# for genes/transcript had not been set. As a result, users were not able to view the human-readable
# names and the description of the genes/transcripts. This is why it was really important to fix this.
# This script takes 1) the host and user details from the source and target db servers and
# 2) a tab delimited file as its input with source db names on the first column and target db names on the second column
# It then uses fix_display_xrefs.pl to output the SQL commands necessary to fix the display_xref_ids missing from
# from the target DBs. It's also using the copy_descriptions.sh script to print the SQL commands need to copy
# gene/trascript descriptions from the source DBs to the target DBs.

#Input
host1=$1
host2=$2
user2=$3
port2=$4
patchout=$5
db_correspondance_file=$6

#Input Examples:
# host1=mysql-ps-staging-2
# host2=mysql-ens-sta-3
# user2=ensro
# port2=4160
# patchout=/homes/digri/e106_wb_patches This is the output files where the SQL commands will be written.
# db_correspondance_file=~/handovered_dbs_sta3_3.txt

# db_correspondance_file EXAMPLE:
#brugia_malayi_core_53_106_279   brugia_malayi_prjna10729_core_16_101_279
#caenorhabditis_brenneri_core_53_106_279 caenorhabditis_brenneri_prjna20035_core_16_101_279
#caenorhabditis_briggsae_core_53_106_279 caenorhabditis_briggsae_prjna10731_core_16_101_279
#caenorhabditis_elegans_core_53_106_279  caenorhabditis_elegans_prjna13758_core_16_101_279
#caenorhabditis_japonica_core_53_106_279 caenorhabditis_japonica_prjna12591_core_16_101_279
#caenorhabditis_remanei_core_53_106_279  caenorhabditis_remanei_prjna53967_core_16_101_279
#onchocerca_volvulus_core_53_106_279     onchocerca_volvulus_prjeb513_core_16_101_279
#pristionchus_pacificus_core_53_106_279  pristionchus_pacificus_prjna12644_core_16_101_279
#strongyloides_ratti_core_53_106_279     strongyloides_ratti_prjeb125_core_16_101_279
#trichuris_muris_core_53_106_279 trichuris_muris_prjeb126_core_16_101_279

export PERL5LIB=${ENSEMBL_CVS_ROOT_DIR}/ensembl/misc-scripts/xref_mapping:$PERL5LIB
cat ${db_correspondance_file} | while read dbname2 dbname1;
  do if [[ $($host2 -Ne "SHOW DATABASES" | grep $dbname2) ]];
    then
      rm $patchout/$dbname2.sql
      touch $patchout/$dbname2.sql
      printf "use ${dbname2};\n" > $patchout/${dbname2}.sql
      perl $WORM_CODE/parasite/scripts/production/xref/fix_display_xrefs.pl ${host2} ${port2} ${user2} ${dbname2} >> $patchout/${dbname2}.sql
      printf "UPDATE transcript tr JOIN xref ON tr.display_xref_id=xref.xref_id SET tr.display_xref_id = NULL WHERE xref.display_label=tr.stable_id;\n" >> $patchout/$dbname2.sql
      sh ${WORM_CODE}/parasite/scripts/production/xref/copy_descriptions.sh $host1 $dbname1 $host2 $dbname2 >> $patchout/${dbname2}.sql
    else
      echo "ERROR. ${dbname2} does not exist in ${host2}"
  fi;
  done
