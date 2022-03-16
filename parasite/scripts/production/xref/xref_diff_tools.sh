diff_query_results(){
 sql=$1
 filter=$2
 species=$3
 core_db=$(perl -MProductionMysql -E 'say ProductionMysql->staging->core_db(@ARGV)' "$species" )
 previous_core_db=$(perl -MProductionMysql -E 'say ProductionMysql->previous_staging->core_db(@ARGV)' "$species" )
 [ "$core_db" ] || return
 [ "$previous_core_db" ] || return
 printf "Old:\n"
 $PREVIOUS_PARASITE_STAGING_MYSQL $previous_core_db -Ne "$sql" | perl -pe "$filter" | sort
 printf "\nNew:\n"
 $PARASITE_STAGING_MYSQL $core_db -Ne "$sql" | perl -pe "$filter" | sort
}
join_query_results(){
 sql=$1
 species=$2
 core_db=$(perl -MProductionMysql -E 'say ProductionMysql->staging->core_db(@ARGV)' "$species" )
 previous_core_db=$(perl -MProductionMysql -E 'say ProductionMysql->previous_staging->core_db(@ARGV)' "$species" )
 [ "$core_db" ] || return
 [ "$previous_core_db" ] || return
 join -t $'\t' -11 -21 \
   <( $PREVIOUS_PARASITE_STAGING_MYSQL $previous_core_db -Ne "$sql" | sort -k1,1) \
   <( $PARASITE_STAGING_MYSQL $core_db -Ne "$sql" | sort -k1,1) \
   | perl -nE 'chomp ; my @F = split "\t"; say unless @F[1] eq @F[2]'
}
diff_xref_table(){
  diff_query_results \
    'select db_name, count(*) from xref join external_db using (external_db_id) group by external_db_id' \
    '' \
    "$@"
}
diff_object_xref_table(){
  diff_query_results \
    'select db_name, count(*) from xref join external_db using (external_db_id) join object_xref using (xref_id) group by external_db_id' \
    '' \
    "$@"
}
 diff_gene_descriptions(){
  diff_query_results \
    'select stable_id, description from gene ' \
    's/\[.*\]//' \
    "$@"
}
join_gene_descriptions(){
  join_query_results \
    'select stable_id, description from gene ' \
    "$@"
}
