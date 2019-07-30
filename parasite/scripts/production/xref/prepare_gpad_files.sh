

get_core_db_and_taxon_list () { 
  for db in $(  $PARASITE_STAGING_MYSQL -Ne "show databases like \"%core_${PARASITE_VERSION}%\"") ; do 
    echo -n "$db" "" 
    $PARASITE_STAGING_MYSQL $db -Ne 'select meta_value from meta where meta_key = "species.taxonomy_id"' 
  done
}
format_taxons_into_filter () { 
  perl -ne 'chomp; my ($core_db, $taxon ) = split /\s+/; print "taxon:$taxon\t\n" if $taxon' "$@" | sort -u 
}

SRC="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [ ! "$1" ] ; then
  >&2 echo "Usage: $0 <target directory>"
  exit 1
fi
DIR=$1

mkdir -pv $DIR/.tmp
rm -vf $DIR/*

test -f "$DIR/.tmp/core_db_and_taxon_list.tsv" || get_core_db_and_taxon_list | pv -ltrbcN "Get taxon list from core databases" > $DIR/.tmp/core_db_and_taxon_list.tsv

cat $DIR/.tmp/core_db_and_taxon_list.tsv \
  | perl -MYAML -e 'my %h; 
      while(<>){
        my ($core_db, $taxon) = split; 
        $h{$taxon}//=[]; 
        $core_db =~ s/_core_.*//; 
        push @{$h{$taxon}}, $core_db; 
      } ; print Dump(\%h);' \
  > $DIR/.tmp/core_db_per_taxon.yaml

if [ -d "/nfs/ftp" ] ; then
  cp -v /nfs/ftp/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz $DIR/.tmp/goa_uniprot_all.gaf.gz
else 
  wget -c ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz -O $DIR/.tmp/goa_uniprot_all.gaf.gz
fi

pv -cN "Pipe out compressed Uniprot file" $DIR/.tmp/goa_uniprot_all.gaf.gz \
  | zgrep -F -f <( format_taxons_into_filter $DIR/.tmp/core_db_and_taxon_list.tsv ) \
  | pv -lcN "Pick up lines with ParaSite taxons" \
  | $SRC/gaf_to_gpad.pl \
  | pv -lcN "Convert into gpad format, and save" \
  | $SRC/separate_gpad_by_species.pl - $DIR/.tmp/core_db_per_taxon.yaml $DIR

num_files=$(ls -1 $DIR | wc -l)
num_core_dbs=$( wc -l < $DIR/.tmp/core_db_and_taxon_list.tsv)
>&2 echo "Finished! Created files: $num_files for $num_core_dbs core databases" 

if [ "$num_files" -ne "$num_core_dbs" ] ; then
  >&2 echo "The two numbers are not the same. Blimey! "
  >&2 echo "Diff directory content vs core dbs: "
  >&2 diff \
    <( ls -1 $DIR | perl -pe 's/annotations_ensembl-(.*).gpa/$1/' | sort )  \
    <( perl -pe 's/_core.*//' < $DIR/.tmp/core_db_and_taxon_list.tsv | sort )
  exit 1
else
  rm -v $DIR/.tmp/**
fi
