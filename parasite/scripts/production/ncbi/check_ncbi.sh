#!/usr/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
rm -rf $DIR/tmp
mkdir $DIR/tmp
$DIR/assemblies_for_taxon.pl 6157 6231 \
  | tee $DIR/tmp/assemblies_for_taxon.out \
  | $DIR/doc_for_assembly.pl \
  | tee $DIR/tmp/doc_for_assembly.out\
  | $DIR/compare_docs_against_current_staging.pl \
  | tee $DIR/tmp/compare_docs_against_current_staging.out \
  | perl -MHash::Merge -MYAML -e 'local $/; my @docs = Load(<>); print Dump (Hash::Merge::merge @docs);' \
  | tee $DIR/tmp/result.out
