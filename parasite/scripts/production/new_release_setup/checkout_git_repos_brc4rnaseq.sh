#!/usr/bin/env bash

if [ -z RNASEQ_ENSEMBL_VERSION ]; then
  echo "Need to define RNASEQ_ENSEMBL_VERSION_ENSEMBL_VERSION in the environment - have you loaded the parasite_brc4rnaseq_106 module file?"
  exit
fi

if [ -z RNASEQ_CVS_ROOT_DIR ]; then
  echo "Need to define RNASEQ_CVS_ROOT_DIR in the environment - have you loaded the parasite_brc4rnaseq_106 module file?"
  exit
fi

echo "Making ${RNASEQ_CVS_ROOT_DIR}"
mkdir -pv ${RNASEQ_CVS_ROOT_DIR}

checkout_ensembl_branch(){
 repo_name=$1
 branch_name=$2
 repo_owner=${3:-ensembl}
 where_to_clone=${RNASEQ_CVS_ROOT_DIR}/${repo_name}
 if [ -d $where_to_clone ] ; then
    echo "$where_to_clone already exists - removing"
    rm -rf $where_to_clone
 fi
 git clone --single-branch -b $branch_name --depth 1 https://github.com/${repo_owner}/${repo_name}.git ${where_to_clone}
}

checkout_ensembl_branch ensembl release/106
checkout_ensembl_branch ensembl-hive version/2.6
checkout_ensembl_branch ensembl-production-imported main
checkout_ensembl_branch ensembl-production-imported-private main
checkout_ensembl_branch ensembl-taxonomy main
checkout_ensembl_branch ensembl-io main


