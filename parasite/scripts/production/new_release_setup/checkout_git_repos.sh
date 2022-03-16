#!/usr/bin/env bash

if [ -z $ENSEMBL_VERSION ]; then
  echo "Need to define ENSEMBL_VERSION in the environment - have you loaded the parasite production module file?"
  exit
fi

if [ -z $ENSEMBL_CVS_ROOT_DIR ]; then
  echo "Need to define ENSEMBL_CVS_ROOT_DIR in the environment - have you loaded the parasite production module file?"
  exit
fi

echo "Making ${ENSEMBL_CVS_ROOT_DIR}"
mkdir -pv ${ENSEMBL_CVS_ROOT_DIR}

checkout_ensembl_branch(){
 repo_name=$1
 branch_name=$2
 repo_owner=${3:-ensembl}
 where_to_clone=${ENSEMBL_CVS_ROOT_DIR}/${repo_name}
 if [ -d $where_to_clone ] ; then
    echo "$where_to_clone already exists - removing"
    rm -rf $where_to_clone
 fi
 git clone --single-branch -b $branch_name --depth 1 https://github.com/${repo_owner}/${repo_name}.git ${where_to_clone}
}

checkout_ensembl_branch ensembl-analysis dev/hive_master
checkout_ensembl_branch ensembl-pipeline master
checkout_ensembl_branch ensembl-taxonomy master
checkout_ensembl_branch ensembl-funcgen main
checkout_ensembl_branch ensembl release/${ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-compara release/${ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-io release/${ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-orm release/${ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-production release/${ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-variation release/${ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-vep release/${ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-metadata release/${ENSEMBL_VERSION}
checkout_ensembl_branch ensj-healthcheck release/${ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-datacheck master
checkout_ensembl_branch ensembl-py main
checkout_ensembl_branch ensembl-hive version/2.6
checkout_ensembl_branch parasite-static master EnsemblGenomes
checkout_ensembl_branch eg-web-common master EnsemblGenomes
checkout_ensembl_branch ensembl-production-imported trunk
checkout_ensembl_branch ensembl-production-imported-private trunk
checkout_ensembl_branch ensembl-taxonomy main
checkout_ensembl_branch GIFTS main
