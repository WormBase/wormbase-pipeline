#!/usr/bin/env bash

if [ -z $GENEBUILD_ENSEMBL_VERSION ]; then
  echo "Need to define GENEBUILD_ENSEMBL_VERSION in the environment - have you loaded the parasite production module file?"
  exit
fi

if [ -z $ ]; then
  echo "Need to define GENEBUILD_CVS_ROOT_DIR in the environment - have you loaded the parasite production module file?"
  exit
fi

echo "Making ${GENEBUILD_CVS_ROOT_DIR}"
mkdir -pv ${GENEBUILD_CVS_ROOT_DIR}

checkout_ensembl_branch(){
 repo_name=$1
 branch_name=$2
 repo_owner=${3:-ensembl}
 where_to_clone=${GENEBUILD_CVS_ROOT_DIR}/${repo_name}
 if [ -d $where_to_clone ] ; then
    echo "$where_to_clone already exists - removing"
    rm -rf $where_to_clone
 fi
 git clone --single-branch -b $branch_name --depth 1 https://github.com/${repo_owner}/${repo_name}.git ${where_to_clone}
}

checkout_ensembl_branch ensembl release/${GENEBUILD_ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-analysis ${ENSEMBL_ANALYSIS_BRANCH}
checkout_ensembl_branch ensembl-hive version/2.5
checkout_ensembl_branch ensembl-production release/${GENEBUILD_ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-killlist main
checkout_ensembl_branch ensembl-io release/${GENEBUILD_ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-variation release/${GENEBUILD_ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-funcgen release/${GENEBUILD_ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-taxonomy main
checkout_ensembl_branch ensembl-metadata release/${GENEBUILD_ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-datacheck release/${GENEBUILD_ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-orm release/${GENEBUILD_ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-compara release/${COMPARA_ENSEMBL_VERSION}
checkout_ensembl_branch ensembl-genes main



