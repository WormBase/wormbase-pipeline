export GB_HOME=/nfs/production/flicek/ensembl/genebuild
export ENSEMBL_ROOT_DIR=/hps/software/users/ensembl/repositories/$USER
GB_REPO=/hps/software/users/ensembl/repositories/genebuild

if [[ -n "$LOCAL_PYENV" ]] && [[ -e "/hps/software/users/ensembl/ensw/latest/envs/minimal.sh" ]]; then
  . /hps/software/users/ensembl/ensw/latest/envs/minimal.sh
elif [[ -e "/hps/software/users/ensembl/ensw/latest/envs/essential.sh" ]]; then
  . /hps/software/users/ensembl/ensw/latest/envs/essential.sh
fi
export GB_SCRATCH=/hps/nobackup/flicek/ensembl/genebuild
export BLASTDB_DIR=$GB_SCRATCH/blastdb
export REPEATMODELER_DIR=$GB_SCRATCH/custom_repeat_libraries/repeatmodeler
export LSB_DEFAULTQUEUE="production"
if [[ -n "$LINUXBREW_HOME" ]];then
  if [[ -z "$WISECONFIGDIR" ]]; then
    export PATH="/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/opt/exonerate09/bin:$PATH"
  fi
  export BIOPERL_LIB="/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/opt/bioperl-169/libexec"
  export WISECONFIGDIR="/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/share/genewise"
  export GBLAST_PATH="/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin"
fi

if [[ -d "/nfs/production/panda/ensembl/production/ensemblftp/data_files" ]];then
  export FTP_DIR="/nfs/production/panda/ensembl/production/ensemblftp/data_files"
fi

export HIVE_EMAIL="$USER@ebi.ac.uk"
export ENSCODE=$ENSEMBL_ROOT_DIR
export GENEBUILDER_ID=0

# Tokens for different services
# webhooks for Slack to send notification to any channel
export SLACK_GENEBUILD='T0F48FDPE/B9B6N48LR/0zjnSpXiK4OlLKB39NutLGCP'
export GSHEETS_CREDENTIALS="$GB_REPO/ensembl-common/private/credentials.json"

alias dbhc="$GB_REPO/ensembl-common/scripts/dbhc.sh"
alias dbcopy="$GB_REPO/ensembl-common/scripts/dbcopy.sh"
alias run_testsuite="$GB_REPO/ensembl-common/scripts/run_testsuite.sh"

alias mkgbdir="mkdir -m 2775"

reload_ensembl_release() {
  EVERSION=`mysql-ens-meta-prod-1 ensembl_metadata -NB -e "SELECT ensembl_version FROM data_release WHERE is_current = 1 ORDER BY ensembl_version DESC LIMIT 1"`
  if [[ $EVERSION -gt ${ENSEMBL_RELEASE:-0} ]]; then
    export ENSEMBL_RELEASE=$EVERSION
  elif [[ $EVERSION -lt $ENSEMBL_RELEASE ]];then
    echo "Something is wrong: ENSEMBL_RELEASE=$ENSEMBL_RELEASE and ensembl_production_$EVERSION"
    return 1
  fi
}

reload_ensembl_release
