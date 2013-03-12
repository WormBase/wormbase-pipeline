# set up Ensembl database environment

# this variable name is required by Compara
setenv ENSEMBL_CVS_ROOT_DIR "${WORM_PACKAGES}/ensembl"
setenv COMPARA "${WORM_PACKAGES}/ensembl"

setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/bioperl-live
setenv PERL5LIB ${COMPARA}/ensembl-compara/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-variation/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-functgenomics/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-pipeline/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-pipeline/scripts:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-analysis/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl/modules:${PERL5LIB}

# used when setting up a new Ensembl Genomes db
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-config/nematodes/celegans/generic-elegans:${PERL5LIB}

# for the Compara analysis pipeline
if ($SANGER) then
    set path = (${ENSEMBL_CVS_ROOT_DIR}/eHive/scripts ${WORM_BIN}/wublast/ ${WORM_BIN} /software/bin ${path})
else
    set path = (${ENSEMBL_CVS_ROOT_DIR}/ensembl-hive/scripts ${WORM_PACKAGES}/wublast/ ${WORM_BIN} ${path})
endif

setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-config/celegans/generic-elegans:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/bioperl-run:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-hive/modules:${PERL5LIB}
setenv PERL5LIB ${WORM_SW_ROOT}/lib:${PERL5LIB}
# this is used in the Compara Build docs section on troubleshooting
setenv COMPARA_URL "mysql://wormadmin:worms@${WORM_DBHOST}:${WORM_DBPORT}/worm_compara_homology_63WS${WORMBASE_RELEASE}"


# for the EFuncGen database
source ${CVS_DIR}/ENSEMBL/etc/efg.config

# stick this stuff at the end
setenv PERL5LIB ${PERL5LIB}:.
setenv PERL5LIB ${PERL5LIB}:${WORM_SW_ROOT}/lib/site_perl
setenv PERL5LIB ${PERL5LIB}:${CVS_DIR}
setenv PERL5LIB ${PERL5LIB}:${CVS_DIR}/ENSEMBL/lib

