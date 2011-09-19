# set up Ensembl database environment

# this variable name is required by Compara
setenv ENSEMBL_CVS_ROOT_DIR "/software/worm/ensembl"

setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/bioperl-live:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-compara/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-variation/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-functgenomics/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-pipeline/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-pipeline/scripts:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-analysis/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl/modules:${PERL5LIB}

# used when setting up a new Ensembl Genomes db
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-config/celegans/WS220:${PERL5LIB}

# for the Compara analysis pipeline
set path = (${ENSEMBL_CVS_ROOT_DIR}/ensembl-hive/scripts /software/worm/bin/wublast/ /software/worm/bin /software/bin $path)
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-config/celegans/generic-elegans:${PERL5LIB}
setenv PERL5LIB /software/worm/ensembl/bioperl-run:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-hive/modules:${PERL5LIB}
setenv PERL5LIB /software/worm/lib:${PERL5LIB}
# this is used in the Compara Build docs section on troubleshooting
setenv COMPARA_URL "mysql://wormadmin:dbworms@farmdb1:3306/worm_compara_homology_63WS${WORMBASE_RELEASE}"


# for the EFuncGen database
source ~/wormbase/scripts/ENSEMBL/etc/efg.config

# stick this stuff at the end
setenv PERL5LIB ${PERL5LIB}:.
setenv PERL5LIB ${PERL5LIB}:/software/worm/lib/site_perl
setenv PERL5LIB ${PERL5LIB}:/nfs/wormpub/wormbase/scripts/
setenv PERL5LIB ${PERL5LIB}:~/wormbase/scripts/ENSEMBL/lib

