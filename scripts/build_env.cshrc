###############
# Convenience environment for Build
###############

# indicate that we are on the Sanger system, not the EBI
setenv SANGER 1

setenv WORMPUB ~wormpub
setenv WORMBASE $WORMPUB # wormpub and wormbase are interchangeable terms
setenv WORM_SW_ROOT /software/worm
setenv PERLBREW_ROOT ${WORMPUB}/software/packages/perlbrew
setenv WORM_PACKAGES $WORM_SW_ROOT
setenv WORM_BIN ${WORM_SW_ROOT}/bin
setenv WORM_LIB ${WORM_SW_ROOT}/lib
setenv CVS_DIR $WORMPUB/wormbase-pipeline/scripts 
setenv SCRIPTS $CVS_DIR
setenv BUILD_HOME $WORMPUB
setenv PIPELINE /lustre/scratch109/ensembl/wormpipe

#setenv LSFPATHS $PIPELINE # Resources hint for LSF
setenv LSB_DEFAULTQUEUE normal
setenv ACEDB_NO_BANNER

###############
# Ensembl database details
###############
setenv WORM_DBHOST farmdb1
setenv WORM_DBPORT 3306
