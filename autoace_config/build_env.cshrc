####################################################################################
#
# This is the WormBase common .cshrc that should be sourced by everyone in WormBase
#
####################################################################################

umask 002

unset autologout
set nobeep
xset b 0 0 0 
set noclobber

# set an environment variable to point to this directory for easy access
setenv WORMPUB /nfs/panda/ensemblgenomes/wormbase
setenv WORMBASE $WORMPUB # wormpub and wormbase are interchangeable terms
setenv WORM_SW_ROOT ${WORMPUB}/software
setenv PERLBREW_ROOT ${WORMPUB}/software/packages/perlbrew

###############
# Convenience environment for Build
###############

setenv CVS_DIR $WORMPUB/wormbase-pipeline/scripts 
setenv SCRIPTS $CVS_DIR
setenv BUILD_HOME $WORMPUB

setenv WORM_PACKAGES ${WORM_SW_ROOT}/packages
setenv WORM_BIN ${WORM_SW_ROOT}/bin
setenv PIPELINE /nfs/nobackup/ensemblgenomes/wormbase/BUILD/pipeline
setenv ACEDB_NO_BANNER

setenv LSB_DEFAULTQUEUE production-rh7

####################
# perl 5 libraries
####################
setenv PERL5LIB
setenv PERL5LIB ${PERL5LIB}:${WORM_SW_ROOT}/packages/bioperl/bioperl-1.2.3
setenv PERL5LIB ${PERL5LIB}:${WORM_SW_ROOT}/packages/bioperl/bioperl-run
setenv PERL5LIB ${PERL5LIB}:${WORM_SW_ROOT}/packages/bioperl/bioperl-run/lib

# some build scripts connect to the ensembl databases in a limited way, but
# do not need the full complement of modules for these
setenv PERL5LIB ${PERL5LIB}:${WORM_SW_ROOT}/packages/ensembl/ensembl/modules
setenv PERL5LIB ${PERL5LIB}:${WORM_SW_ROOT}/packages/ensembl/ensembl-compara/modules
setenv PERL5LIB ${PERL5LIB}:${SCRIPTS}

setenv WORM_SW_ROOT ${WORMPUB}/software
setenv WORM_BIN ${WORM_SW_ROOT}/bin
setenv EG_PACKAGES /nfs/panda/ensemblgenomes/external
setenv EG_BIN ${EG_PACKAGES}/bin

# PATH
set path  = (${EG_BIN} ${EG_BIN}/blat ${EG_BIN}/exonerate-2/bin ${WORM_BIN} ${WORM_PACKAGES}/bamtools/bin ${WORM_PACKAGES}/bowtie ${WORM_PACKAGES}/cufflinks ${WORM_PACKAGES}/samtools ${WORM_PACKAGES}/tophat ${EG_BIN}/EMBOSS/bin $path)

