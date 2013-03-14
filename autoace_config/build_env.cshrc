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
setenv PIPELINE /nfs/nobackup2/ensemblgenomes/wormbase/BUILD/pipeline
setenv WORM_PACKAGES ${WORM_SW_ROOT}/packages
setenv WORM_BIN ${WORM_SW_ROOT}/bin
setenv EG_PACKAGES /nfs/panda/ensemblgenomes/external
setenv EG_BIN ${EG_PACKAGES}/bin
setenv LSFPATHS $PIPELINE # Resources hint for LSF
setenv ACEDB_NO_BANNER

setenv LSB_DEFAULTQUEUE production-rh6

###############
# Ensembl MySQL database details
###############
setenv WORM_DBHOST mysql-wormbase-pipelines
setenv WORM_DBPORT 4331

####################
# perl 5 libraries
####################
setenv PERL5LIB ${WORM_SW_ROOT}/lib/perl5/site_perl
setenv PERL5LIB ${PERL5LIB}:${WORM_SW_ROOT}/packages/bioperl/bioperl-live

# for worm ensembl
setenv PERL5LIB ${PERL5LIB}:${WORM_SW_ROOT}/packages/ensembl/ensembl/modules
setenv PERL5LIB ${PERL5LIB}:${WORM_SW_ROOT}/packages/ensembl/ensembl-compara/modules
setenv PERL5LIB ${PERL5LIB}:${WORM_SW_ROOT}/packages/ensembl/ensembl-pipeline/modules
setenv PERL5LIB ${PERL5LIB}:${WORM_SW_ROOT}/packages/ensembl/ensembl-analysis/modules
setenv PERL5LIB ${PERL5LIB}:${SCRIPTS}

setenv ENSEMBL_REGISTRY $SCRIPTS/ENSEMBL/etc/E_registry.pl


# PATH
set path  = (${EG_BIN} ${WORM_BIN} ${WORM_PACKAGES}/bamtools/bin ${WORM_PACKAGES}/bowtie ${WORM_PACKAGES}/cufflinks ${WORM_PACKAGES}/samtools ${WORM_PACKAGES}/tophat $path)

# WU-BLAST
setenv BLASTDB ${PIPELINE}/blastdb/Ensembl/
setenv BLASTMAT ${EG_PACKAGES}/wublast/matrix
setenv BLASTFILTER ${EG_PACKAGES}/wublast/filter

