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

###############
# Convenience environment for Build
###############

setenv CVS_DIR $WORMPUB/wormbase-pipeline/scripts 
setenv SCRIPTS $CVS_DIR
setenv BUILD_HOME $WORMPUB
setenv ACEDB_NO_BANNER





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


# To do:
# PATH
set worm_packages = "/net/isilon3/production/panda/ensemblgenomes/wormbase/software/packages"
set path  = ($worm_packages/bamtools/bin $worm_packages/bowtie $worm_packages/cufflinks $worm_packages/samtools $worm_packages/tophat $path)

# PERL5LIB
# BLASTMAT etc
