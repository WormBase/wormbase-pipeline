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

###############
# Convenience environment for Build
###############

setenv CVS_DIR $WORMPUB/wormbase-pipeline/scripts 
setenv SCRIPTS $CVS_DIR
setenv BUILD_HOME $WORMPUB
setenv ACEDB_NO_BANNER

# To do:
# PATH
# PERL5LIB
# BLASTMAT etc
