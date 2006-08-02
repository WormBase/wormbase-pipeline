#!/bin/csh
# sets the paths to ensembl conform directories
unalias cd
unalias ls
setenv CVSROOT :ext:wormpipe@cvs.sanger.ac.uk:/nfs/ensembl/cvsroot
setenv CVS_RSH ssh
setenv ZOE /usr/local/ensembl/Zoe
setenv PERL5LIB /nfs/acari/wormpipe/ensembl/ensembl/misc-scripts/xref_mapping/:/nfs/acari/wormpipe/ensembl/ensembl-config/celegans/WB158:/nfs/acari/wormpipe/ensembl/ensembl-pipeline/scripts:.:/nfs/acari/wormpipe/ensembl/ensembl/modules:/nfs/acari/wormpipe/ensembl/ensembl-pipeline/modules:/nfs/acari/wormpipe/ensembl/ensembl/modules:/nfs/acari/wormpipe/ensembl/bioperl-live
