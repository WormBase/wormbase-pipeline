# set up Ensembl database environment

if ($?prompt) then
    if ( $SHELL:t == "lstcsh" || $SHELL:t == "tcsh" ) then
        set prompt='%m[%$USER Build %$SPECIES]\!: '
    else
        set prompt="`hostname -s`[Build]\!: "
    endif
endif

#############################################
# are we are on the Sanger system, or the EBI
#############################################
if (-e "/software/worm") then
    setenv SANGER 1                

    setenv WORMPUB ~wormpub

    setenv WORM_SW_ROOT /software/worm
    
    setenv WORM_PACKAGES $WORM_SW_ROOT
    
    setenv PIPELINE /lustre/scratch109/ensembl/wormpipe
    
    setenv LSB_DEFAULTQUEUE normal
    
    # Ensembl MySQL database details
    setenv WORM_DBHOST farmdb1
    setenv WORM_DBPORT 3306

    # the filesystem's group name for WormBase
    setenv WORM_GROUP_NAME worm

    # LSF submit multiple for -M
    setenv LSF_SUBMIT_MULTIPLE "000"

    # WU-BLAST
    setenv BLASTDB /data/blastdb/Ensembl/
    setenv BLASTMAT /usr/local/ensembl/data/blastmat/
    setenv BLASTFILTER /usr/local/ensembl/bin
    setenv WU_BLAST_PATH /usr/local/ensembl/bin/

    setenv COILSDIR /usr/local/ensembl/data/coils
    setenv LSB_DEFAULTPROJECT wormbase

    # NCBI-BLAST
    setenv NCBI_BLAST_PATH /software/bin/

    # PATH
    setenv EG_BIN /usr/local/ensembl/bin
    setenv WORM_BIN ${WORM_SW_ROOT}/bin
    set path  = (${EG_BIN} $path)
    # for the Compara analysis pipeline
    setenv ENSEMBL_CVS_ROOT_DIR "${WORM_PACKAGES}/ensembl"
    set path = (${ENSEMBL_CVS_ROOT_DIR}/eHive/scripts ${WORM_BIN}/wublast/ ${WORM_BIN} /software/bin ${path})

    # PERL5LIB
    setenv PERL5LIB ${WORM_SW_ROOT}/lib/site_perl
    setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/bioperl-1.2.3:${PERL5LIB}
    setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/bioperl-run:${PERL5LIB}


else
##########################
# we are on the EBI system
##########################


    setenv WORMPUB /nfs/panda/ensemblgenomes/wormbase

    setenv WORM_SW_ROOT ${WORMPUB}/software
    
    setenv WORM_PACKAGES ${WORM_SW_ROOT}/packages
    
    setenv EG_PACKAGES /nfs/panda/ensemblgenomes/external/

    setenv CVS_DIR $WORMPUB/wormbase-pipeline/scripts 

    setenv PIPELINE /gpfs/nobackup/ensembl_genomes/wormbase/BUILD/pipeline
    
    setenv LSB_DEFAULTQUEUE production-rh6

    # Ensembl MySQL database details
    setenv WORM_DBHOST mysql-wormbase-pipelines
    setenv WORM_DBPORT 4331

    # the filesystem's group name for WormBase
    setenv WORM_GROUP_NAME nucleotide

    # LSF submit multiple for -M  (the -M is the same as the rusage value on the EBI)
    setenv LSF_SUBMIT_MULTIPLE ""

    # WU-BLAST
    setenv BLASTDB ${PIPELINE}/blastdb/Ensembl/
    setenv BLASTMAT ${EG_PACKAGES}/wublast/matrix
    setenv BLASTFILTER ${EG_PACKAGES}/wublast/filter
    setenv WU_BLAST_PATH ${EG_PACKAGES}/wublast/

    # NCBI-BLAST
    setenv NCBI_BLAST_PATH ${WORM_PACKAGES}/ncbi-blast/bin

    # PATH
    setenv EG_BIN ${EG_PACKAGES}/bin
    setenv WORM_BIN ${WORM_SW_ROOT}/bin
    set path  = (${EG_BIN} ${WORM_BIN} $path)
    # for the Compara analysis pipeline
    setenv ENSEMBL_CVS_ROOT_DIR "${WORM_PACKAGES}/ensembl"
    set path = (${ENSEMBL_CVS_ROOT_DIR}/ensembl-hive/scripts ${WORM_PACKAGES}/wublast/ ${WORM_BIN} ${path})
    # to use the ensembl-geneomes version of perl which knows about Oracle (but not ACE.pm)
    set path = (/nfs/panda/ensemblgenomes/perl/perlbrew/perls/5.14.2/bin/ ${path})

    # use java version 1.6
    setenv JAVA_HOME /sw/arch/pkg/jdk1.6
    setenv ANT_HOME /nfs/panda/ensemblgenomes/external/apache-ant
    setenv PATH ${JAVA_HOME}/bin:${ANT_HOME}/bin:${PATH}

    # ruby and rake
    set path = ( /net/isilon3/production/panda/ensemblgenomes/wormbase/software/packages/gem/bin/ $path )
    setenv RUBYLIB /net/isilon3/production/panda/ensemblgenomes/wormbase/software/packages/gem/lib/
    setenv GEM_HOME /net/isilon3/production/panda/ensemblgenomes/wormbase/software/packages/gem

    # PERL5LIB
    setenv PERL5LIB ${WORM_PACKAGES}/bioperl/bioperl-live
    setenv PERL5LIB ${WORM_PACKAGES}/bioperl/bioperl-run/lib:${PERL5LIB}


endif

############################################
# these are now system-independent locations
############################################

setenv LSFPATHS $PIPELINE # Resources hint for LSF

setenv WORMBASE $WORMPUB # wormpub and wormbase are interchangeable terms
setenv WORM_BIN ${WORM_SW_ROOT}/bin
setenv WORM_LIB ${WORM_SW_ROOT}/lib

####################
# perl 5 libraries
####################

setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-compara/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-variation/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-funcgen/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-pipeline/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-pipeline/scripts:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-analysis/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-production/modules:${PERL5LIB}
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl/modules:${PERL5LIB}

# used when setting up a new Ensembl Genomes db
setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-config/nematodes/generic_nematode:${PERL5LIB}

setenv PERL5LIB ${ENSEMBL_CVS_ROOT_DIR}/ensembl-hive/modules:${PERL5LIB}
setenv PERL5LIB ${WORM_SW_ROOT}/lib:${PERL5LIB}

setenv WORMBASE_RELEASE `grep "NAME WS" ${WORMPUB}/BUILD/autoace/wspec/database.wrm | perl -ne 's/NAME WS//;print $_'`
setenv COMPARA_URL "mysql://wormadmin:worms@${WORM_DBHOST}:${WORM_DBPORT}/worm_compara_homology_WS${WORMBASE_RELEASE}"

setenv ENSEMBL_REGISTRY $CVS_DIR/ENSEMBL/etc/Registry.pm

# for the EFuncGen database
#source ${CVS_DIR}/ENSEMBL/etc/efg.config

# stick this stuff at the end
setenv PERL5LIB ${PERL5LIB}:.
setenv PERL5LIB ${PERL5LIB}:${CVS_DIR}
setenv PERL5LIB ${PERL5LIB}:${CVS_DIR}/ENSEMBL/lib

