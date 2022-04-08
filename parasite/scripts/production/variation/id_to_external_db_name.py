#!/hps/software/users/wormbase/parasite/shared/.pyenv/versions/production-tools/bin/python

idregex_to_db = {
    '^WBVar\d{8}$'  : 'WORMBASE_VARIATION',
    '^WBRNAi\d{8}$' : 'WORMBASE_RNAI',
    '^RNAi$'        : 'RNAi_STUDY'
}