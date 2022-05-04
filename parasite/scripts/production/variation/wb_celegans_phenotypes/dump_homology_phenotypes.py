import os
from optparse import OptionParser
import pandas as pd
from ProductionMysql import *
from ProductionUtils import *

# Parse user input
parser = OptionParser(usage='usage: %prog [options] arguments')
parser.add_option("-o", "--output_dir", dest="OUTPUT_DIR",
                  help="Output directory for the phenotype homology dumps (Required)")
parser.add_option("-f", "--force", dest="force", action="store_true", default=False,
                  help="Force to override any pre-existing dumps (Optional) (Default: Not enabled)")
parser.add_option("-s", "--species", dest="SPECIES",
                  help="Only export homology pheontype dumps for the core dbs in the list, "
                       "e.g. acanthocheilonema_viteae_prjeb1697_core_17_105_1,acrobeloides_nanus_prjeb26554_core_17_105_1 (Optional)")
(options, args) = parser.parse_args()

OUTPUT_DIR = options.OUTPUT_DIR

# Check Input
if OUTPUT_DIR is None:
    print(dtnow() + ": ERROR - No output directory has been specified. Exiting.")
    raise ValueError
if options.SPECIES is not None:
    SPECIES = [x.strip() for x in options.SPECIES.split(",") if x.strip() in staging.core_databases]
    for pre_spe in options.SPECIES.split(','):
        if pre_spe.strip() not in staging.core_databases:
            print(dtnow() + ": INFO - "+pre_spe+" is not a core db in the staging server. Skipping...")
    if len(SPECIES) == 0:
        print(dtnow() + ": ERROR - None of the specified species (--species) is a core db. Exiting.")
        raise ValueError
else:
    SPECIES = None

VARIATION_DBS = staging.variation_dbs("variation_"+PARASITE_VERSION+"_"+ENSEMBL_VERSION)
COMPARATORS = {Variation(staging.host, x).species():Variation(staging.host, x).core() for x in VARIATION_DBS}

for cdb in staging.core_databases:
    print_info("Working on "+cdb)
    for comparator in COMPARATORS:
        comp_cdb = COMPARATORS[comparator]
        if cdb == comp_cdb: continue
        pheno_file = OUTPUT_DIR+"/"+Core(staging.host, comp_cdb).ftp_filename_n_filename()
        if 



# pheno_file="/homes/digri/schistosoma_mansoni.PRJEA36577.WBPS17.phenotypes.wb.gz"
# ortho_file="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/acanthocheilonema_viteae/PRJEB1697/acanthocheilonema_viteae.PRJEB1697.WBPS16.orthologs.tsv.gz"
#
# pheno_df = pd.read_csv(pheno_file, comment="!", header=None, sep='\t')
# ortho_pre_df = pd.read_csv(ortho_file, header=0, sep='\t', iterator=True, chunksize=1000)
# ortho_df = pd.concat([chunk[chunk['ortholog_species_name'] == comparator] for chunk in ortho_pre_df])
#
# merged_df = pd.merge(ortho_df[["gene_id","ortholog_gene_id"]],pheno_df, left_on='ortholog_gene_id', right_on=1)
