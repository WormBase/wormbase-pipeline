import os
from optparse import OptionParser
import pandas as pd
from ProductionMysql import *
from ProductionUtils import *
import gzip
import shutil

# Parse user input
parser = OptionParser(usage='usage: %prog [options] arguments')
parser.add_option("-o", "--output_dir", dest="OUTPUT_DIR",
                  help="Output directory for the phenotype homology dumps (Required)")
parser.add_option("-f", "--force", dest="force", action="store_true", default=False,
                  help="Force to override any pre-existing dumps (Optional) (Default: Not enabled)")
parser.add_option("-s", "--species", dest="SPECIES",
                  help="Only export homology phenotype dumps for the core dbs in the list, "
                       "e.g. acanthocheilonema_viteae_prjeb1697_core_17_105_1,acrobeloides_nanus_prjeb26554_core_17_105_1 (Optional)")
(options, args) = parser.parse_args()

OUTPUT_DIR = options.OUTPUT_DIR

# Check Input
if OUTPUT_DIR is None:
    print(dtnow() + ": ERROR - No output directory has been specified. Exiting.")
    raise ValueError
else:
    OUTPUT_DIR = OUTPUT_DIR.rstrip('/')

print(options.SPECIES)

if options.SPECIES is not None:
    SPECIES = [x.strip() for x in options.SPECIES.split(",") if x.strip() in staging.core_databases]
    for pre_spe in options.SPECIES.split(','):
        if pre_spe.strip() not in staging.core_databases:
            print(dtnow() + ": INFO - "+pre_spe+" is not a core db in the staging server. Skipping...")
    if len(SPECIES) == 0:
        print(dtnow() + ": ERROR - None of the specified species (--species) is a core db. Exiting.")
        raise ValueError
else:
    SPECIES = staging.core_dbs("core_"+os.environ["PARASITE_VERSION"]+"_"+os.environ["ENSEMBL_VERSION"])

VARIATION_DBS = staging.variation_dbs("variation_"+PARASITE_VERSION+"_"+ENSEMBL_VERSION)
COMPARATORS = {Variation(staging.host, x).species_name():Variation(staging.host, x).core_dbname() for x in VARIATION_DBS}

for cdb in SPECIES:
    ortho_file = OUTPUT_DIR + "/" + Core(staging.host, cdb).ftp_filename_n_filename() + ".orthologs.tsv.gz"
    output_file = OUTPUT_DIR + "/" + Core(staging.host, cdb).ftp_filename_n_filename() + ".orthology-inferred_phenotypes.gaf"

    if os.path.isfile(output_file):
        if options.force is False:
            print_info(output_file + " already exists and the --force option has not been enabled. Skipping.")
            continue


    taxon_id = Core(staging.host, cdb).meta_value("species.taxonomy_id").strip()
    if taxon_id == "":
        exit_with_error("Couldn't find a taxon_id for: "+cdb+". Exiting")
    if not os.path.isfile(ortho_file):
        print_info("There's no orthologs file for " + cdb + " ("+ortho_file+"). Skipping.")
        continue

    print_info("Dumping orthology-inferred phenotypes for: " + cdb)

    if os.path.isfile(output_file+".gz"):
        if options.force is False:
            print_info(output_file + " already exists and the --force option has not been enabled. Skipping.")
            continue

    print_info("Orthologs file detected: "+ortho_file)
    print_info("Taxon_id: "+taxon_id)

    print_info("Dumping to: " + output_file)
    # Open output file connection:
    if os.path.isfile(output_file+".gz"): os.remove(output_file+".gz")
    if os.path.isfile(output_file): os.remove(output_file)
    f = open(output_file, 'a')

    ## Write top-file comments:
    f.write("!gaf-version 2.1\n")
    f.write("!GO Annotation File (GAF) 2.1 Description: http://geneontology.org/docs/go-annotation-file-gaf-format-2.1\n")
    f.write("!generated-by: WormBase ParaSite\n")
    f.write("!project-URL: https://parasite.wormbase.org\n")
    f.write("!project-release: WBPS"+os.environ['PARASITE_VERSION']+"\n")
    f.write("!columns descriptions:\n")
    f.write("! 1) DB\n")
    f.write("! 2) DB Object ID (Gene ID)\n")
    f.write("! 3) DB Object Symbol (Gene Symbol)\n")
    f.write("! 4) Qualifier (NOT = It has been demonstrated that the gene is not associated with the phenotype)\n")
    f.write("! 5) Phenotype ID\n")
    f.write("! 6) DB:Reference\n")
    f.write("! 7) Evidence Code (http://geneontology.org/page/guide-go-evidence-codes): Inferred from Sequence "
            "Orthology (ISO)\n")
    f.write("! 8) Orthologous Species:Gene that the gene association has been inferred from\n")
    f.write("! 9) Aspect (P = biological process)\n")
    f.write("! 10) DB Object Type\n")
    f.write("! 11) Taxon\n")
    f.write("! 12) Release Date\n")
    f.write("! 13) Assigned By\n")
    f.write("! 14) Annotation Extension (Phenotype description provided by https://github.com/obophenotype/c-elegans-phenotype-ontology)\n")

    # Find comparators (mansoni and/or elegans)
    COMPARE_AGAINST = {x:COMPARATORS[x] for x in COMPARATORS if COMPARATORS[x]!=cdb}
    rows_counter = 0
    for comparator in COMPARE_AGAINST:

        comp_cdb = COMPARE_AGAINST[comparator]
        comparator_name = Core(staging.host, comp_cdb).meta_value("species.display_name").strip()

        print_info("Comparing against: " + comparator_name)

        pheno_file = OUTPUT_DIR + "/" + Core(staging.host, comp_cdb).ftp_filename_n_filename() + ".phenotypes.gaf.gz"
        pheno_df = pd.read_csv(pheno_file, comment="!", header=None, sep='\t')

        ortho_pre_df = pd.read_csv(ortho_file, header=0, sep='\t', iterator=True, chunksize=1000)
        ortho_df = pd.concat([chunk[chunk['ortholog_species_name'].str.startswith(comparator_name)] for chunk in ortho_pre_df])

        if len(set(ortho_df['ortholog_species_name'].to_list())) !=1:
            exit_with_error(ortho_file+" has multiple different ortholog_species_name values starting with "+comparator_name+". Exiting.")

        merged_df = pd.merge(ortho_df[["gene_id","ortholog_gene_id"]],pheno_df, left_on='ortholog_gene_id', right_on=1)
        merged_df["with_from"] = comparator_name.replace(' ','_') + ":" + merged_df["ortholog_gene_id"]
        merged_df["code"] = "ISO"
        merged_df["taxon_id"] = "taxon:"+taxon_id

        final_df = merged_df[[0,"gene_id","gene_id",3,4,5,"code","with_from",7,8,10,11,12]]

        ## Write final table
        rows_counter += len(final_df.index)
        final_df.to_csv(f, sep="\t", header=None, index=False)

    f.close()

    # Count lines
    print_info("Wrote "+str(rows_counter)+" lines to output.")
    count = 0
    with open(output_file, 'r') as fp:
        for line in fp:
            if line.startswith("!"): continue
            count += 1
            data = line.split("\t")
            if len(data)!=13:
                exit_with_error("This line has "+str(len(data))+" instead of 13 fields as expected. Exiting...\n"+line)
    if count < 1000:
        print_warning("The output file has less than 1000 lines. Is this ok?")

    print_info("Gzipping output")
    with open(output_file, 'rb') as f_in:
        with gzip.open(output_file+".gz", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    if os.path.isfile(output_file+".gz"):
        os.remove(output_file)

    print_info("Done\n")
