import pandas as pd

comparator = "Schistosoma mansoni (PRJEA36577)"

pheno_file="/homes/digri/schistosoma_mansoni.PRJEA36577.WBPS17.phenotypes.wb.gz"
ortho_file="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/acanthocheilonema_viteae/PRJEB1697/acanthocheilonema_viteae.PRJEB1697.WBPS16.orthologs.tsv.gz"

pheno_df = pd.read_csv(pheno_file, comment="!", header=None, sep='\t')
ortho_pre_df = pd.read_csv(ortho_file, header=0, sep='\t', iterator=True, chunksize=1000)
ortho_df = pd.concat([chunk[chunk['ortholog_species_name'] == comparator] for chunk in ortho_pre_df])

