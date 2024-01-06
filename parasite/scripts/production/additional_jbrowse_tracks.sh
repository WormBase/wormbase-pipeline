#!/usr/bin/env bash 

JBROWSE_INSTALL_DIR="/hps/software/users/wormbase/parasite/software/jbrowse"
JBROWSE_OUT_DIR="$PARASITE_SCRATCH/jbrowse/WBPS${PARASITE_VERSION}/out"
TEMP_DIR="$PARASITE_SCRATCH/temp"
ISOSEQ_DIR="/nfs/production/flicek/wormbase/parasite/data/isoseq"
# on noah, this data is in /nfs/panda/ensemblgenomes/wormbase/parasite/projects/isoseq

# genomes in this list had their annotation updated in WBPS release 18
# so we add WBPS17 gene models as an additional track
genomes_updated_17=(					\
heterodera_glycines_prjna381081				\
necator_americanus_prjna72135 			\
strongyloides_stercoralis_prjeb528
)

mkdir -p $TEMP_DIR

# for genome in ${genomes_updated_17[@]}; do
# 	VERSION=17
# 	species=$(echo $genome | grep -o "^[^_]\+_[^_]\+")
# 	bioproject_lc=$(echo $genome | grep -o "[^_]\+$")
# 	bioproject_uc=${bioproject_lc^^}
# 	if [[ -f "$TEMP_DIR/${species}.${bioproject_uc}.WBPS${VERSION}.annotations.gff3" ]]; then
# 		echo "Found $TEMP_DIR/${species}.${bioproject_uc}.WBPS${VERSION}.annotations.gff3 - will not re-download"
# 	else
# 		gff=$(wget -q -O $TEMP_DIR/${species}.${bioproject_uc}.WBPS${VERSION}.annotations.gff3.gz ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS${VERSION}/species/${species}/${bioproject_uc}/${species}.${bioproject_uc}.WBPS${VERSION}.annotations.gff3.gz)

# 		echo "Unzipping $TEMP_DIR/${species}.${bioproject_uc}.WBPS${VERSION}.annotations.gff3.gz"

# 		gunzip $TEMP_DIR/${species}.${bioproject_uc}.WBPS${VERSION}.annotations.gff3.gz
# 	fi

# 	cmd="perl $JBROWSE_INSTALL_DIR/bin/flatfile-to-json.pl
# 	--gff $TEMP_DIR/${species}.${bioproject_uc}.WBPS${VERSION}.annotations.gff3
# 	--trackLabel WBPS${VERSION}_genemodels
# 	--key 'WBPS${VERSION} gene models'
# 	--out $JBROWSE_OUT_DIR/${genome}
# 	--trackType CanvasFeatures
# 	--type gene
# 	--config \"{ \\\"menuTemplate\\\" : [
# 			 { \\\"action\\\" : \\\"newWindow\\\" ,
# 			   \\\"label\\\" : \\\"View gene on WormBase ParaSite Release ${VERSION} archive site\\\" ,
# 			   \\\"url\\\" : \\\"https://release-${VERSION}.parasite.wormbase.org/Gene/Summary?g={name}\\\"}]}\"
# 	--metadata \"{ \\\"Track\\\" : \\\"WBPS${VERSION} gene models\\\",
#                        \\\"Category\\\" : \\\"Genome annotation\\\",
# 		       \\\"Description\\\" : \\\"Gene models from release ${VERSION} of WormBase ParaSite\\\" }\" "

# 	echo "Running $cmd"

# 	eval $cmd

# 	index_cmd="perl $JBROWSE_INSTALL_DIR/bin/generate-names.pl
# 		--out $JBROWSE_OUT_DIR/${genome}
# 		--tracks WBPS${VERSION}_genemodels
# 		--compress
# 		--incremental"

# 	echo "Running $index_cmd"

# 	eval $index_cmd
# done

# # genomes in this list had their annotation updated in WBPS release 19
# # so we add WBPS18 gene models as an additional track
# genomes_updated_18=(					\
# ancylostoma_ceylanicum_prjna231479				\
# )

# mkdir -p $TEMP_DIR

# for genome in ${genomes_updated_18[@]}; do
# 	VERSION=18
# 	species=$(echo $genome | grep -o "^[^_]\+_[^_]\+")
# 	bioproject_lc=$(echo $genome | grep -o "[^_]\+$")
# 	bioproject_uc=${bioproject_lc^^}
# 	if [[ -f "$TEMP_DIR/${species}.${bioproject_uc}.WBPS${VERSION}.annotations.gff3" ]]; then
# 		echo "Found $TEMP_DIR/${species}.${bioproject_uc}.WBPS${VERSION}.annotations.gff3 - will not re-download"
# 	else
# 		gff=$(wget -q -O $TEMP_DIR/${species}.${bioproject_uc}.WBPS${VERSION}.annotations.gff3.gz ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS${VERSION}/species/${species}/${bioproject_uc}/${species}.${bioproject_uc}.WBPS${VERSION}.annotations.gff3.gz)

# 		echo "Unzipping $TEMP_DIR/${species}.${bioproject_uc}.WBPS${VERSION}.annotations.gff3.gz"

# 		gunzip $TEMP_DIR/${species}.${bioproject_uc}.WBPS${VERSION}.annotations.gff3.gz
# 	fi

# 	cmd="perl $JBROWSE_INSTALL_DIR/bin/flatfile-to-json.pl
# 	--gff $TEMP_DIR/${species}.${bioproject_uc}.WBPS${VERSION}.annotations.gff3
# 	--trackLabel WBPS${VERSION}_genemodels
# 	--key 'WBPS${VERSION} gene models'
# 	--out $JBROWSE_OUT_DIR/${genome}
# 	--trackType CanvasFeatures
# 	--type gene
# 	--config \"{ \\\"menuTemplate\\\" : [
# 			 { \\\"action\\\" : \\\"newWindow\\\" ,
# 			   \\\"label\\\" : \\\"View gene on WormBase ParaSite Release ${VERSION} archive site\\\" ,
# 			   \\\"url\\\" : \\\"https://release-${VERSION}.parasite.wormbase.org/Gene/Summary?g={name}\\\"}]}\"
# 	--metadata \"{ \\\"Track\\\" : \\\"WBPS${VERSION} gene models\\\",
#                        \\\"Category\\\" : \\\"Genome annotation\\\",
# 		       \\\"Description\\\" : \\\"Gene models from release ${VERSION} of WormBase ParaSite\\\" }\" "

# 	echo "Running $cmd"

# 	eval $cmd

# 	index_cmd="perl $JBROWSE_INSTALL_DIR/bin/generate-names.pl
# 		--out $JBROWSE_OUT_DIR/${genome}
# 		--tracks WBPS${VERSION}_genemodels
# 		--compress
# 		--incremental"

# 	echo "Running $index_cmd"

# 	eval $index_cmd
# done
# ### Other special tracks
# # by special request we also added an additional track for alternative N americanus gene models
# # The GFF lives on the FTP site

# if [[ -f "$TEMP_DIR/NAME_17JAN2018.genome.gff" ]]; then
# 	echo "Found $TEMP_DIR/NAME_17JAN2018.genome.gff - will not redownload"
# else
# 	gff=$(wget -q -O $TEMP_DIR/NAME_17JAN2018.genome.gff.gz ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/datasets/logan_2020_32453752/NAME_17JAN2018.genome.gff.gz)

# 	echo "Unzipping $TEMP_DIR/NAME_17JAN2018.genome.gff.gz"

# 	gunzip $TEMP_DIR/NAME_17JAN2018.genome.gff.gz
# fi

# cmd="perl $JBROWSE_INSTALL_DIR/bin/flatfile-to-json.pl
#     --gff $TEMP_DIR/NAME_17JAN2018.genome.gff
#         --trackLabel alt_gene_models_logan
#         --key 'Alternative gene models (Logan et al., 2020)'
#         --out $JBROWSE_OUT_DIR/necator_americanus_prjna72135
#         --trackType CanvasFeatures
#         --type gene
#         --metadata \"{ \\\"Track\\\" : \\\"Alternative gene models (Logan et al., 2020)\\\",
# 		       \\\"Category\\\" : \\\"Genome annotation\\\",
#                        \\\"Study\\\" : \\\"<a href=\"https://pubmed.ncbi.nlm.nih.gov/32453752/\">Comprehensive analysis of the secreted proteome of adult Necator americanus hookworms (Logan et al., 2020)</a>\\\"}\" "

# echo "Running $cmd"

# eval $cmd

# index_cmd="perl $JBROWSE_INSTALL_DIR/bin/generate-names.pl
# 	--out $JBROWSE_OUT_DIR/necator_americanus_prjna72135
# 	--tracks alt_gene_models_logan
# 	--compress
# 	--incremental"

# echo "Running $index_cmd"

# eval $index_cmd

# Another alternative set for necator

# if [[ -f "$TEMP_DIR/WashU_ES_genome_2023.01.12.gff" ]]; then
# 	echo "Found $TEMP_DIR/WashU_ES_genome_2023.01.12.gff - will not redownload"
# else
# 	gff=$(wget -q -O $TEMP_DIR/WashU_ES_genome_2023.01.12.gff https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/datasets/PRJNA72135_alt/WashU_ES_genome_2023.01.12.gff)

# fi

# cmd="perl $JBROWSE_INSTALL_DIR/bin/flatfile-to-json.pl
#     --gff $TEMP_DIR/WashU_ES_genome_2023.01.12.gff
#         --trackLabel alt_gene_models_uzoechi
#         --key 'Alternative gene models (Uzoechi et al., 2023)'
#         --out $JBROWSE_OUT_DIR/necator_americanus_prjna72135
#         --trackType CanvasFeatures
#         --type gene
#         --metadata \"{ \\\"Track\\\" : \\\"Alternative gene models (Uzoechi et al., 2023)\\\",
# 		       \\\"Category\\\" : \\\"Genome annotation\\\",
#                        \\\"Study\\\" : \\\"<a href=\"https://pubmed.ncbi.nlm.nih.gov/36678443/\">Excretory/Secretory Proteome of Females and Males of the Hookworm Ancylostoma ceylanicum (Uzoechi et al., 2023)</a>\\\"}\" "

# echo "Running $cmd"

# eval $cmd

# index_cmd="perl $JBROWSE_INSTALL_DIR/bin/generate-names.pl
# 	--out $JBROWSE_OUT_DIR/necator_americanus_prjna72135
# 	--tracks alt_gene_models_uzoechi
# 	--compress
# 	--incremental"

# echo "Running $index_cmd"

# eval $index_cmd

# # IsoSeq for D. immitis

# mkdir -p $JBROWSE_OUT_DIR/dirofilaria_immitis_prjeb1797/bam

# cp $ISOSEQ_DIR/SRR12046914.sort.bam* $JBROWSE_OUT_DIR/dirofilaria_immitis_prjeb1797/bam
# cp $ISOSEQ_DIR/SRR12046915.sort.bam* $JBROWSE_OUT_DIR/dirofilaria_immitis_prjeb1797/bam

# bam_cmd="perl $JBROWSE_INSTALL_DIR/bin/add-bam-track.pl
# 	--in $JBROWSE_OUT_DIR/dirofilaria_immitis_prjeb1797/trackList.json
# 	--out $JBROWSE_OUT_DIR/dirofilaria_immitis_prjeb1797/trackList.json
# 	--bam_url bam/SRR12046914.sort.bam
# 	--label 'SRR12046914: Adult female'
# 	--config \"{ \\\"metadata\\\" : {
# 			 \\\"Track\\\" : \\\"SRR12046915: Adult female\\\",
# 			 \\\"Category\\\" : \\\"IsoSeq\\\",
# 			 \\\"Developmental stage\\\" : \\\"adult\\\",
# 			 \\\"ENA BioProject\\\" : \\\"<a href=\"https://www.ebi.ac.uk/ena/browser/view/PRJNA640410\">Study page: PRJNA640410</a>\\\",
# 			 \\\"ENA study\\\" : \\\"<a href=\"https://www.ebi.ac.uk/ena/data/view/SRP267883\">Study page: SRP267883</a>\\\",
# 			 \\\"Sex\\\" : \\\"female\\\",
# 			 \\\"Submitting centre\\\" : \\\"University of Wisconsin-Madison\\\" }} \" "

# echo "Running $bam_cmd"

# eval $bam_cmd

# bam_cmd="perl $JBROWSE_INSTALL_DIR/bin/add-bam-track.pl
#         --in $JBROWSE_OUT_DIR/dirofilaria_immitis_prjeb1797/trackList.json
#         --out $JBROWSE_OUT_DIR/dirofilaria_immitis_prjeb1797/trackList.json
#         --bam_url bam/SRR12046915.sort.bam
#         --label 'SRR12046915: Adult male'
#         --config \"{ \\\"metadata\\\" : {
#                          \\\"Track\\\" : \\\"SRR12046915: Adult male\\\",
#                          \\\"Category\\\" : \\\"IsoSeq\\\",
#                          \\\"Developmental stage\\\" : \\\"adult\\\",
#                          \\\"ENA BioProject\\\" : \\\"<a href=\"https://www.ebi.ac.uk/ena/browser/view/PRJNA640410\">Study page: PRJNA640410</a>\\\",
#                          \\\"ENA study\\\" : \\\"<a href=\"https://www.ebi.ac.uk/ena/data/view/SRP267883\">Study page: SRP267883</a>\\\",
#                          \\\"Sex\\\" : \\\"male\\\",
#                          \\\"Submitting centre\\\" : \\\"University of Wisconsin-Madison\\\" }} \" "

# echo "Running $bam_cmd"

# eval $bam_cmd

# Low confidence gene models for Schmidtea 1
if [[ -f "$TEMP_DIR/schmidtea_mediterranea_s2f19h1.low_confidence.gff3" ]]; then
	echo "Found $TEMP_DIR/schmidtea_mediterranea_s2f19h1.low_confidence.gff3 - will not redownload"
else
	gff=$(wget -q -O $TEMP_DIR/schmidtea_mediterranea_s2f19h1.low_confidence.gff3 https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/datasets/S2F19H1_PRJNA885486/schmidtea_mediterranea_s2f19h1.low_confidence.gff3)
fi

cmd="perl $JBROWSE_INSTALL_DIR/bin/flatfile-to-json.pl
    --gff $TEMP_DIR/schmidtea_mediterranea_s2f19h1.low_confidence.gff3
        --trackLabel low_confidence_genes
        --key 'Low confidence gene models'
        --out $JBROWSE_OUT_DIR/schmidtea_mediterranea_s2f19h1prjna885486
        --trackType CanvasFeatures
        --type gene
        --metadata \"{ \\\"Track\\\" : \\\"Low confidence gene models\\\",
		       \\\"Category\\\" : \\\"Genome annotation\\\"}\" "

echo "Running $cmd"

eval $cmd

index_cmd="perl $JBROWSE_INSTALL_DIR/bin/generate-names.pl
	--out $JBROWSE_OUT_DIR/schmidtea_mediterranea_s2f19h1prjna885486
	--tracks low_confidence_genes
	--compress
	--incremental"

echo "Running $index_cmd"

eval $index_cmd

# Low confidence gene models for Schmidtea 2
if [[ -f "$TEMP_DIR/schmidtea_mediterranea_s2f19h2.low_confidence.gff3" ]]; then
	echo "Found $TEMP_DIR/schmidtea_mediterranea_s2f19h2.low_confidence.gff3 - will not redownload"
else
	gff=$(wget -q -O $TEMP_DIR/schmidtea_mediterranea_s2f19h2.low_confidence.gff3 https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/datasets/S2F19H1_PRJNA885486/schmidtea_mediterranea_s2f19h2.low_confidence.gff3)
fi

cmd="perl $JBROWSE_INSTALL_DIR/bin/flatfile-to-json.pl
    --gff $TEMP_DIR/schmidtea_mediterranea_s2f19h2.low_confidence.gff3
        --trackLabel low_confidence_genes
        --key 'Low confidence gene models'
        --out $JBROWSE_OUT_DIR/schmidtea_mediterranea_s2f19h2prjna885486
        --trackType CanvasFeatures
        --type gene
        --metadata \"{ \\\"Track\\\" : \\\"Low confidence gene models\\\",
		       \\\"Category\\\" : \\\"Genome annotation\\\"}\" "

echo "Running $cmd"

eval $cmd

index_cmd="perl $JBROWSE_INSTALL_DIR/bin/generate-names.pl
	--out $JBROWSE_OUT_DIR/schmidtea_mediterranea_s2f19h2prjna885486
	--tracks low_confidence_genes
	--compress
	--incremental"

echo "Running $index_cmd"

eval $index_cmd