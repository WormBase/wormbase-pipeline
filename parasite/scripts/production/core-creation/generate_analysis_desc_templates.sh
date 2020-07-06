## generate a template for analysis descriptions of new genomes, prior to import into the production db

cp ${PARASITE_ANALYSIS_DIR}/analysis_descriptions.wbps{$PREVIOUS_PARASITE_VERSION}.txt ${PARASITE_ANALYSIS_DIR}/analysis_descriptions.wbps{$PARASITE_VERSION}.template.txt

for genome in $(ls ${PARASITE_DATA}); do 
	echo ${genome}_gene$'\t' \
	'Gene models produced by <a href="?">?</a>, as described by <a href="?">?</a>'$'\t' \
	'Coding'$'\t' \
	'1'$'\t' \
	'1'$'\t' \
	'103'$'\t' \
	'NULL' >> ${PARASITE_ANALYSIS_DIR}/analysis_descriptions.wbps{$PARASITE_VERSION}.template.txt

	echo ${genome}_non_coding$'\t' \
	'Gene models produced by <a href="?">?</a>, as described by <a href="?">?</a>'$'\t' \
	'Non coding'$'\t' \
	'1'$'\t' \
	'1'$'\t' \
	'11'$'\t' \
	'NULL' >> ${PARASITE_ANALYSIS_DIR}/analysis_descriptions.wbps{$PARASITE_VERSION}.template.txt

done
