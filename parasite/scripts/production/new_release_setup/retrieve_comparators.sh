#!/usr/bin/env bash 

# find comparator dbs on E! and EG FTP sites and import to $PARASITE_STAGING_MYSQL.
# do not delete or overwrite if a database already exists for a species.

comparators=( \
amphimedon_queenslandica	\
capitella_teleta		\
ciona_intestinalis		\
crassostrea_gigas		\
danio_rerio			\
drosophila_melanogaster		\	
homo_sapiens			\	
ixodes_scapularis		\	
mus_musculus			\
nematostella_vectensis		\
rattus_norvegicus		\
saccharomyces_cerevisiae	\
schizosaccharomyces_pombe	\
trichoplax_adhaerens	)

mkdir -p ${PARASITE_SCRATCH2}/comparators_TMP

for comparator in ${comparators[@]}; do 
	mkdir -p ${PARASITE_SCRATCH2}/comparators_TMP/${comparator}
	
	# try ensembl first	

	ftp=$(rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-${ENSEMBL_VERSION}/mysql/${comparator}_core_${ENSEMBL_VERSION}_*/ ${PARASITE_SCRATCH2}/comparators_TMP/${comparator}/ )
	if [ $? -eq 0 ]; then 
		echo "Found ${comparator}_core_${ENSEMBL_VERSION}"
		continue	
	fi
	# then ensembl metazoa

	ftp=$(rsync -av rsync://ftp.ensemblgenomes.org/all/pub/metazoa/release-${EG_VERSION}/mysql/${comparator}_core_${EG_VERSION}_${ENSEMBL_VERSION}_*/ ${PARASITE_SCRATCH2}/comparators_TMP/${comparator}/ )

        if [ $? -eq 0 ]; then
                echo "Found ${comparator}_core_${EG_VERSION}_${ENSEMBL_VERSION}"
                continue
	fi
	# finally ensembl fungi

        ftp=$(rsync -av rsync://ftp.ensemblgenomes.org/all/pub/fungi/release-${EG_VERSION}/mysql/${comparator}_core_${EG_VERSION}_${ENSEMBL_VERSION}_*/ ${PARASITE_SCRATCH2}/comparators_TMP/${comparator}/ )

        if [ $? -eq 0 ]; then
                echo "Found ${comparator}_core_${EG_VERSION}_${ENSEMBL_VERSION}"
                continue
        fi

	echo "Could not find ${comparator}!"
	exit
done

for comparator in ${comparators[@]}; do
	
	sql=$(ls ${PARASITE_SCRATCH2}/comparators_TMP/${comparator}/*.sql.gz)
	regex="([^\/]*)\.sql.gz"
	if [[ $sql =~ $regex ]]; then
		db=${BASH_REMATCH[1]}
	else 
		echo "could not infer database name for ${comparator}"
		exit
	fi

	existing=$($PARASITE_STAGING_MYSQL -Ne "SHOW DATABASES LIKE \"%${db}%\" ")	

	if [[ "$existing" ]]; then

		echo "$existing already exists on $PARASITE_STAGING_MYSQL - will not update"
		continue
	fi

	echo "Unzipping ${comparator} files.."

	gunzip ${PARASITE_SCRATCH2}/comparators_TMP/${comparator}/*.gz

	$PARASITE_STAGING_MYSQL-ensrw -Ne "CREATE DATABASE $db"

	$PARASITE_STAGING_MYSQL-ensrw $db <  ${PARASITE_SCRATCH2}/comparators_TMP/${comparator}/${db}.sql

	echo "Populating ${db}.."

	$PARASITE_STAGING_MYSQL-ensrw mysqlimport -L ${db} ${PARASITE_SCRATCH2}/comparators_TMP/${comparator}/*.txt

	echo "deleting tmp directory.."

	rm -r ${PARASITE_SCRATCH2}/comparators_TMP/${comparator}/
	
	echo "Finished populating $db"
done


