#!/usr/bin/env bash


usage(){
cat << EOF
 
Usage: $0 -d at8_species_name_version -s species_name -b species_abbreviation -g directory/path/genome/FASTA/file -a directory/path/annotation/GFF3/file -t nnnnn -x 2012-11 -v 1.0

This script loads the specified genome into Ensembl
and then runs GenBlastg Runnable
You need to know the latest WormBase release number, check here:
/nfs/wormpub/BUILD/WORMPEP/wormpep...

Options:
 -h   Show this message
 -d   Database name
 -e   Ensembl branch folder
 -b   Species abbreviation (tag)
 -g   Full directory location of genome fasta file
 -a   Full directory location of annotation *gff3 file
 -f   Working Folder path

EOF
}

pipeline_fold=`dirname $0`
current_fold='/nfs/panda/ensemblgenomes/wormbase/parasite/scripts/core-creation'

while getopts "hd:b:g:a:f`dirname $0`:e:" OPTION
do
  case $OPTION in
    h) usage
      exit 1 ;;
    d) database_name=$OPTARG ;;
    b) spp_abb=$OPTARG ;;
    g) genome=$OPTARG ;;
    a) annotation=$OPTARG ;;
    f) current_fold=$OPTARG ;;    
    e) GIT=$OPTARG ;;
    ?) usage
      exit ;;
  esac
done


mysqlcmd="mysql -uensro -hmysql-eg-devel-1.ebi.ac.uk -P4126"


echo "Check database $database_name upload."
echo "=== Assembly: ==="
echo "Dump and compare scaffolds from the database"

export DATA=$current_fold
export ENSEMBL_ANALYSIS=$GIT/ensembl-analysis/ 
export AS=$ENSEMBL_ANALYSIS/scripts
mkdir -p $DATA/dump_out

#echo "perl $AS/sequence_dump.pl -dbuser ensrw -dbhost mysql-ps-staging-1 -dbport 4451 -dbpass scr1b3PSs1 -dbname $database_name -coord_system_name scaffold -onefile -output_dir $DATA/dump_out"
#perl $AS/sequence_dump.pl -dbuser ensrw -dbhost mysql-eg-devel-1.ebi.ac.uk  -dbport 4126 -dbpass scr1b3PSs1 -dbname $database_name -coord_system_name scaffold -onefile -output_dir $DATA/dump_out
echo "perl $AS/sequence_dump.pl -dbuser ensrw -dbhost mysql-eg-devel-1.ebi.ac.uk -dbport 4126 -dbpass scr1b3d1  -dbname $database_name -coord_system_name scaffold -onefile -output_dir $DATA/dump_out"
perl $AS/sequence_dump.pl -dbuser ensrw -dbhost mysql-eg-devel-1.ebi.ac.uk  -dbport 4126 -dbpass scr1b3d1  -dbname $database_name -coord_system_name scaffold -onefile -output_dir $DATA/dump_out
mv $DATA/dump_out/scaffold.fa $DATA/dump_out/$database_name\_scaffold.fa

echo "The number of sequences in database $database_name are: " 
grep -c '^>' $DATA/dump_out/$database_name\_scaffold.fa

echo "The number of sequences in fasta file are: " 
grep -c '^>' $genome 

echo "" 
echo "Check $database_name annotation upload" 
echo "Number of genes to be uploaded from the GFF file" 
echo `awk '{if ($3 == "gene")print $0}' $annotation | grep -c ""` 

echo "" 
echo "Number of genes uploaded in $database_name" 
echo "$mysqlcmd -D $database_name -e 'select count(*) from gene;'" 
$mysqlcmd -D $database_name -e 'select count(*) from gene;'  

echo ""
echo "Number of transcripts to be uploaded from the GFF file"
echo `awk '{if ($3 == "transcript")print $0}' $annotation | grep -c ""`
echo "Number of mRNAs to be uploaded from the GFF file"
echo `awk '{if ($3 == "mRNA")print $0}' $annotation | grep -c ""`

echo ""
echo "Total number of transcripts (transcript+mRNA) uploaded in $database_name"
echo "$mysqlcmd -D $database_name -e 'select count(*) from transcript;'"
$mysqlcmd -D $database_name -e 'select count(*) from transcript;'

echo "" 
echo "Final check:" 
echo "Protein sequences: ensure these match the one used to load (taking into account the biotype protein_coding vs non_coding_cds):"  			
echo "$mysqlcmd -D $database_name -e 'select biotype, count(*) from gene group by biotype;'" 
echo "Endend the upload of assembly and annotation in $database_name" 
$mysqlcmd -D $database_name -e 'select biotype, count(*) from gene group by biotype;' 

echo ""
echo "Dump protein fasta file:"
#echo "perl $GIT/ensembl-analysis/scripts/protein/dump_translations.pl -dbname $database_name -dbhost mysql-ps-staging-1 -dbport 4451 -dbpass scr1b3PSs1 -dbuser ensrw -stable_id 2> $DATA/dump_out/$database_name\_error.txt > $DATA/dump_out/$database_name\_protein.fa"
#perl $GIT/ensembl-analysis/scripts/protein/dump_translations.pl -dbname $database_name -dbhost mysql-ps-staging-1 -dbport 4451 -dbpass scr1b3PSs1 -dbuser ensrw -stable_id 2> $DATA/dump_out/$database_name\_error.txt > $DATA/dump_out/$database_name\_protein.fa
echo "perl dump_translations.pl -dbname $database_name -dbhost mysql-eg-devel-1.ebi.ac.uk  -dbport 4126 -dbpass scr1b3d1 -dbuser ensrw -stable_id -file_err $DATA/dump_out/$database_name\_error.txt -file $DATA/dump_out/$database_name\_protein.fa"
perl $pipeline_fold/dump_translations.pl -dbname $database_name -dbhost mysql-eg-devel-1.ebi.ac.uk  -dbport 4126 -dbpass scr1b3d1  -dbuser ensrw -stable_id -file_err $DATA/dump_out/$database_name\_error.fa -file $DATA/dump_out/$database_name\_protein.fa

echo "The number of coding protein sequences is"
grep -c ">" $DATA/dump_out/$database_name\_protein.fa
echo "While the number of non coding sequences (with stop codon inside) are:"
grep -c ">" $DATA/dump_out/$database_name\_error.fa

echo ""
echo "Dump transcript fasta file:"
#echo "perl $WA/src/ensembl-other/scripts/dump_proteins.pl -database $database_name -transcript_only > $DATA/dump_out/$database_name\_transcript.fa"
#echo "perl $WORM_CODE/scripts/ENSEMBL/scripts/dump_proteins.pl -dbhost mysql-ps-staging-1 -dbport 4451 -dbuser ensro —dbpass scr1b3PSs1 -dbname $database_name -outfile $DATA/dump_out/$database_name\_transcript.fa -transcript_only";
#perl $WORM_CODE/scripts/ENSEMBL/scripts/dump_proteins.pl -dbhost mysql-ps-staging-1 -dbport 4451 -dbuser ensro —dbpass scr1b3PSs1 -dbname $database_name -outfile $DATA/dump_out/$database_name\_transcript.fa -transcript_only

echo "perl $WORM_CODE/scripts/ENSEMBL/scripts/dump_proteins.pl -host mysql-eg-devel-1.ebi.ac.uk  -port 4126 -user ensro —pass scr1b3d1  -dbname $database_name -outfile $DATA/dump_out/$database_name\_transcript.fa -transcript_only";
perl $WORM_CODE/scripts/ENSEMBL/scripts/dump_proteins.pl -host mysql-eg-devel-1.ebi.ac.uk  -port 4126 -user ensro —pass scr1b3d1  -dbname $database_name -outfile $DATA/dump_out/$database_name\_transcript.fa -transcript_only

#perl $WA/src/ensembl-other/scripts/dump_proteins.pl -database $database_name -transcript_only > $DATA/dump_out/$database_name\_transcript.fa


echo "" 
echo "Ensure every gene has a canonical_transcript_id:" 
echo "1) the result from this mysql is empty" 
echo "$mysqlcmd -D $database_name -e 'select count(*) from gene where canonical_transcript_id = 0;'" 
$mysqlcmd -D $database_name -e 'select count(*) from gene where canonical_transcript_id = 0;'  
