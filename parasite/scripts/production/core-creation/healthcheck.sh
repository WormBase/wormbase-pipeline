#!/usr/bin/env bash


usage(){
cat << EOF
 
Usage: $0 -u ensro -s mysql-eg-devel-1.ebi.ac.uk -p 4126 -d ms41_echinococcus_canadensis_prjeb8992_core_10_88_1 -g GCA_900004735.1_ECANG7_genomic.fna -a EcanG7_V2_2.gff3 

Options:
 -h   Show this message
 -d   Database name
 -g   Full directory location of genome fasta file
 -a   Full directory location of annotation *gff3 file
 -f   Working Folder path
 -u   User
 -p   Port
 -s   MySQL server
EOF
}

pipeline_fold=`dirname $0`
current_fold='/nfs/panda/ensemblgenomes/wormbase/parasite/scripts/core-creation'

while getopts "hd:b:g:a:f:u:p:s:" OPTION
do
  case $OPTION in
    h) usage
      exit 1 ;;
    d) database_name=$OPTARG ;;
    b) spp_abb=$OPTARG ;;
    g) genome=$OPTARG ;;
    a) annotation=$OPTARG ;;
    f) current_fold=$OPTARG ;;    
    u) user=$OPTARG ;;
    p) port=$OPTARG ;;
    s) server=$OPTARG ;;
  esac
done

mysqlcmd="mysql -u$user -h$server -P$port"

echo "Check database $database_name upload."
echo "=== Assembly: ==="
echo "Dump and compare scaffolds from the database"

export DATA=$current_fold
mkdir -p $DATA/dump_out

echo "perl $ENSEMBL_CVS_ROOT_DIR/ensembl-analysis/scripts/sequence_dump.pl -dbuser $user -dbhost $server -dbport $port -dbname $database_name -coord_system_name scaffold -onefile -output_dir $DATA/dump_out"
perl $ENSEMBL_CVS_ROOT_DIR/ensembl-analysis/scripts/sequence_dump.pl -dbuser ensrw -dbuser $user -dbhost $server -dbport $port  -dbname $database_name -coord_system_name scaffold -onefile -output_dir $DATA/dump_out
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
echo "perl dump_translations.pl -dbname $database_name -dbhost mysql-eg-devel-1.ebi.ac.uk  -dbuser $user -dbhost $server -dbport $port -stable_id -file_err $DATA/dump_out/$database_name\_error.txt -file /dev/null"
perl $pipeline_fold/dump_translations.pl -dbname $database_name -dbuser $user -dbhost $server -dbport $port -stable_id -file_err $DATA/dump_out/$database_name\_error.fa -file /dev/null 

echo "The number of coding protein sequences is"
grep -c ">" $DATA/dump_out/$database_name\_protein.fa
echo "While the number of non coding sequences (with stop codon inside) are:"
grep -c ">" $DATA/dump_out/$database_name\_error.fa

echo "" 
echo "Ensure every gene has a canonical_transcript_id:" 
echo "1) the result from this mysql is empty" 
echo "$mysqlcmd -D $database_name -e 'select count(*) from gene where canonical_transcript_id = 0;'" 
$mysqlcmd -D $database_name -e 'select count(*) from gene where canonical_transcript_id = 0;'  
