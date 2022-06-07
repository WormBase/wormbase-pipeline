#Quick script to generate alpha_mappings.txt files for species with AlphaFold data
#Nishadi De Silva, Feb 2022

use strict;
use warnings;
use JSON;
#use Data::Dumper;
use Bio::EnsEMBL::MetaData::DBSQL::MetaDataDBAdaptor;
use Cwd;


#wget http://ftp.ebi.ac.uk/pub/databases/alphafold/download_metadata.json
#perl generate_mapping_files.pl download_metadata.json

#Default options
my $metadata_file = shift;
my $afdb_path = "/nfs/ftp/public/databases/alphafold/latest/";
my $working_dir = "$ENV{PARASITE_SCRATCH}/alphafold/WBPS$ENV{PARASITE_VERSION}";
my $uniprot_api = "https://rest.uniprot.org/beta/proteomes";
my $tmp_dir = "${working_dir}/tmp";

#Exclude these species. They are 1) massive and 2) already done by mag
my @exclude_species = ("arabidopsis_thaliana" , "danio_rerio"  ,"glycine_max" , "homo_sapiens"  ,"mus_musculus"  ,"rattus_norvegicus","zea_mays");

#Read metadata file provided by AFDB
open my $FH, "<", $metadata_file or die "could not open $metadata_file \n";
my $json = JSON->new;
my $data = $json->decode(<$FH>);

#Connect to Ensembl metadata DB
my $metadatadba=Bio::EnsEMBL::MetaData::DBSQL::MetaDataDBAdaptor->new(
						-group => 'core',
						-dbname => 'ensembl_metadata_'.$ENV{ENSEMBL_VERSION},
						-host => 'mysql-ps-staging-1',
						-port => '4451',
						-user => 'ensro',
						-driver => 'mysql');
#Get metadata adaptors
my $gcdba = $metadatadba->get_GenomeOrganismInfoAdaptor();

#for each afdb species, find uniprot proteome ID, extract GCA and create alpha_mappings file (if assembly exists in ensembl)

foreach (@$data){

	my $ref_proteome = $_->{'reference_proteome'}; #uniprot proteome ID
	my $archive_name = $_->{'archive_name'}; #name of tar file in AlphaFold

	if($ref_proteome) {
        mkdir $tmp_dir unless (-d $tmp_dir);
		if(!-e "${tmp_dir}/${ref_proteome}.json"){
			`wget ${uniprot_api}/${ref_proteome}.json -P $tmp_dir`;
		}
		my $uniprot_json = JSON->new;
		open my $UNIPROT_FH, "<", "${tmp_dir}/${ref_proteome}.json" or die "Could not open uniprot file \n";
		my $uniprot_file = $uniprot_json->decode(<$UNIPROT_FH>);
		my $assembly_details = $uniprot_file->{'taxonomy'};
		my $gca = $assembly_details->{'scientificName'};
		my $genome = $gcdba->fetch_by_scientific_name($gca);
		if($genome){
			my $prod_name = $genome->name();
			if ( !grep( /^$prod_name$/, @exclude_species ) && !-e "${working_dir}/${prod_name}/alpha_mappings.txt" ) { #Currently does not overwrite existing mappings
				print "Working on: ", $prod_name, "\n";
				mkdir "${working_dir}/${prod_name}" unless (-d "${working_dir}/${prod_name}");
				`tar xaf $afdb_path/$archive_name *.pdb.gz --to-command='zgrep --line-buffered DBREF | sed -e s,\$"DBREF ","\${TAR_FILENAME%.*}:DBREF",g' > ${working_dir}/${prod_name}/alpha_mappings.txt`;
			}
		}
        unlink "$working_dir/tmp/$ref_proteome.json";
	}
}
