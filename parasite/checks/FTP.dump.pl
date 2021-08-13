
# core _${PARASITE_VERSION}_ | while read -r i ; 
# do perl FTP.dump.pl $($PARASITE_STAGING_MYSQL details script) 
#  -dumpPath /hps/nobackup/production/ensemblgenomes/parasite/production/dumps/WBPS15/FTP
# -dbname $i -wbps_version $PARASITE_VERSION ; done


#!/usr/bin/env perl
use IO::Uncompress::Gunzip;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Bio::SeqIO;
use File::Find;


#my $input = shift or die ; 
my ($host, $user, $db, $pass, $port,$dp, $dump_file, $wbps_version);
GetOptions
(
'host=s'       => \$host,
'dbname=s'     => \$db,
'pass=s'       => \$pass,
'port=i'       => \$port,
'user=s'       => \$user,
'db=s'       => \$db,
'dumpPath=s'   => \$dp,
'wbps_version=i' => \$wbps_version,
#'dump_file=s'  => \$dump_file,
)|| die ("check command line params\n");

my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor
(
-host   => $host,
-user   => $user,
-dbname => $db,
-pass   => $pass,
-port   => $port,
);

my $dbh = $dba->dbc->db_handle;
my $sql_dna = "SELECT COUNT(*) FROM seq_region LEFT JOIN seq_region_attrib
               ON seq_region.seq_region_id = seq_region_attrib.seq_region_id 
	       LEFT JOIN attrib_type
	       ON attrib_type.attrib_type_id = seq_region_attrib.attrib_type_id
               WHERE attrib_type.code = 'toplevel';"; # count number of dna in database
my $sql_pep = "SELECT COUNT(translation_id) FROM translation;";# count number of protein in database
my $sql_trans = "SELECT COUNT(stable_id) FROM transcript;"; # count number of transcript in database
my $sql_1 = "SELECT meta_value FROM meta WHERE meta_key = 'species.url';";

# sub routine to count the sequence in dump
sub mycount 
{
my $file_path = shift ;
my $fa = Bio::SeqIO ->new ( '-format' =>'fasta' , '-file' => "$file_path" );
my $count_seq=0;
while (my $seq = $fa->next_seq())
 {
    my $sequence = $seq->seq();
    my $id       = $seq->display_id();
    $count_seq++;
    
 }
return $count_seq ; 
}


my $sth_trans = $dbh->prepare($sql_trans);          # prepare the query
$sth_trans->execute();                        # execute the query

my $sth_pep = $dbh->prepare($sql_pep);          # prepare the query
$sth_pep->execute();                        # execute the query

my $sth_dna = $dbh->prepare($sql_dna);          # prepare the query
$sth_dna->execute();                        # execute the query


my $sth_1 = $dbh->prepare($sql_1);          # prepare the query
$sth_1->execute();                        # execute the query


# parse the query to variable 
my $value_trans; 
$value_trans = $sth_trans->fetchrow_array; 

my $value_pep; 
$value_pep = $sth_pep->fetchrow_array; 

my $value_dna; 
$value_dna = $sth_dna->fetchrow_array; 

my $value1; 
$value1 = $sth_1->fetchrow_array; 


my @val = split (/_/, $value1);
my $spe_name = lc("$val[0]_$val[1]");
my $bioproject = uc($val[2]);
$bioproject =~ s/(.)PRJ/$1_PRJ/;

# Getting path for different files in dump in variable for both ziped and unzipped file as I have mix of them
my ( $pepgz,$dna_genomicgz, $dna_genomicSMgz, $dna_genomicMgz,$transgz, $pep,$dna_genomic, $dna_genomicSM, $dna_genomicM,$trans);
my $path = "$dp/$spe_name/$bioproject/$spe_name.$bioproject.WBPS$wbps_version";

$pepgz = "$path.protein.fa.gz" ;
$dna_genomicgz = "$path.genomic.fa.gz" ;
$dna_genomicMgz = "$path.genomic_masked.fa.gz" ;
$dna_genomicSMgz = "$path.genomic_softmasked.fa.gz" ;
$transgz = "$path.mRNA_transcripts.fa.gz" ;

#path for files 

my $pep_path = "$path.protein.fa";
my $genomic_path = "$path.genomic.fa";
my $genomicM_path = "$path.genomic_masked.fa";
my $genomicSM_path = "$path.genomic_softmasked.fa";
my $trans_path = "$path.mRNA_transcripts.fa";


print "working on: $db\n";
#working on protein
my ($protein_result, $genomic_result, $genomicSM_result, $genomicM_result, $trans_result); 

if (-e $pepgz)
{
system("gunzip $pepgz");
$protein_result = mycount ($pep_path) ;           
if ($value_pep != $protein_result)
 {
print "$db\tpep_in_db:$value_pep\tpep_in_dump:$protein_result\n";

 }

}

else {"$db\tprotein file is not found\n";}


system("gzip $pep_path");
#print "$db\tpep_in_db:$value_pep\tpep_in_dump:$protein_result\n";





#working on genomic DNA
if (-e $dna_genomicgz)
{
system("gunzip $dna_genomicgz");
$genomic_result = mycount ($genomic_path) ;
if ($value_dna != $genomic_result)
 {
print "$db\tgenomeDNA_in_db:$value_dna\tgenomeDNA_in_dump:$genomic_result\n";

 }
}

else {"$db\tgenomic DNA file is not found\n";}

system("gzip $genomic_path");



#working on masked genomic DNA

if (-e $dna_genomicMgz)
{
system("gunzip $dna_genomicMgz");
$genomicM_result = mycount ($genomicM_path) ;
if ($value_dna != $genomicM_result)
 {
print "$db\tgenomeDNA_in_db:$value_dna\tgenomeDNAmasked_in_dump:$genomicM_result\n";

 }
}

else {"$db\tgenomic Masked DNA file is not found\n";}

system("gzip $genomicM_path");





#working on soft masked genomic DNA

if (-e $dna_genomicSMgz)
{
system("gunzip $dna_genomicSMgz");
$genomicSM_result = mycount ($genomicSM_path) ;
if ($value_dna != $genomicSM_result)
 {
print "$db\tgenomeDNA_in_db:$value_dna\tgenomesoftmaskedDNA_in_dump:$genomicSM_result\n";

 }


}

else {"$db\tgenomic SoftMasked DNA file is not found\n";}

system("gzip $genomicSM_path");


#working on soft transcript

if (-e $transgz)
{
system("gunzip $transgz");
$trans_result = mycount ($trans_path) ;

if ($value_trans != $trans_result)
 {
print "$db\tmRNA_in_db:$value_trans\tmRNA_in_dump:$trans_result\n";

 }
}

else {"$db\tgenomic SoftMasked DNA file is not found\n";}

system("gzip $trans_path");

