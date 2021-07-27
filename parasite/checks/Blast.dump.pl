#!/usr/bin/env perl
# run the script as this: core _${PARASITE_VERSION}_ | while read -r i ; do perl Blast.dump.pl $($PARASITE_STAGING_MYSQL details script) 
# -dumpPath /hps/nobackup/production/ensemblgenomes/parasite/production/dumps/WBPS15/BLAST -dbname $i ; done
use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Bio::SeqIO;

#my $input = shift or die ; 
my ($host, $user, $db, $pass, $port,$dp, $dump_file);
GetOptions
(
'host=s'       => \$host,
'dbname=s'     => \$db,
'pass=s'       => \$pass,
'port=i'       => \$port,
'user=s'       => \$user,
'db=s'       => \$db,
'dumpPath=s'   => \$dp,
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
# connect to database
my $dbh = $dba->dbc->db_handle;
my $sql_dna = "SELECT COUNT(distinct(name)) FROM seq_region WHERE coord_system_id = 1;";# count number of DNA in database
my $sql_pep = "SELECT COUNT(translation_id) FROM translation;";# count number of protein in database
my $sql_trans = "SELECT COUNT(stable_id) FROM transcript;"; # count number of transcript in database
my $sql_1 = "SELECT meta_value FROM meta WHERE meta_key = 'species.url';"; # retreive the species.url from meta table
my $sql_2 = "SELECT meta_value FROM meta WHERE meta_key = 'assembly.default';";

#my @row = $dbh->selectrow_array($sql)  ;
#print "@row\n";
my $sth_trans = $dbh->prepare($sql_trans);          # prepare the query
$sth_trans->execute();                        # execute the query

my $sth_pep = $dbh->prepare($sql_pep);          # prepare the query
$sth_pep->execute();                        # execute the query

my $sth_dna = $dbh->prepare($sql_dna);          # prepare the query
$sth_dna->execute();                        # execute the query


my $sth_1 = $dbh->prepare($sql_1);          # prepare the query
$sth_1->execute();                        # execute the query

my $sth_2 = $dbh->prepare($sql_2);          # prepare the query
$sth_2->execute();                        # execute the query

#my @row;
# parse the query to variable 
my $value_trans; 
$value_trans = $sth_trans->fetchrow_array; 

my $value_pep; 
$value_pep = $sth_pep->fetchrow_array; 

my $value_dna; 
$value_dna = $sth_dna->fetchrow_array; 

my $value1; 
$value1 = $sth_1->fetchrow_array; 

my $value2; 
$value2 = $sth_2->fetchrow_array; 

#while (@row = $sth->fetchrow_array) 
#{  
 #  print join(", ", @row), "\n";        # retrieve one row
#}

# preparing path for dump files
my $cDNA ; 
$cDNA = "$dp/$value1\.$value2\.cdna.all.fa" ;

my $pep ; 
$pep = "$dp/$value1\.$value2\.pep.all.fa" ;

my $dna_rm ; 
$dna_rm = "$dp/$value1\.$value2\.dna_rm.toplevel.fa" ;

my $dna ; 
$dna = "$dp/$value1\.$value2\.dna.toplevel.fa" ;

# match the number of cDNA in dump and database



if(! -e $cDNA)
{
   print "$db\tfile doesn't exit\n$cDNA\n";

}
else
{
my $fa = Bio::SeqIO ->new ( '-format' =>'fasta' , '-file' => "$cDNA" );
my $count_seq=0;
while (my $seq = $fa->next_seq())
 {
    my $sequence = $seq->seq();
    my $id       = $seq->display_id();
    $count_seq++;
 }
 if ($value_trans != $count_seq)
 {
print "$db\tcDNA_in_db:$value_trans\tcDNA_in_dump:$count_seq\n";

 }
}

# match the number of protein in dump and database


if(! -e $pep)
{
   print "$db\tfile doesn't exit\n"

}
else

{
my $fa = Bio::SeqIO ->new ( '-format' =>'fasta' , '-file' => "$pep" );
my $count_seq=0;
while (my $seq = $fa->next_seq())
 {
    my $sequence = $seq->seq();
    my $id       = $seq->display_id();
    $count_seq++;
 }
 
if ($value_pep != $count_seq)
 {
print "$db\tpep_in_db:$value_pep\tpep_in_dump:$count_seq\n";

 }
}

# match the number of DNA rm in dump and database



if(! -e $dna)
{print "$db\tfile doesn't exit\n"

}
else

{


my $fa = Bio::SeqIO ->new ( '-format' =>'fasta' , '-file' => "$dna_rm" );
my $count_seq=0;
while (my $seq = $fa->next_seq())
 {
    my $sequence = $seq->seq();
    my $id       = $seq->display_id();
    $count_seq++;
 }
 
 if ($value_dna != $count_seq)
 {
print "$db\tDNA_rm_in_db:$value_dna\tDNA_rm_in_dump:$count_seq\n";

 }
}

# match the number of DNA in dump and database



if(! -e $dna)
{print "$db\tfile doesn't exit\n"

}
else

{


my $fa = Bio::SeqIO ->new ( '-format' =>'fasta' , '-file' => "$dna" );
my $count_seq=0;
while (my $seq = $fa->next_seq())
 {
    my $sequence = $seq->seq();
    my $id       = $seq->display_id();
    $count_seq++;
 }
 
 if ($value_dna != $count_seq)
 {
print "$db\tDNA_in_db:$value_dna\tDNA_in_dump:$count_seq\n";

 }
}
