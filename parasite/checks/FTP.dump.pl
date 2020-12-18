#!/usr/bin/env perl

use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Bio::SeqIO;
use File::Find;


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

my $dbh = $dba->dbc->db_handle;
my $sql_dna = "SELECT COUNT(distinct(name)) FROM seq_region;";# count number of dna in database
my $sql_pep = "SELECT COUNT(translation_id) FROM translation;";# count number of protein in database
my $sql_trans = "SELECT COUNT(stable_id) FROM transcript;"; # count number of transcript in database
my $sql_1 = "SELECT meta_value FROM meta WHERE meta_key = 'species.url';";
my $sql_2 = "SELECT meta_value FROM meta WHERE meta_key = 'species.bioproject_id';"; 

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



my @val = split (/_/, $value1);
my $spe_name = lc("$val[0]_$val[1]");


my ( $pep,$dna_genomic, $dna_genomicSM, $dna_genomicM,$trans  );
$pep = "$dp/$spe_name/$value2/$spe_name\.$value2\.WBPS15.protein.fa.gz" ;
$dna_genomic = "$dp/$spe_name/$value2/$spe_name\.$value2\.WBPS15.genomic.fa.gz" ;
$dna_genomicM = "$dp/$spe_name/$value2/$spe_name\.$value2\.WBPS15.genomic_masked.fa.gz" ;
$dna_genomicSM = "$dp/$spe_name/$value2/$spe_name\.$value2\.WBPS15.genomic_softmasked.fa.gz" ;
$trans = "$dp/$spe_name/$value2/$spe_name\.$value2\.WBPS15.mRNA_transcripts.fa.gz" ;

#working on protein

if(! -e $pep)
{
print "$db\tfile doesn't exit\n";
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
 

 #{
print "$db\tpep_in_db:$value_pep\tpep_in_dump:$count_seq\n";

 #}
}

#working on genomic DNA

if(! -e $dna_genomic)
{
print "$db\tfile doesn't exit\n";
}
else

{

my $fa = Bio::SeqIO ->new ( '-format' =>'fasta' , '-file' => "$dna_genomic" );
my $count_seq=0;
while (my $seq = $fa->next_seq())
 {
    my $sequence = $seq->seq();
    my $id       = $seq->display_id();
    $count_seq++;
 }
 
#if ($value != $count_seq)
 #{
print "$db\tpep_in_db:$value_dna\tpep_in_dump:$count_seq\n";

 #}
}


#working on masked genomic DNA

if(! -e $dna_genomicM)
{
print "$db\tfile doesn't exit\n";
}
else

{

my $fa = Bio::SeqIO ->new ( '-format' =>'fasta' , '-file' => "$dna_genomicM" );
my $count_seq=0;
while (my $seq = $fa->next_seq())
 {
    my $sequence = $seq->seq();
    my $id       = $seq->display_id();
    $count_seq++;
 }
 
#if ($value != $count_seq)
 #{
print "$db\tpep_in_db:$value_dna\tpep_in_dump:$count_seq\n";

 #}
}

#working on soft masked genomic DNA

if(! -e $dna_genomicSM)
{
print "$db\tfile doesn't exit\n";
}
else

{

my $fa = Bio::SeqIO ->new ( '-format' =>'fasta' , '-file' => "$dna_genomicSM" );
my $count_seq=0;
while (my $seq = $fa->next_seq())
 {
    my $sequence = $seq->seq();
    my $id       = $seq->display_id();
    $count_seq++;
 }
 
#if ($value != $count_seq)
 #{
print "$db\tpep_in_db:$value_dna\tpep_in_dump:$count_seq\n";

 #}
}

#working on soft transcript

if(! -e $$trans)
{
print "$db\tfile doesn't exit\n";
}
else

{

my $fa = Bio::SeqIO ->new ( '-format' =>'fasta' , '-file' => "$trans" );
my $count_seq=0;
while (my $seq = $fa->next_seq())
 {
    my $sequence = $seq->seq();
    my $id       = $seq->display_id();
    $count_seq++;
 }
 
#if ($value != $count_seq)
 #{
print "$db\tpep_in_db:$value_trans\tpep_in_dump:$count_seq\n";

 #}
}