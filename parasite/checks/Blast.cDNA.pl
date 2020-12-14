#!/usr/bin/env perl

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

my $dbh = $dba->dbc->db_handle;
my $sql = "SELECT COUNT(stable_id) FROM transcript;";
my $sql_1 = "SELECT meta_value FROM meta WHERE meta_key = 'species.url';";
my $sql_2 = "SELECT meta_value FROM meta WHERE meta_key = 'assembly.default';";

#my @row = $dbh->selectrow_array($sql)  ;
#print "@row\n";
my $sth = $dbh->prepare($sql);          # prepare the query
$sth->execute();                        # execute the query

my $sth_1 = $dbh->prepare($sql_1);          # prepare the query
$sth_1->execute();                        # execute the query

my $sth_2 = $dbh->prepare($sql_2);          # prepare the query
$sth_2->execute();                        # execute the query

#my @row;
my $value; 
$value = $sth->fetchrow_array; 

my $value1; 
$value1 = $sth_1->fetchrow_array; 

my $value2; 
$value2 = $sth_2->fetchrow_array; 

#while (@row = $sth->fetchrow_array) 
#{  
 #  print join(", ", @row), "\n";        # retrieve one row
#}
my $cDNA ; 
$cDNA = "$dp/$value1\.$value2\.cdna.all.fa" ;

if(! -e $cDNA)
{print "$db\tfile doesn't exit\n$cDNA\n";

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
 if ($value != $count_seq)
 {
print "$db\tcDNA_in_db:$value\tcDNA_in_dump:$count_seq\n";

 }
}
