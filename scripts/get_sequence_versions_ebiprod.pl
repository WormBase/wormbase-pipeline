#!/usr/bin/env perl -w

use DBI;
use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my ($bioproject_id, $ena_cred, $verbose, $all_records); 

&GetOptions(
  'enacred=s'       => \$ena_cred,
  'bioprojectid=s'  => \$bioproject_id, 
  'verbose'         => \$verbose,
  'nonpublic'       => \$all_records,
    );


if (not defined $bioproject_id) {
  die "Incorrect invocation\n";
}


##############
# Query ENA
#############

my $ena_dbh = &get_ena_dbh($ena_cred);

#
# statusid = 4 => only live public records
#
my $ena_sql = 
    "   SELECT primaryacc#, b.version"
    . "   FROM dbentry d, bioseq b"
    . "   WHERE d.bioseqid = b.seqid"
    . "   AND study_id = '$bioproject_id'"
    . "   AND dataclass <> 'SET'";
if (not $all_records) {
  $ena_sql .=  "   AND statusid = 4";
}

my $ena_sth = $ena_dbh->dbc->prepare($ena_sql);

print STDERR "Doing primary lookup of sequence versions in ENA ORACLE database...\n" if $verbose;

$ena_sth->execute || die "Cannot execute statement: $DBI::errstr";

my (%results);

while ( ( my @results ) = $ena_sth->fetchrow_array ) {
  $results{$results[0]} = $results[1];
}

die $ena_sth->errstr if $ena_sth->err;
$ena_sth->finish;
$ena_dbh->dbc->disconnect_if_idle;

foreach my $acc (sort keys %results) {
  printf("%s\t%s\n", $acc, $results{$acc});
}


exit(0);



#####################
sub get_ena_dbh {
  my ($cred_file) = @_;

  my ($dbname, $user, $pass, $host, $port);

  open(my $cfh, $cred_file) or die "Could not open $cred_file for reading\n";
  while(<$cfh>) {
    /^\#/ and next;

    /^DBNAME:(\S+)/ and $dbname = $1;
    /^USER:(\S+)/ and $user = $1;
    /^PASS:(\S+)/ and $pass = $1;
    #/^HOST:(\S+)/ and $host = $1;
    #/^PORT:(\S+)/ and $port = $1;    
  }

  my $dbh = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    #-host => $host,
    #-port => $port,
    -user => $user,
    -pass => $pass,
    -dbname => $dbname,
    -driver => 'Oracle');

  return $dbh;
}

