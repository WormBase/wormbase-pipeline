#!/usr/bin/env perl -w

use DBI;
use strict;
use Getopt::Long;

my ($org_id, $ena_db, $verbose); 

&GetOptions(
  'enadb=s'         => \$ena_db,
  'orgid=s'         => \$org_id,
  'verbose'         => \$verbose,
    );


if (not defined $org_id or
    not defined $ena_db) {
  die "Incorrect invocation\n";
}


my %attr   = ( PrintError => 0,
               RaiseError => 0,
               AutoCommit => 0 );


##############
# Query ENA
#############

my $ena_dbh = &get_ena_dbh();

my $ena_sql =  "SELECT d.primaryacc#, b.version"
    . " FROM dbentry d, bioseq b"
    . " WHERE d.primaryacc# IN ("
    . "   SELECT primaryacc#"
    . "   FROM dbentry" 
    . "   JOIN sourcefeature USING (bioseqid)"
    . "   WHERE organism = $org_id"
    . "   AND project# = 1"
    . "   AND statusid = 4)"
    . " AND d.bioseqid = b.seqid";
    
my $ena_sth = $ena_dbh->prepare($ena_sql);

print STDERR "Doing primary lookup of CDS entries in ENA ORACLE database...\n" if $verbose;

$ena_sth->execute || die "Can't execute statement: $DBI::errstr";

my (%results);

while ( ( my @results ) = $ena_sth->fetchrow_array ) {
  $results{$results[0]} = $results[1];
}

die $ena_sth->errstr if $ena_sth->err;
$ena_sth->finish;
$ena_dbh->disconnect;

foreach my $acc (sort keys %results) {
  printf("%s\t%s\n", $acc, $results{$acc});
}


exit(0);



#####################
sub get_ena_dbh {

  my $dbh = DBI->connect("dbi:Oracle:$ena_db", 
                         'ena_reader', 
                         'reader', 
                         \%attr)
      or die "Can't connect to ENA database: $DBI::errstr";

  return $dbh;
}

