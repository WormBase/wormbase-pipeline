#!/usr/bin/env perl -w

use DBI;
use strict;
use Getopt::Long;

my ($org_id, $ena_db, $uniprot_db, $verbose); 

&GetOptions(
  'enadb=s'         => \$ena_db,
  'uniprotdb=s'     => \$uniprot_db,
  'orgid=s'         => \$org_id,
  'verbose'         => \$verbose,
    );


if (not defined $uniprot_db or
    not defined $org_id or
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

# locus_tag = 84
# pseudo = 28
# pseudogene = 114
# standard name = 23
# gene = 12
# Get most info for each PID with a /locus_tag + featID was gene (#12)
my $ena_sql =  "SELECT d.primaryacc#, b.version, c.PROTEIN_ACC, c.version, c.chksum, fq.text, c.featid, d.project#, d.statusid"
    . " FROM cdsfeature c, dbentry d, bioseq b, feature_qualifiers fq"
    . " WHERE d.primaryacc# IN ("
    . "   SELECT primaryacc#"
    . "   FROM dbentry" 
    . "   JOIN sourcefeature USING (bioseqid)"
    . "   WHERE organism = $org_id"
    . "   AND project# = 1"
    . "   AND statusid = 4)"
    . " AND c.bioseqid = d.bioseqid"
    . " AND d.bioseqid = b.seqid"
    . " AND fq.featid  = c.featid"
    . " AND not exists (SELECT 1"
    . "   FROM feature_qualifiers a"
    . "   WHERE a.featid = c.featid"
    . "   AND a.fqualid = 114)"
    . " AND fq.fqualid   = 23";

    
my $ena_sth = $ena_dbh->prepare($ena_sql) or die "Can't prepare statement: $DBI::errstr";

print STDERR "Doing primary lookup of CDS entries in ENA ORACLE database...\n" if $verbose;

$ena_sth->execute or die "Can't execute statement: $DBI::errstr";

my (@resultsArr, %resultsHash);

while ( ( my @results ) = $ena_sth->fetchrow_array ) {
  push @{$resultsHash{$results[2]}}, {  
    NT_AC       => $results[0],
    NT_version  => $results[1],
    AA_PID      => $results[2],
    AA_version  => $results[3],
    AA_checksum => $results[4],
    AA_text     => $results[5],
    AA_ID       => $results[6],
  };
}
die $ena_sth->errstr if $ena_sth->err;
$ena_sth->finish;
$ena_dbh->disconnect;

#################
# query uniprot
#################

my $uniprot_dbh = &get_uniprot_dbh();

my $uniprot_sql =  "SELECT e.accession, e.name, p.protein_id "
      . "FROM sptr.dbentry e, sptr.embl_protein_id p "
      . "WHERE p.dbentry_id = e.dbentry_id "
      . "AND e.deleted='N' "
      . "AND e.merge_status <> 'R' "
      . "AND e.entry_type in ('0', '1') "
      . "AND e.tax_id = $org_id";


my $uniprot_sth = $uniprot_dbh->prepare($uniprot_sql);
print STDERR "Reading Uniprot database to get accessioans and ids...\n" if $verbose;
$uniprot_sth->execute();

while( (my @results) = $uniprot_sth->fetchrow_array) {
  if (exists $resultsHash{$results[2]}) {
    foreach my $el (@{$resultsHash{$results[2]}}) {
      $el->{SWALL_AC} = $results[0];
      $el->{SWALL_ID} = $results[1];
    }
  }
}

die $uniprot_sth->errstr if $uniprot_sth->err;
$uniprot_sth->finish;
$uniprot_dbh->disconnect;

##################
# output results
##################

@resultsArr = sort { $a->{NT_AC} cmp $b->{NT_AC} or $a->{AA_PID} cmp $b->{AA_PID} } map { @$_ } values(%resultsHash);
foreach my $entry (@resultsArr) {
  
  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
         $entry->{NT_AC},
         $entry->{NT_version},
         $entry->{AA_PID},
         $entry->{AA_version},
         $entry->{AA_checksum},
         $entry->{AA_text},
         exists($entry->{SWALL_AC}) ? $entry->{SWALL_AC} : "UNDEFINED",
         exists($entry->{SWALL_ID}) ? $entry->{SWALL_ID} : "UNDEFINED");

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


#####################
sub get_uniprot_dbh {
  print STDERR "Connecting to $uniprot_db\n" if $verbose;

   my $dbh = DBI->connect("dbi:Oracle:$uniprot_db", 
                         "spselect", 
                         "spselect", 
                          \%attr) 
      or die "Cannot connect to Uniprot database $DBI::errstr\n";
   return $dbh;
}
