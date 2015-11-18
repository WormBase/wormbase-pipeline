#!/usr/bin/env perl -w

use DBI;
use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my ($org_id, $ena_cred, $uniprot_cred, $bioproject_id, $verbose); 

&GetOptions(
  'orgid=s'         => \$org_id,
  'bioprojectid=s'  => \$bioproject_id, 
  'enacred=s'       => \$ena_cred,
  'uniprotcred=s'   => \$uniprot_cred,
  'verbose'         => \$verbose,
    );


if (not defined $bioproject_id or
    not defined $ena_cred or
    not defined $uniprot_cred) {
  die "Incorrect invocation: you must supply -orgid, -enacred and -uniprotcred\n";
}


##############
# Query ENA
#############

#my $test_uniprot_dbh = &get_uniprot_dbh($uniprot_cred);

my $ena_dbh = &get_ena_dbh($ena_cred);

# locus_tag = 84
# pseudo = 28
# pseudogene = 114
# standard name = 23
# gene = 12
# Get most info for each PID with a /locus_tag + featID was gene (#12)
# statusid = 4 => public/active records only (i.e. not supressed)
my $ena_sql =  "SELECT d.primaryacc#, b.version, c.PROTEIN_ACC, c.version, c.chksum, fq.text, c.featid, d.project#, d.statusid"
    . " FROM cdsfeature c, dbentry d, bioseq b, feature_qualifiers fq"
    . " WHERE d.primaryacc# IN ("
    . "   SELECT primaryacc#"
    . "   FROM dbentry" 
    . "   WHERE study_id = '$bioproject_id'"
    . "   AND statusid = 4)"
    . " AND c.bioseqid = d.bioseqid"
    . " AND d.bioseqid = b.seqid"
    . " AND fq.featid  = c.featid"
    . " AND not exists (SELECT 1"
    . "   FROM feature_qualifiers a"
    . "   WHERE a.featid = c.featid"
    . "   AND a.fqualid = 114)"
    . " AND fq.fqualid   = 23";

    
my $ena_sth = $ena_dbh->dbc->prepare($ena_sql) or die "Can't prepare statement: $DBI::errstr";

print STDERR "Doing primary lookup of CDS entries in ENA ORACLE database...\n" if $verbose;

$ena_sth->execute or die "Can't execute statement: $DBI::errstr";

my (@resultsArr, %resultsHash, %uniparc_mapping);

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
$ena_dbh->dbc->disconnect_if_idle;

#################
# query uniprot 1 - get basic protein_id info
#################

my $uniprot_dbh = &get_uniprot_dbh($uniprot_cred);

my $uniprot_sql =  "SELECT e.accession, e.name, p.protein_id "
      . "FROM sptr.dbentry e, sptr.embl_protein_id p "
      . "WHERE p.dbentry_id = e.dbentry_id "
      . "AND e.deleted='N' "
      . "AND e.merge_status <> 'R' "
      . "AND e.entry_type in ('0', '1') "
      . "AND e.tax_id = $org_id";


my $uniprot_sth = $uniprot_dbh->dbc->prepare($uniprot_sql);
print STDERR "Reading Uniprot database to get accessions and ids...\n" if $verbose;
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
$uniprot_dbh->dbc->disconnect_if_idle;

##################
# query uniprot 2 - use UniParc to get Uniprot isoform identifiers
###################

$uniprot_dbh = &get_uniprot_dbh($uniprot_cred);

# dbid = 1 => ENA protein_id
# dbid = 24 => Uniprot isoform id
$uniprot_sql = 'SELECT e.UPI, e.AC, e.DBID '
    . 'FROM uniparc.xref@uapro e '
    . "WHERE taxid = $org_id "
    . "AND  deleted = 'N' "
    . 'AND  dbid in (1,24)';

$uniprot_sth = $uniprot_dbh->dbc->prepare($uniprot_sql);
print STDERR "Reading UniParc for mapping of ENA protein ids to UP isoform ids\n" if $verbose;
$uniprot_sth->execute;
while (my @results = $uniprot_sth->fetchrow_array) {
  my ($upi, $acc, $db_id) = @results;
  if ($db_id == 1) {
    # ena protein id; ignore if it's not one of ours
    if (exists $resultsHash{$acc}) {
      $uniparc_mapping{$upi}->{protein_id} = $acc;
    }
  } else {
    $uniparc_mapping{$upi}->{isoform_id} = $acc;
  }
}
die $uniprot_sth->errstr if $uniprot_sth->err;
$uniprot_sth->finish;
$uniprot_dbh->dbc->disconnect_if_idle;

foreach my $upi (keys %uniparc_mapping) {
  if (exists $uniparc_mapping{$upi}->{protein_id} and 
      exists $uniparc_mapping{$upi}->{isoform_id}) {
    my ($pid) = $uniparc_mapping{$upi}->{protein_id};
    my ($iid) = $uniparc_mapping{$upi}->{isoform_id};
    foreach my $el (@{$resultsHash{$pid}}) {
      $el->{SWALL_ISOFORM} = $iid;
    }
  }
}

##################
# query uniprot 2 - get EC numbers
###################

$uniprot_dbh = &get_uniprot_dbh($uniprot_cred);
#
# subcategory_type_id = 3 => EC number
#
$uniprot_sql = "SELECT e.accession, cs.catg_type, cs.subcatg_type, ds.descr "
    . "FROM dbentry e "
    . "  JOIN dbentry_2_desc ds on (ds.dbentry_id = e.dbentry_id) "
    . "  JOIN cv_desc cs on (cs.desc_id = ds.desc_id) "
    . "WHERE cs.catg_type = 'RecName' " 
    . "AND cs.subcatg_type in ('EC', 'Full') "
    . "AND e.accession = ?";

$uniprot_sth = $uniprot_dbh->dbc->prepare($uniprot_sql);
print STDERR "Reading UniParc for mapping of ENA protein ids to EC numbers\n" if $verbose;
foreach my $entry_a (values %resultsHash) {
  foreach my $entry (@$entry_a) {
    if (exists $entry->{SWALL_AC}) {
      $uniprot_sth->execute($entry->{SWALL_AC});
      while (my ($acc, $type, $subtype, $de_text) = $uniprot_sth->fetchrow_array) {
        if ($subtype eq 'EC') {
          $entry->{EC} = $de_text;
        } elsif ($type eq 'RecName' and $subtype eq 'Full') {
          $entry->{DE} = $de_text if $de_text !~ /^Uncharacterized protein/;
        }
      }
    }
  }
}
die $uniprot_sth->errstr if $uniprot_sth->err;
$uniprot_sth->finish;
$uniprot_dbh->dbc->disconnect_if_idle;

##################
# output results
##################

@resultsArr = sort { $a->{NT_AC} cmp $b->{NT_AC} or $a->{AA_PID} cmp $b->{AA_PID} } map { @$_ } values(%resultsHash);
foreach my $entry (@resultsArr) {
  
  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
         $entry->{NT_AC},
         $entry->{NT_version},
         $entry->{AA_PID},
         $entry->{AA_version},
         $entry->{AA_checksum},
         $entry->{AA_text},
         exists($entry->{SWALL_AC}) ? $entry->{SWALL_AC} : ".",
         exists($entry->{SWALL_ID}) ? $entry->{SWALL_ID} : ".",
         exists($entry->{SWALL_ISOFORM}) ? $entry->{SWALL_ISOFORM} : ".",
         exists($entry->{EC}) ? $entry->{EC} : ".",
         exists($entry->{DE}) ? $entry->{DE} : ".",
      );

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


#####################
sub get_uniprot_dbh {
  my ($cred_file) = @_;

  my ($sid, $host, $port, $user, $pass);

  open(my $cfh, $cred_file) or die "Could not open $cred_file for reading\n";
  while(<$cfh>) {
    /^\#/ and next;

    /^DBNAME:(\S+)/ and $sid = $1;
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
    -dbname => $sid,
    -driver => 'Oracle');

  return $dbh;
}
