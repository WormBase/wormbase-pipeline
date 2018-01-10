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
    not defined $uniprot_cred or
    not defined $org_id) {
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
# statusid = 4 => public/active records only (i.e. not suppressed)
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

my (%resultsHash, %uniparc_mapping, $gcrp_release, %gcrp_members);

while ( ( my @results ) = $ena_sth->fetchrow_array ) {
  $resultsHash{$results[2]} =  {  
    NT_AC       => $results[0],
    NT_version  => $results[1],
    AA_PID      => $results[2],
    AA_version  => $results[3],
    AA_checksum => $results[4],
    AA_stdname  => $results[5],
    AA_ID       => $results[6],
  };
}
die $ena_sth->errstr if $ena_sth->err;
$ena_sth->finish;
$ena_dbh->dbc->disconnect_if_idle;


my ($uniprot_dbh, $uniprot_sql, $uniprot_sth);

#################
# query uniprot 1 - get reference proteome information
#################
$uniprot_dbh = &get_uniprot_dbh($uniprot_cred);

$uniprot_sql = "SELECT MAX(release) from sptr.gene_centric_entry WHERE tax_id = $org_id";
$uniprot_sth = $uniprot_dbh->dbc->prepare($uniprot_sql);
$uniprot_sth->execute();

($gcrp_release) = $uniprot_sth->fetchrow_array;

die $uniprot_sth->errstr if $uniprot_sth->err;
$uniprot_sth->finish;

$uniprot_sql = "SELECT DISTINCT(accession) from sptr.gene_centric_entry WHERE tax_id = $org_id and is_canonical = 1 AND release = '$gcrp_release'";
$uniprot_sth = $uniprot_dbh->dbc->prepare($uniprot_sql);
$uniprot_sth->execute();
while( (my @results) = $uniprot_sth->fetchrow_array) {
  my ($acc) = @results;
  $gcrp_members{$acc} = 1;
}

die $uniprot_sth->errstr if $uniprot_sth->err;
$uniprot_sth->finish;
$uniprot_dbh->dbc->disconnect_if_idle;

#################
# query uniprot 2 - get basic protein_id info
#################

$uniprot_dbh = &get_uniprot_dbh($uniprot_cred);

$uniprot_sql =  "SELECT e.accession, et.descr, p.protein_id "
    . "FROM sptr.dbentry e, sptr.embl_protein_id p, sptr.cv_entry_type et "
    . "WHERE p.dbentry_id = e.dbentry_id "
    . "AND e.deleted='N' "
    . "AND e.merge_status <> 'R' "
    . "AND e.entry_type = et.entry_type_id " 
    . "AND et.descr in ('Swiss-Prot', 'TrEMBL') " 
    . "AND e.tax_id = $org_id";


$uniprot_sth = $uniprot_dbh->dbc->prepare($uniprot_sql);
print STDERR "Reading Uniprot database to get accessions and ids...\n" if $verbose;
$uniprot_sth->execute();

while( (my @results) = $uniprot_sth->fetchrow_array) {
  my ($uniacc,  $unitype, $pid) = @results;

  if (exists $resultsHash{$pid}) {
    my $el = $resultsHash{$pid};
    $el->{SWALL_AC}   = $uniacc;
    $el->{SWALL_TYPE} = ($unitype eq 'Swiss-Prot') ? "SP" : "TR";
    $el->{SWALL_GCRP} = $gcrp_release if exists $gcrp_members{$uniacc};
  } elsif (exists $gcrp_members{$uniacc}) {
    # these are entries in the GCRP, but not associated with a WormBase genome project annotation
    # (either because they are on the Mitochondrion, or because they are independently curated
    # by SwissProt. We will add these to the file, but without an associated parent sequence
    $resultsHash{$pid} = {
      AA_PID       => $pid,
      SWALL_AC     => $uniacc,
      SWALL_TYPE   => ($unitype eq 'Swiss-Prot') ? "SP" : "TR",
      SWALL_GCRP   => $gcrp_release, 
    };
  }
}

die $uniprot_sth->errstr if $uniprot_sth->err;
$uniprot_sth->finish;
$uniprot_dbh->dbc->disconnect_if_idle;

##################
# query uniprot 3 - use UniParc to get Uniprot isoform identifiers
###################

$uniprot_dbh = &get_uniprot_dbh($uniprot_cred);

# dbid = 1 => ENA protein_id
# dbid = 58 => ENA protein_id on CON records
# dbid = 24 => Uniprot isoform id
$uniprot_sql = 'SELECT e.UPI, e.AC, e.DBID '
    . 'FROM uniparc.xref@uapro e '
    . "WHERE taxid = $org_id "
    . "AND  deleted = 'N' "
    . 'AND  dbid in (1,24, 58)';

$uniprot_sth = $uniprot_dbh->dbc->prepare($uniprot_sql);
print STDERR "Reading UniParc for mapping of ENA protein ids to UP isoform ids\n" if $verbose;
$uniprot_sth->execute;
while (my @results = $uniprot_sth->fetchrow_array) {
  my ($upi, $acc, $db_id) = @results;
  if ($db_id == 1 or $db_id == 58) {
    # ena protein id; ignore if it's not one of ours
    if (exists $resultsHash{$acc}) {
      $uniparc_mapping{$upi}->{protein_id} = $acc;
    } else {
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
    if (exists $resultsHash{$pid}) {
      $resultsHash{$pid}->{SWALL_ISOFORM} = $iid;
    }
  }
}

##################
# query uniprot 4 - get EC numbers
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
foreach my $entry (values %resultsHash) {
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
die $uniprot_sth->errstr if $uniprot_sth->err;
$uniprot_sth->finish;
$uniprot_dbh->dbc->disconnect_if_idle;

##################
# Putting it all together
##################

my @WBresults = grep { exists $_->{AA_stdname} } values(%resultsHash);
my @nonWBresults =  grep { not exists $_->{AA_stdname} } values(%resultsHash);
my %wb_swall;

foreach my $entry (sort { $a->{NT_AC} cmp $b->{NT_AC} or $a->{AA_PID} cmp $b->{AA_PID} } @WBresults ) {
  
  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
         exists($entry->{NT_AC}) ? $entry->{NT_AC} : ".",
         exists($entry->{NT_version}) ? $entry->{NT_version} : ".",
         exists($entry->{AA_PID}) ? $entry->{AA_PID} : ".",
         exists($entry->{AA_version}) ? $entry->{AA_version} : ".",
         exists($entry->{AA_checksum}) ? $entry->{AA_checksum} : ".",
         exists($entry->{AA_stdname}) ? $entry->{AA_stdname} : ".",
         exists($entry->{SWALL_AC}) ? $entry->{SWALL_TYPE} . ":" . $entry->{SWALL_AC} : ".",
         exists($entry->{SWALL_ISOFORM}) ? $entry->{SWALL_ISOFORM} : ".",
         exists($entry->{EC}) ? $entry->{EC} : ".",
         exists($entry->{DE}) ? $entry->{DE} : ".",
         exists($entry->{SWALL_GCRP}) ? $entry->{SWALL_GCRP} : ".",
      );
  $wb_swall{$entry->{SWALL_AC}} = 1 if exists $entry->{SWALL_AC};

}
foreach my $entry (sort { $a->{SWALL_AC} cmp $b->{SWALL_AC} or  $a->{AA_PID} cmp $b->{AA_PID} } @nonWBresults ) {
  next if exists $wb_swall{$entry->{SWALL_AC}};

  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
         ".",
         ".",
         $entry->{AA_PID}, 
         ".",
         ".",
         ".",
         exists($entry->{SWALL_AC}) ? $entry->{SWALL_TYPE} . ":" . $entry->{SWALL_AC} : ".",
         exists($entry->{SWALL_ISOFORM}) ? $entry->{SWALL_ISOFORM} : ".",
         exists($entry->{EC}) ? $entry->{EC} : ".",
         exists($entry->{DE}) ? $entry->{DE} : ".",
         exists($entry->{SWALL_GCRP}) ? $entry->{SWALL_GCRP} : ".",
      );
  $wb_swall{$entry->{SWALL_AC}} = 1;
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
