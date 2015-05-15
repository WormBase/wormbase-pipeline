#!/usr/bin/env perl -w

use DBI;
use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my ($bioproject_id, $ena_cred, $verbose); 

&GetOptions(
  'enacred=s'       => \$ena_cred,
  'bioprojectid=s'  => \$bioproject_id,
  'verbose'         => \$verbose,
    );


if (
    not defined $bioproject_id or
    not defined $ena_cred) {
  die "Incorrect invocation: you must supply -enacred and -orgid\n";
}


##############
# Query ENA
#############

my $ena_dbh = &get_ena_dbh($ena_cred);

# /locus_tag = 84
# /standard name (wormbase transcript/CDS/pseudogene name) = 23
# /gene = 12
my $ena_sql =  "SELECT d.primaryacc#, sf.featid, fq.fqualid, fq.text"
    . " FROM dbentry d, cv_dataclass dc,  bioseq b, seqfeature sf, feature_qualifiers fq"
    . " WHERE d.primaryacc# IN ("
    . "   SELECT primaryacc#"
    . "   FROM dbentry" 
    . "   WHERE study_id = '$bioproject_id'"
    . "   AND statusid = 4)"
    . " AND d.dataclass = dc.dataclass"
    . " AND d.bioseqid = b.seqid"
    . " AND b.seqid = sf.bioseqid"
    . " AND sf.featid  = fq.featid"
    . " AND fq.fqualid IN (23, 84, 12)";

    
my $ena_sth = $ena_dbh->dbc->prepare($ena_sql) or die "Can't prepare statement: $DBI::errstr";

print STDERR "Doing primary lookup of locus_tag entries in ENA ORACLE database...\n" if $verbose;

$ena_sth->execute or die "Can't execute statement: $DBI::errstr";

my (%feats, %g_data);

while ( my ($clone_acc, $feat_id, $qual_id, $qual_val ) = $ena_sth->fetchrow_array ) {
  $feats{$feat_id}->{$qual_id} = $qual_val;
  $feats{$feat_id}->{clone_acc} = $clone_acc;
}
die $ena_sth->errstr if $ena_sth->err;
$ena_sth->finish;
$ena_dbh->dbc->disconnect_if_idle;

foreach my $fid (keys %feats) {
  if (exists $feats{$fid}->{"23"}) {
    if (exists $feats{$fid}->{"84"}) {
      $g_data{$feats{$fid}->{"23"}}->{locus_tag} = $feats{$fid}->{"84"};
    } 
    if (exists $feats{$fid}->{"12"}) {
      $g_data{$feats{$fid}->{"23"}}->{gene} = $feats{$fid}->{"12"};
    }
    if (exists $feats{$fid}->{clone_acc}) {
      $g_data{$feats{$fid}->{"23"}}->{clone_acc} = $feats{$fid}->{clone_acc};
    }
  }
}

foreach my $k (keys %g_data) {
  print "$k";
  printf "\t%s", (exists $g_data{$k}->{clone_acc}) ?  $g_data{$k}->{clone_acc} : "."; 
  printf "\t%s", (exists $g_data{$k}->{gene}) ?  $g_data{$k}->{gene} : "."; 
  printf "\t%s", (exists $g_data{$k}->{locus_tag}) ?  $g_data{$k}->{locus_tag} : "."; 
  print "\n";
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
