#!/usr/bin/env perl -w

use DBI;
use strict;
use Getopt::Long;

my ($org_id, $ena_db, $spidx, $tridx, $verbose); 

&GetOptions(
            'swissidx=s'      => \$spidx,
            'tremblidx=s'     => \$tridx,
            'enadb=s'         => \$ena_db,
            'orgid=s'         => \$org_id,
            'verbose'         => \$verbose,
    );

if (not defined $spidx or
    not defined $tridx or
    not defined $org_id or
    not defined $ena_db) {
  die "Incorrect invocation\n";
}


my %attr   = ( PrintError => 0,
               RaiseError => 0,
               AutoCommit => 0 );

my $dbh = DBI->connect("dbi:Oracle:$ena_db", 'ena_reader', 'reader', \%attr)
    || die "Can't connect to database: $DBI::errstr";

# locus_tag = 84
# pseudo = 28
# standard name = 23
# gene = 12
# Get most info for each PID with a /locus_tag + featID was gene (#12)
my $sql1 =  "SELECT d.primaryacc#, d.version#, c.PROTEIN_ACC, c.version, c.chksum, fq.text, c.featid, d.project#, d.statusid"
    . " FROM cdsfeature c, dbentry d, feature_qualifiers fq"
    . " WHERE d.primaryacc# IN ("
    . "   SELECT primaryacc#"
    . "   FROM dbentry" 
    . "   JOIN sourcefeature USING (bioseqid)"
    . "   WHERE organism = $org_id"
    . "   AND project# = 1"
    . "   AND statusid = 4)"
    . " AND c.bioseqid = d.bioseqid"
    . " AND fq.featid  = c.featid"
    . " AND not exists (SELECT 1"
    . "   FROM feature_qualifiers a"
    . "   WHERE a.featid = c.featid"
    . "   AND a.fqualid = 28)"
    . " AND fq.fqualid   = 12";
my $sql2 = "SELECT fq.text FROM feature_qualifiers fq WHERE fq.featID = ? AND fq.fqualid = 23";
    
my $sth1 = $dbh->prepare($sql1);
my $sth2 = $dbh->prepare($sql2);

print STDERR "Doing primary lookup of CDS entries in ENA ORACLE database...\n" if $verbose;

$sth1->execute || die "Can't execute statement: $DBI::errstr";

my (@resultsArr, %resultsHash);

while ( ( my @results ) = $sth1->fetchrow_array ) {
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
die $sth1->errstr if $sth1->err;
$sth1->finish;

print STDERR "Reading Uniprot idx files to get accessions...\n" if $verbose;

foreach my $file ($spidx, $tridx) {
  open(my $fh, $file) or die "Could not open $file for reading\n";

  print STDERR " reading through file $file...\n" if $verbose;

  while(<$fh>) {
    /^(\S+)\.(\d+)\s+(\S+)/ and do {
      my ($pidacc, $pidver, $uniprot_acc) = ($1, $2, $3);

      if (exists $resultsHash{$pidacc}) {
        foreach my $el (@{$resultsHash{$pidacc}}) {
          $el->{SWALL_AC}->{$uniprot_acc} = 1;
        }
      }
    }
  }
}

@resultsArr = sort { $a->{NT_AC} cmp $b->{NT_AC} or $a->{AA_PID} cmp $b->{AA_PID} } map { @$_ } values(%resultsHash);
foreach my $entry (@resultsArr) {
  my @standard_names;
  
  $sth2->bind_param( 1, $entry->{AA_ID} );
  $sth2->execute || die "Can't execute statement: $DBI::errstr";
  while ( ( my @results ) = $sth2->fetchrow_array ) {
    push @standard_names, @results;
  }
  
  # not interested in things without a standard_name
  #next if not @standard_names;
  
  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
         $entry->{NT_AC},
         $entry->{NT_version},
         $entry->{AA_PID},
         $entry->{AA_version},
         $entry->{AA_checksum},
         $entry->{AA_text},
         exists($entry->{SWALL_AC}) ? join(";", sort keys %{$entry->{SWALL_AC}}) : "UNDEFINED",
         @standard_names ? join(";", @standard_names) : "UNDEFINED");
}

die $sth2->errstr if $sth2->err;
$sth2->finish;
$dbh->disconnect;

exit(0);

