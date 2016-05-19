#!/usr/bin/env perl

use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SeqRegionSynonym;

my (
  $dbname,
  $dbhost,
  $dbuser,
  $dbport,
  $dbpass,
  $test
    );


&GetOptions(
  'dbname=s' => \$dbname,
  'user=s' => \$dbuser,
  'host=s' => \$dbhost,
  'port=s' => \$dbport,
  'pass=s' => \$dbpass,
  'test'   => \$test,
);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	'-dbname' => $dbname,
	'-host' => $dbhost,
	'-user' => $dbuser,
	'-port' => $dbport,
	'-pass' => $dbpass
);

my (%extdb_ids, @all_syns);

while(<>) {
  /^\#/ and next;
  my ($coord_sys, $seq_region_name, $synonym, $extdb) = split(/\s+/, $_);

  my $slice = $db->get_SliceAdaptor->fetch_by_region($coord_sys, 
                                                      $seq_region_name);

  if (not defined $slice) {
    die "Failed to fetch slice $coord_sys:$seq_region_name\n";
  }

  if (not exists $extdb_ids{$extdb}) {
    my $dbid = $db->get_DBEntryAdaptor->get_external_db_id($extdb, undef, 1);
    if (not defined $dbid) {
      die "Failed to find external_db entry for name '$extdb'\n";
    }
    $extdb_ids{$extdb} = $dbid;
  }

  my $sr_id = $slice->get_seq_region_id;
  my $edb_id = $extdb_ids{$extdb};

  my $sr_syn = Bio::EnsEMBL::SeqRegionSynonym->new(-seq_region_id => $sr_id,
                                                   -external_db_id => $edb_id,
                                                   -synonym => $synonym);

  push @all_syns, $sr_syn;
  
}

if ($test) {
  printf "Would have stored %d seq-region synonyms\n", scalar(@all_syns);
} else {
  foreach my $syn (@all_syns) {
    $db->get_SeqRegionSynonymAdaptor->store($syn);
  }
}
