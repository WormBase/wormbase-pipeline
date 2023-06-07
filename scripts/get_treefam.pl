#!/usr/bin/env perl
# get a TreeFam family ids from the TreeFam compara database
# needs to be run at the EBI
#
# Version: TreeFam 9 / 69

use DBI;
use Wormbase;
use Log_files;
use Getopt::Long;
use Storable;
use IO::File;
use strict;

my ($debug,$store,$species,$test,$acefile,$no_load);
GetOptions(
  'debug=s'   => \$debug,
  'store=s'   => \$store,
  'species=s' => \$species,
  'test'      => \$test,
  'acefile=s' => \$acefile,
  'noload'    => \$no_load,
)||die($@);


my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species
      );
}
my $log = Log_files->make_build_log($wormbase);

$acefile ||= $wormbase->acefiles."/treefam.ace";
my $out = IO::File->new(">$acefile") or $log->log_and_die("Cannot write to $acefile\n");

my $dsn = 'DBI:mysql:database=treefam_production_9_69;host=mysql-ens-compara-prod-4.ebi.ac.uk;user=ensro;port=4401';

$log->write_to("connecting to: $dsn\n");
my $dbh = DBI->connect($dsn,'treefam_ro')||$log->log_and_die($@);

my $sth=$dbh->prepare('SELECT t.stable_id FROM member m JOIN gene_tree_member USING(member_id) JOIN gene_tree_node USING(node_id) JOIN gene_tree_root t USING(root_id) WHERE source_name="ENSEMBLPEP" AND m.stable_id=? AND t.stable_id IS NOT NULL;')||$log->log_and_die($@);
my $sth2=$dbh->prepare('SELECT t.stable_id FROM member m JOIN gene_tree_member USING(member_id) JOIN gene_tree_node USING(node_id) JOIN gene_tree_root t USING(root_id) WHERE source_name="ENSEMBLPEP" AND m.stable_id like ? AND t.stable_id IS NOT NULL;')||$log->log_and_die($@);


# actual bit that does things

my %cds2wormpep = $wormbase->FetchData('cds2wormpep');

while (my($cds,$protein_id)=each %cds2wormpep){
 my @ids;
 $sth->execute($cds);
 while (my @id = $sth->fetchrow_array()) {
           push @ids,"$id[0]";
 }
 unless ($ids[0]){
   $sth2->execute("$cds.%");
   while (my @id = $sth2->fetchrow_array()) {
           push @ids,"$id[0]";
   }
 }

 if ($ids[0]){
  print $out "Protein : $protein_id\n";
  map {print $out "Database TREEFAM TREEFAM_ID $_\n"} @ids;
  print $out "\n";
  print $out "CDS : $cds\n";
  map {print $out "Database TREEFAM TREEFAM_ID $_\n"} @ids;
  print $out "\n";
 }
}

$out->close;

$wormbase->load_to_database($wormbase->autoace, $acefile, 'treefam', $log) unless ($no_load or $test);

$log->mail;
