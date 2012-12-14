#!/usr/bin/env perl

use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my (
  $dbname,
  $dbhost,
  $dbuser,
  $dbport,
  $dbpass,
  @seq_region,
  $codon_table,
    );

&GetOptions(
  'dbname=s' => \$dbname,
  'dbuser=s' => \$dbuser,
  'dbhost=s' => \$dbhost,
  'dbport=s' => \$dbport,
  'dbpass=s' => \$dbpass,
  'codontable=s' => \$codon_table,
    );


my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  '-dbname' => $dbname,
  '-host' => $dbhost,
  '-user' => $dbuser,
  '-port' => $dbport,
  '-pass' => $dbpass
);


$codon_table = 5 if not defined $codon_table;
if (not @ARGV) {
  my $sl = $db->get_SliceAdaptor->fetch_by_region('toplevel',
                                                  'MtDNA');
  push @seq_region, $sl;
} else {

  foreach my $id (@ARGV) {
    my $sl = $db->get_SliceAdaptor->fetch_by_name($id);
    if (not defined $sl) {
      die "Could not fetch slice $id from db\n";
    }
    push @seq_region, $sl;
  }
}


my $attrib_type_id = &fetch_attrib_type_id_by_code('codon_table');

if (defined $attrib_type_id) {
  &store_on_slices($attrib_type_id, $codon_table, @seq_region);
}


###################################

sub fetch_attrib_type_id_by_code {
  my $code = shift;

  my $attrib_type_id;

  my $sth = $db->dbc->prepare("SELECT attrib_type_id from attrib_type where code = ?");
  $sth->execute($code);

  if($sth->rows()) {
    ($attrib_type_id) = $sth->fetchrow_array();
  }
  $sth->finish();


  return $attrib_type_id;
}



sub store_on_slices {
  my ($attrib_type_id, $value, @slices) = @_;

  my $sth = $db->dbc->prepare("INSERT into seq_region_attrib VALUES (?, $attrib_type_id, '$value')");
  foreach my $slice (@slices) {
    my $sr_id = $slice->get_seq_region_id;
    $sth->execute($sr_id);
  }
  $sth->finish;

}
