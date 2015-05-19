#!/usr/bin/perl -w

use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBEntry;

my (
  $dbname,
  $dbhost,
  $dbuser,
  $dbport,
  $dbpass,
  $xref_file,
  $no_write,
);

$dbuser = 'ensro';

&GetOptions(
  'dbname=s' => \$dbname,
  'user=s'   => \$dbuser,
  'host=s'   => \$dbhost,
  'port=s'   => \$dbport,
  'pass=s'   => \$dbpass,
  'xref=s'   => \$xref_file,
  'nowrite'  => \$no_write,
);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	'-dbname' => $dbname,
	'-host' => $dbhost,
	'-user' => $dbuser,
	'-port' => $dbport,
	'-pass' => $dbpass
  
);

my (%gene2gseq, %gene2wbgene, %gene2locus, %transcript2tran, %translation2pep);

my $xref_fh;
if ($xref_file =~ /\.gz$/) {
  open( $xref_fh, "gunzip -c $xref_file |") or die "Could not open gzip stream on $xref_file\n";
} else {
  open( $xref_fh, $xref_file) or die "Could not open $xref_file for reading\n";
}

while(<$xref_fh>) {
  my ($gseq_id, $wbgene_id, $locus_id, $transcript_id, $wormpep_id) = split(/\s+/, $_);

  $gene2wbgene{$wbgene_id} = Bio::EnsEMBL::DBEntry->new(
    -primary_id  => $wbgene_id,
    -display_id  => $wbgene_id,
    -dbname      => 'wormbase_gene',
    -version     => 0,
    -info_type   => 'DIRECT',
    -description => "",
    -info_text   => "");
  
  if ($gseq_id ne '.') {
    $gene2gseq{$wbgene_id} = Bio::EnsEMBL::DBEntry->new(
      -primary_id  => $wbgene_id,
      -display_id  => $gseq_id,
      -dbname      => 'wormbase_gseqname',
      -version     => 0,
      -info_type   => 'DIRECT',
      -description => "",
      -info_text   => "");      
  } 
  if ($locus_id ne '.') {
    $gene2locus{$wbgene_id} = Bio::EnsEMBL::DBEntry->new(
      -primary_id  => $wbgene_id,
      -display_id  => $locus_id,
      -dbname      => 'wormbase_locus',
      -version     => 0,
      -info_type   => 'DIRECT',
      -description => "",
      -info_text   => "");
    
  }
  if ($transcript_id ne '.') {
    $transcript2tran{$transcript_id} = Bio::EnsEMBL::DBEntry->new(
      -primary_id  => $transcript_id,
      -display_id  => $transcript_id,
      -dbname      => 'wormbase_transcript',
      -version     => 0,
      -info_type   => 'DIRECT',
      -description => "",
      -info_text   => "");
  }
  if ($wormpep_id ne '.') {
    $translation2pep{$transcript_id} = Bio::EnsEMBL::DBEntry->new(
      -primary_id  => $wormpep_id,
      -display_id  => $wormpep_id,
      -dbname      => 'wormpep_id',
      -version     => 0,
      -info_type   => 'DIRECT',
      -description => "",
      -info_text   => "");
  }
}

my (@to_write);

foreach my $g (@{$db->get_GeneAdaptor->fetch_all}) {
  if (exists $gene2wbgene{$g->stable_id}) {
    push @to_write, [$gene2wbgene{$g->stable_id}, $g->dbID, 'Gene'];
  }
  if (exists $gene2gseq{$g->stable_id}) {
    push @to_write, [$gene2gseq{$g->stable_id}, $g->dbID, 'Gene'];
  }
  if (exists $gene2locus{$g->stable_id}) {
    push @to_write, [$gene2locus{$g->stable_id}, $g->dbID, 'Gene'];
  }
}


foreach my $t (@{$db->get_TranscriptAdaptor->fetch_all}) {
  if (exists $transcript2tran{$t->stable_id}) {
    push @to_write, [$transcript2tran{$t->stable_id}, $t->dbID, 'Transcript'];
  }
}

foreach my $tr (@{$db->get_TranslationAdaptor->fetch_all}) {
  if (exists $translation2pep{$tr->stable_id}) {
    push @to_write, [$translation2pep{$tr->stable_id}, $tr->dbID, 'Translation'];
  }
}


if ($no_write) {
  foreach my $obj (@to_write) {
    my ($dbentry, $ensid, $enstype) = @$obj;
    printf "Would write %s %s %d %s\n", $dbentry->primary_id, $dbentry->display_id, $ensid, $enstype;
  }
} else {
  foreach my $dbname ("wormbase_gene", "wormbase_gseqname", "wormbase_locus", "wormbase_transcript", "wormpep_id") {
    my $del_sql = "DELETE object_xref.*, xref.* "
        . "FROM external_db, xref, object_xref "
        . "WHERE external_db.external_db_id = xref.external_db_id "
        . "AND xref.xref_id = object_xref.xref_id "
        . "AND db_name = '$dbname'";

    my $sth = $db->dbc->prepare($del_sql);
    $sth->execute();
    $sth->finish;
  }

  foreach my $obj (@to_write) {
    my ($dbentry, $ensid, $enstype) = @$obj;
    $db->get_DBEntryAdaptor->store($dbentry, $ensid, $enstype, 1);
  }
}
