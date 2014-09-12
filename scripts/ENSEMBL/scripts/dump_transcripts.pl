#!/usr/bin/perl -w

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use IO::File;
use Bio::SeqIO;

my $dbhost = "mysql-wormbase-pipelines";
my $dbuser = "wormro";
my $dbport = "4331";
my $dbpass = "";
    
my ($database,$outfile);

my ($cds, $pep, $mrna) = (0,0,0);

GetOptions( 
  'host=s'        => \$dbhost,
  'port=s'        => \$dbport,
  'user=s'        => \$dbuser,
  'pass=s'        => \$dbpass,
  'dbname=s'      => \$database,
  'mrna'          => \$mrna,
  'cds'           => \$cds, 
  'pep'           => \$pep,
  'outfile=s'     => \$outfile) or die(@!);

if (not $mrna and not $cds and not $pep or
    $cds + $mrna + $pep > 1) {
  die "You must supply exactly one of -pep -cds or -mrna\n";
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host   => $dbhost,
        -user   => $dbuser,
        -dbname => $database,
        -pass   => $dbpass,
        -port   => $dbport,
    );

my $seqio = (defined $outfile) 
    ? Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile")
    : Bio::SeqIO->new(-format => 'fasta', -fh => \*STDOUT);

my $gene_adaptor = $db->get_GeneAdaptor();
my @genes = @{$gene_adaptor->fetch_all()};
foreach my $gene(@genes){
  my $gene_id = $gene->stable_id();
  foreach my $trans (@{$gene->get_all_Transcripts()}) {
    my $trans_id = $trans->stable_id;

    my ($id, $desc_text, $seqstring);
    
    if ($pep) {
      next unless $trans->biotype eq 'protein_coding';

      my $tr = $trans->translation;
      $id = $tr->stable_id;
      $desc_text = "transcript_id=$trans_id gene_id=$gene_id";
      $seqstring = $trans->translation()->seq;
    } elsif ($cds) {
      next unless $trans->biotype eq 'protein_coding';

      $id = $trans_id;
      $seqstring = $trans->translateable_seq();
      $desc_text = "gene_id=$gene_id";
    } elsif ($mrna) {
      $id = $trans_id;
      $seqstring = $trans->seq->seq();
      $desc_text = "gene_id=$gene_id";
    }

    my $seq = Bio::PrimarySeq->new(-id => $id,
                                   -seq => $seqstring,
                                   -desc => $desc_text);

    $seqio->write_seq($seq);

  }
}

$seqio->close();

exit(0);
