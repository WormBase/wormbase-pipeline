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
my ($ebi_header_prefix, $species_string, @test_genes);

GetOptions( 
  'host=s'           => \$dbhost,
  'port=s'           => \$dbport,
  'user=s'           => \$dbuser,
  'pass=s'           => \$dbpass,
  'dbname=s'         => \$database,
  'mrna'             => \$mrna,
  'cds'              => \$cds, 
  'pep'              => \$pep,
  'ebiblastheader=s' => \$ebi_header_prefix,
  'testgene=s'       => \@test_genes,
  'outfile=s'        => \$outfile) or die(@!);

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

if (defined $ebi_header_prefix) {
  $species_string = ucfirst( $db->get_MetaContainer->get_production_name );
}

my $seqio = (defined $outfile) 
    ? Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile")
    : Bio::SeqIO->new(-format => 'fasta', -fh => \*STDOUT);

my $gene_adaptor = $db->get_GeneAdaptor();
my @genes;

if (@test_genes) {
  map { push @genes, $gene_adaptor->fetch_by_stable_id($_) } @test_genes;
} else {
  @genes = @{$gene_adaptor->fetch_all()};
}

foreach my $gene(@genes){
  my $gene_id = $gene->stable_id();
  foreach my $trans (@{$gene->get_all_Transcripts()}) {
    my $trans_id = $trans->stable_id;
    my $slice_id = $trans->feature_Slice->name;
    my $status = lc($trans->status);

    my ($id, $desc_text, $seqstring);
    
    if ($pep) {
      next unless $trans->biotype eq 'protein_coding';

      my $tr = $trans->translation;
      
      my $pep_desc;
      if (defined $trans->description && $trans->description ne '' && $trans->description ne $trans->stable_id) {
        $pep_desc = $trans->description;
      } elsif (defined $gene) {
        $pep_desc = $gene->description;
      }
      if($pep_desc) {
        $pep_desc = $1 if $pep_desc =~ /(.*).\s\[.*\]/;
        $pep_desc =~ s/\"/\\"/g;
        $pep_desc = sprintf('"%s"', $pep_desc);
      }
      
      if (defined $ebi_header_prefix) {
        $id = join(":", $ebi_header_prefix, $species_string, $tr->stable_id);
        $desc_text = sprintf("pep:%s %s gene:%s transcript:%s species:%s", 
                             $status, $slice_id, $gene_id, $trans_id, $species_string);
        $desc_text .= "description:$pep_desc" if $pep_desc;
      } else {
        $id = $tr->stable_id;
        $desc_text = "transcript=$trans_id gene=$gene_id";
      }
      $seqstring = $trans->translation()->seq;
    } elsif ($cds) {
      next unless $trans->biotype eq 'protein_coding';

      if (defined $ebi_header_prefix) {
        $id = join(":", $ebi_header_prefix, $species_string, $trans_id);
        $desc_text = sprintf("cds:%s %s gene:%s transcript:%s species:%s", 
                             $status, $slice_id, $gene_id, $trans_id, $species_string);
      } else {
        $id = $trans_id;
        $desc_text = "gene=$gene_id";
      }      
      $seqstring = $trans->translateable_seq();
    } elsif ($mrna) {
      if (defined $ebi_header_prefix) {
        $id = join(":", $ebi_header_prefix, $species_string, $trans_id);
        $desc_text = sprintf("cdna:%s %s gene:%s transcript:%s species:%s", 
                             $status, $slice_id, $gene_id, $trans_id, $species_string);
      } else {
        $id = $trans_id;
        $desc_text = "gene=$gene_id";
      }

      $seqstring = $trans->seq->seq();
    }

    my $seq = Bio::PrimarySeq->new(-id => $id,
                                   -seq => $seqstring,
                                   -desc => $desc_text);

    $seqio->write_seq($seq);

  }
}

$seqio->close();

exit(0);
