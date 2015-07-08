#!/usr/bin/env perl 

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use IO::File;

my $dbhost = "mysql-wormbase-pipelines";
my $dbuser = "wormro";
my $dbport = "4331";
my $dbpass = "";
    
my ($database,$dna,$transcript_only,$outfile);

GetOptions( 
  'host=s'        => \$dbhost,
  'port=s'        => \$dbport,
  'user=s'        => \$dbuser,
  'pass=s'        => \$dbpass,
  'dbname=s'        => \$database,
  'dna'             => \$dna,
  'transcript_only' => \$transcript_only, 
  'outfile=s'       => \$outfile) or die(@!);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host   => $dbhost,
        -user   => $dbuser,
        -dbname => $database,
        -pass   => $dbpass,
        -port   => $dbport,
    );

my $fh;
if (defined $outfile) {
  open($fh, ">$outfile") or die "Could not open $outfile for writing\n";
} else {
  $fh = \*STDOUT;
}


my $gene_adaptor = $db->get_GeneAdaptor();
my @genes = @{$gene_adaptor->fetch_all()};
foreach my $gene(@genes){
  my $geneId=$gene->stable_id();
  foreach my $trans (@{$gene->get_all_Transcripts()}) {
    my $protein=$trans->translation();
    my $transcriptId=$trans->stable_id();
    my $proteinId=$protein?$protein->stable_id():$transcriptId;
    my $geneID=$gene->stable_id();
    my $seq;
    if($dna){
      next unless $trans->biotype eq 'protein_coding';
      $seq=$trans->translateable_seq();
    }elsif($transcript_only){
      $seq=$trans->seq->seq();
    }else{
      next unless $trans->biotype eq 'protein_coding';
      $seq=$protein->seq()
    }
    printf $fh ">$proteinId\t$transcriptId\t$geneId\n%s\n",reformat($seq);
  }
}

sub reformat {
  my ($seq)=@_;
  my @bases = split(//,$seq);
  my $count=1;
  my $reformatted;
  foreach my $b(@bases){
    $reformatted.=$b;
    $reformatted.="\n" if ($count++ % 60 == 0 && $count < scalar(@bases));
  }
  return $reformatted;
}
