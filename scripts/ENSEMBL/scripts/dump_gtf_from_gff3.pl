#!/usr/bin/env perl
#
use strict;
use warnings;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Getopt::Long;

use WormBase2Ensembl;
use Bio::SeqIO;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Utils::IO::GTFSerializer;

my ($debug, $wb_version, $genome_file,  
    %genes_by_slice, $out_file, $out_fh, $gff3);


GetOptions(
  'genome=s'     => \$genome_file,
  'gff3=s'       => \$gff3,
  'outfile:s'    => \$out_file,
  'wbversion:s'  => \$wb_version,
  'debug'        => \$debug,
          )or die ("Couldn't get options");

my $slices = &make_slices($genome_file);

if ($gff3) {
  my $wb2ens = WormBase2Ensembl->new(    
    -slices  => $slices,
    -debug   => ($debug) ? 1 : 0,
    -verbose => 1);
  
  my $anadummy =  Bio::EnsEMBL::Analysis->new();

  my $genes = $wb2ens->parse_genes_gff3( $gff3, 
                                         $anadummy,
                                         $anadummy,
                                         $anadummy,
                                         { WormBase => 1, WormBase_imported => 1 });
                
  while(my $g = shift @$genes) {
    push @{$genes_by_slice{$g->slice->seq_region_name}}, $g;
  }
}


if ($out_file) {
  open( $out_fh, ">$out_file") or die "Could not open $out_file for writing\n";
} else {
  $out_fh = \*STDOUT;
}

print $out_fh "#!genebuild-version $wb_version\n";

my $serializer = Bio::EnsEMBL::Utils::IO::GTFSerializer->new($out_fh);
foreach my $slice (values %$slices) {
  if (exists $genes_by_slice{$slice->seq_region_name}) {
    foreach my $g (sort { $a->start <=> $b->start } @{$genes_by_slice{$slice->seq_region_name}}) {
      $serializer->print_Gene($g);
    }
  }
}

exit(0);


#######################################
sub make_slices {
  my ($genome) = @_;

  my $gfh;
  if ($genome =~ /\.gz$/) {
    open($gfh, "gunzip -c $genome |") or die "Could not open gunzip stream to $genome\n";
  } else {
    open($gfh, $genome) or die "Could not open $genome for reading\n";
  }
    
  my $seqio = Bio::SeqIO->new(-fh => $gfh,
                              -format => 'fasta');

  my %slices;

  while(my $seq = $seqio->next_seq) {
    my $slice = Bio::EnsEMBL::Slice->new(
      -seq               => $seq->seq,
      -seq_region_name   => $seq->id, 
      -start             => 1,
      -end               => $seq->length,
      -seq_region_length => $seq->length,
      -strand            => 1);
    $slices{$seq->id} = $slice;
  }


  return \%slices;
}
