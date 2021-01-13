#!/usr/bin/env perl
#
use strict;
use warnings;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Getopt::Long;
use Storable;

use Wormbase;
use Log_files;
use WormBase2Ensembl;
use Bio::SeqIO;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Utils::IO::GTFSerializer;

my ($debug, $test, $store, $species,  $verbose, $wb, $gene_name,
    $out_file, $out_fh, $gff3, $genome, %genes_by_slice);


&GetOptions(
  "debug=s"     => \$debug,
  "test"        => \$test,
  "verbose"     => \$verbose,
  "store:s"     => \$store,
  "species:s"   => \$species,
  "genome=s"    => \$genome,
  "gff3=s"      => \$gff3,
  "outgtf=s"    => \$out_file,
  'public_name' => \$gene_name,
 )or die ("Couldn't get options");


if ( $store ) {
  $wb = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wb = Wormbase->new( -debug    => $debug,
                       -test     => $test,
                       -organism => $species,
      );
}

my $log = Log_files->make_build_log($wb);


$genome = $wb->genome_seq if not defined $genome;
$gff3 = $wb->processed_GFF_file( 1 ) if not defined $gff3;
if (not -e $gff3) {
  $gff3 .= ".gz";
  if (not -e $gff3) {
    $log->log_and_die("Could not find GFF3 file $gff3\n");
  }
}
$out_file =  sprintf("%s/%s.gtf", $wb->sequences, $wb->species) if not defined $out_file;

my $slices = &make_slices($genome);

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
  open( $out_fh, ">$out_file") or $log->log_and_die("Could not open $out_file for writing\n");
} else {
  $log->log_and_die("$out_file not specified\n");
}

print $out_fh "#!genebuild-version ", $wb->get_wormbase_version_name, "\n";

my $serializer = Bio::EnsEMBL::Utils::IO::GTFSerializer->new($out_fh);
foreach my $slice (values %$slices) {
  if (exists $genes_by_slice{$slice->seq_region_name}) {
    foreach my $g (sort { $a->start <=> $b->start } @{$genes_by_slice{$slice->seq_region_name}}) {
      $serializer->print_Gene($g);
    }
  }
}

$serializer = undef;
close $out_fh;

add_public_name($out_file) if $gene_name;


$log->mail();
exit(0);

###################################
# add public_name to GTF as gene_name 
sub add_public_name{
    my ($out) = @_;
    
    # open acedb connection
    my $db=Ace->connect(-path => $wb->autoace)||die($Ace::Error);
    # open temp file handle
    my $tmpfile='/tmp/gtf.tmp';
    open OUTF, ">$tmpfile";
    
    # read file
    open INF, $out;
    while (<INF>){
      # parse WBGeneID
      if (/WormBase\s+gene\s.*(WBGene\d+)/){
         chomp;
         my $gene = $db->fetch(Gene => "$1")||$log->log_and_die("cannot find gene: $1\n");
         # look up the Public_name
         my $public_name = $gene->Public_name;
         print OUTF "$_ gene_name \"$public_name\";\n";
      }else{
         print OUTF;
      }
    }
    close OUTF;
    close INF;

    # move the temp file to the correct location
    `mv -f $tmpfile $out`;
    
}



#######################################
sub make_slices {
  my ($genome) = @_;

  my $gfh;
  if ($genome =~ /\.gz$/) {
    open($gfh, "gunzip -c $genome |") or $log->log_and_die("Could not open gunzip stream to $genome\n");
  } else {
    open($gfh, $genome) or $log->log_and_die("Could not open $genome for reading\n");
  }
    
  my $seqio = Bio::SeqIO->new(-fh => $gfh,
                              -format => 'fasta');

  my %slices;

  while(my $seq = $seqio->next_seq) {
    my $name = $seq->id;
    if ($wb->species eq 'elegans' and $name =~ /^CHROMOSOME_/) {
      $name =~ s/^CHROMOSOME_//; 
    } elsif ($wb->species eq 'briggsae' and $name =~ /^chr/) {
      $name =~ s/^chr//; 
    }
    my $slice = Bio::EnsEMBL::Slice->new(
      -seq               => uc($seq->seq),
      -seq_region_name   => $name,
      -start             => 1,
      -end               => $seq->length,
      -seq_region_length => $seq->length,
      -strand            => 1);
    $slices{$name} = $slice;
  }


  return \%slices;
}
