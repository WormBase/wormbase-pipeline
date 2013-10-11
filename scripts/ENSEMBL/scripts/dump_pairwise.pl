#!/usr/bin/env perl
# very unsophisticated script to dump the 6 pairwise alignment pairs from EG


use strict;
use Bio::EnsEMBL::Registry; # the registry is release dependent and overriden to 73 below
use Bio::EnsEMBL::Utils::Exception qw(throw);

my $method = 'LASTZ_NET';
Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'mysql.ebi.ac.uk', -user => 'anonymous',-port => 4157,-db_version => 73);

# that below is release dependent
my $comparaDB = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
        -user   =>'anonymous',
        -port   => 4157,
        -dbname => 'ensembl_compara_metazoa_20_73',
        -host   => 'mysql.ebi.ac.uk',
);

my $genome_db_adaptor = $comparaDB->get_GenomeDBAdaptor();
my $mlss_adaptor      = $comparaDB->get_MethodLinkSpeciesSetAdaptor();
my $block_adaptor     = $comparaDB->get_GenomicAlignBlockAdaptor();

my @species = ('caenorhabditis_japonica','caenorhabditis_briggsae','caenorhabditis_brenneri','caenorhabditis_remanei','brugia_malayi');

my $outf = IO::File->new('>c_elegans.genomic_alignment.gff3')||die(@!);
map {dump_me('caenorhabditis_elegans',$_,$outf)} @species;

foreach my $s(@species){
   my $sname = "$1_$2" if  $s=~/^(\w)\w+_(\w+)/;
   my $of = IO::File->new(">$sname.genomic_alignment.gff3") ||die(@!);
   dump_me($s,'caenorhabditis_elegans',$of);
}

sub dump_me {
 my ($species1,$species2,$file)=@_;

 # that is to get the short form name that can be picked up by the build 
 $species2=~/^(\w)\w+(_\w+)/;
 my $short_name = "$1$2";

 my $mlss = $mlss_adaptor->fetch_by_method_link_type_registry_aliases($method,[$species1,$species2]);
 my $slice_adaptor     = Bio::EnsEMBL::Registry->get_adaptor($species1,'core','Slice');

 my $slices = $slice_adaptor->fetch_all('toplevel');
 while(my $slice = shift @$slices){
          my $blocks = $block_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss,$slice);
          while (my $block = shift @$blocks){
             my $aligns = $block->get_all_non_reference_genomic_aligns();
             while (my $frag = shift @$aligns){
               printf $file "%s\t%s\t%s\t%s\t%s\t%s\t%s\t\%s\tTarget=%s %s %s %s\n",
                       $block->reference_slice->seq_region_name,
                       $short_name.'_'.$method,
                       'conserved_region',
                       $block->reference_slice_start,
                       $block->reference_slice_end,
                       '.', #score
                       ($block->reference_slice_strand > 0 ? '+':'-'),
                       '.', # phase
                       $short_name.':'.$frag->dnafrag->name,
                       $frag->dnafrag_start,
                       $frag->dnafrag_end,
                       ($frag->dnafrag_strand>0?'+':'-');
             }
          }
 }
}
