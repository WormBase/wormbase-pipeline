use strict;
use feature 'say';
use ProductionMysql;
use List::MoreUtils qw/first_index last_index/;
my %chosen_bioprojects = (
    'ancylostoma_ceylanicum' => 'prjna231479',
    'ascaris_suum' => 'prjna62057',
    'caenorhabditis_remanei' => 'prjna53967',
    'dictyocaulus_viviparus' => 'prjna72587',
    'echinococcus_granulosus' => 'prjeb121',
    'fasciola_hepatica' => 'prjeb25283',
    'haemonchus_contortus' => 'prjeb506',
    'heligmosomoides_polygyrus' => 'prjeb1203',
    'loa_loa' => 'prjna246086',
    'macrostomum_lignano' => 'prjna371498',
    'meloidogyne_arenaria' => 'prjeb8714',
    'onchocerca_flexuosa' => 'prjna230512',
    'onchocerca_ochengi' => 'prjeb1204',
    'taenia_asiatica' => 'prjna299871',
    'toxocara_canis' => 'prjeb533',
    'trichinella_nativa' => 'prjna179527',
    'trichinella_pseudospiralis' => 'iss13prjna257433',
    'trichinella_spiralis' => 'prjna12603',
    'trichuris_suis' => 'prjna179528',
    'wuchereria_bancrofti' => 'prjeb536',
);

my %species;
for my $species (ProductionMysql->staging->species){
  my ($spe, $cies, $bp) = split "_", $species;
  next unless $bp;
  $species{"${spe}_${cies}"} //= [];
  push $species{"${spe}_${cies}"}, $bp;
}
my @species;
for my $sp (sort keys %species) {
   my @bps = @{$species{$sp}};
   if (@bps > 1 ){
     die "$sp pick BioProject:". join (", ", @bps) unless grep({$_ eq $chosen_bioprojects{$sp} } @bps);
     die unless "$chosen_bioprojects{$sp}";
     push @species, join("_",$sp, $chosen_bioprojects{$sp});
   } else {
     die unless $bps[0];
     push @species, join("_",$sp, $bps[0]);
   }
}
my $ftp_base = "ftp://ftp.ebi.ac.uk/pub/databases/wormbase/collaboration/EBI/GxA/latest";
sub line {
   my %args = @_;
   my ($spe, $cies, $bp) = split "_", $args{core_db};
   $bp =~ s/iss13prjna257433/iss13_prjna257433/;
   $bp = uc($bp);
   return join " ", (
     "${spe}_${cies}",
     $args{taxon},
     "wbps",
     sprintf("$ftp_base/%s.%s.WBPSRELNO.genomic.fa.gz", "${spe}_${cies}", $bp),
     sprintf("$ftp_base/%s.%s.WBPSRELNO.mRNA_transcripts.fa.gz", "${spe}_${cies}", $bp),
     sprintf("$ftp_base/%s.%s.WBPSRELNO.canonical_geneset.gtf.gz", "${spe}_${cies}", $bp),
     $args{assembly}
   );
}

my $inh = *STDIN;
open ($inh, "<", @ARGV[0]) if @ARGV;
my @lines = (<$inh>);
die unless @lines;
my $first_index = first_index { $_ =~ /wbps/} @lines;
my $last_index = last_index { $_ =~ /wbps/} @lines;
my $outh = *STDOUT;
open ($outh, ">" , @ARGV[0]) if @ARGV;

map {print $outh $_} @lines[0 .. $first_index - 1 ];

for my $species (@species){
  my $core_db = ProductionMysql->staging->core_db($species);
  say $outh &line(
    core_db => $core_db,
    taxon => ProductionMysql->staging->meta_value($core_db, "species.taxonomy_id"),
    assembly => ProductionMysql->staging->meta_value($core_db, "assembly.name"),
  );
}
map {print $outh $_} @lines[$last_index +1 .. @lines - 1 ];
