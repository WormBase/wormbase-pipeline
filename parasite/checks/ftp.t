use strict;
use warnings;
use Test::More;
use ProductionMysql;
use SpeciesFtp;
use File::Slurp;

my @our_species = sort ProductionMysql->staging->species("_core_$ENV{PARASITE_VERSION}_");
my @file_types = qw/annotations.gff3 canonical_geneset.gtf CDS_transcripts.fa genomic.fa genomic_masked.fa genomic_softmasked.fa mRNA_transcripts.fa protein.fa/;

ok(SpeciesFtp->current_staging->root_exists, "staging exists") or do {
  done_testing;
  exit;
};
for my $file_type (@file_types){
  subtest "staging $file_type" => sub {
    for my $species (@our_species){
      my $path = SpeciesFtp->current_staging->path_to($species, $file_type);
      ok (-s $path, $path);
    }
  };
}

ok(SpeciesFtp->dot_next->root_exists, ".next exists") or do {
  done_testing;
  exit;
};
for my $file_type (@file_types){
  subtest ".next $file_type" => sub {
    for my $species (@our_species){
      my $path = SpeciesFtp->dot_next->path_to($species, $file_type);
      ok (-s $path, $path);
    }
  };
}

unless( exists $ENV{SKIP_FTP_COMPARA_TESTS} && $ENV{SKIP_FTP_COMPARA_TESTS} ) {
   my @compara_species = grep {/.*_.*_.*/} map {chomp; $_} read_file("$ENV{PARASITE_CONF}/compara.species_list") or do {
     ok(0, "Can read $ENV{PARASITE_CONF}/compara.species_list");
     done_testing;
     exit;
   };
   for my $file_type (("orthologs.tsv", "paralogs.tsv")){
     subtest ".next $file_type" => sub {
       for my $species (@compara_species){
         my $path = SpeciesFtp->dot_next->path_to($species, $file_type);
         ok (-s $path, $path);
       }
     }
   }
}

done_testing;
