#!/usr/local/ensembl/bin/perl

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;
use Getopt::Long;

my $usage = "
getsynteny.pl  --host ensembldb.ensembl.org
               --user anonymous
               --dbname ensembl_compara_41
               --chr_names \"22\"
               --species1 \"Homo sapiens\"
               [--assembly1 NCBI30]
               --species2 \"Mus musculus\"
               [--assembly2 MGSC3]
               [--method_link_type SYNTENY]

";

my $help = 0;
my $host = 'ia64b';
my $user = 'wormro';
my $pass;
my $dbname = 'worm_compara_lagan';

my $species1 = 'Caenorhabditis elegans';
my $species1_assembly;
my $species2 = 'Caenorhabditis briggsae';
my $species2_assembly;
my $species3 = 'Caenorhabditis remanei';
my $species3_assembly;
my $species4 = 'Caenorhabditis brenneri';
my $species4_assembly;


my $method_link_type = "SYNTENY";

my $chr_names = "all";

&GetOptions('help' => \$help,
            'host:s' => \$host,
            'user:s' => \$user,
            'dbname:s' => \$dbname,
            'pass:s' => \$pass,
            'species1:s' => \$species1,
            'assembly1:s' => \$species1_assembly,
            'species2:s' => \$species2,
            'assembly2:s' => \$species2_assembly,
            'species3:s' => \$species3,
            'assembly3:s' => \$species3_assembly,
            'species4:s' => \$species4,
            'assembly4:s' => \$species4_assembly,
            'chr_names=s' => \$chr_names,
            'method_link_type=s' => \$method_link_type);

if ($help) {
  print $usage;
  exit 0;
}

my $dba = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor (-host        => $host,
                                                       -user        => $user,
                                                       -pass        => $pass,
                                                       -dbname      => $dbname);

my $gdba = $dba->get_GenomeDBAdaptor;
my $dfa = $dba->get_DnaFragAdaptor;
my $mlssa = $dba->get_MethodLinkSpeciesSetAdaptor;
my $sra = $dba->get_SyntenyRegionAdaptor;

my $gdb1 = $gdba->fetch_by_name_assembly($species1,$species1_assembly);
my $gdb2 = $gdba->fetch_by_name_assembly($species2,$species2_assembly);
my $gdb3 = $gdba->fetch_by_name_assembly($species3,$species3_assembly);
my $gdb4 = $gdba->fetch_by_name_assembly($species4,$species4_assembly);

my $mlss = $mlssa->fetch_by_method_link_type_GenomeDBs($method_link_type, [$gdb1, $gdb2,$gdb3,$gdb4]);

my $dfgs;

if (defined $chr_names and $chr_names ne "all") {
  my @chr_names = split /,/, $chr_names;
  foreach my $chr_name (@chr_names) {
    push @{$dfgs}, $dfa->fetch_by_GenomeDB_and_name($gdb1, $chr_name);
  }
} else {
  $dfgs = $dfa->fetch_all_by_GenomeDB_region($gdb1);
}

my $total_nb_syntenies = 0;
foreach my $df (@{$dfgs}) {
  my $syntenies = $sra->fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss, $df);

  next unless (scalar @{$syntenies});
  print STDERR "For DnaFrag ".$df->name.", length ",$df->length,", ";
  print STDERR "got features " . scalar @{$syntenies} . "\n";
  $total_nb_syntenies += scalar @{$syntenies};

  my $dfname = $df->name;

  foreach my $sr (@{$syntenies}) {
    my ($species1_dfr_string, $species2_dfr_string , $species3_dfr_string,$species4_dfr_string);
    foreach my $dfr (@{$sr->get_all_DnaFragRegions}) {
      my $strand = "+";
      my ($gff_start,$gff_stop)=($dfr->dnafrag_start,$dfr->dnafrag_end);
      
      if ($dfr->dnafrag_strand < 0) {
        $strand = "-";
        ($gff_stop,$gff_start)=($dfr->dnafrag_start,$dfr->dnafrag_end);
      }
      
      if ($dfr->dnafrag->genome_db->name eq $species1) {
        $species1_dfr_string = $dfr->dnafrag->name . "\t" .
          "synteny\t" .
          "similarity\t" .
          $dfr->dnafrag_start . "\t" .
          $dfr->dnafrag_end . "\t" .
          "0.0" . "\t" .
          $strand . "\t" .
          ".\t" ;
      } elsif ($dfr->dnafrag->genome_db->name eq $species2) {
        $species2_dfr_string = 'Target "Sequence:'.$dfr->dnafrag->name . "\"\t".
          "$gff_start\t$gff_stop\tSpecies: \"$species2\"";
      } elsif ($dfr->dnafrag->genome_db->name eq $species3) {
        $species3_dfr_string = 'Target "Sequence:'.$dfr->dnafrag->name . "\"\t".
        "$gff_start\t$gff_stop\tSpecies: \"$species3\"";
      } elsif ($dfr->dnafrag->genome_db->name eq $species4) {
        $species4_dfr_string = 'Target "Sequence:'.$dfr->dnafrag->name . "\"\t".
        "$gff_start\t$gff_stop\tSpecies: \"$species4\"";
      }


    }
    print $species1_dfr_string . $species2_dfr_string."\n" if $species2_dfr_string;
    print $species1_dfr_string . $species3_dfr_string."\n" if $species3_dfr_string;
    print $species1_dfr_string . $species4_dfr_string."\n" if $species4_dfr_string;
  }
}

print STDERR "Total number of_synteny regions $total_nb_syntenies\n";
