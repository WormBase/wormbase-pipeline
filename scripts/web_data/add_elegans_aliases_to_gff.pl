#!/usr/bin/env perl

use strict;
use Getopt::Long;

use lib $ENV{CVS_DIR};
use Wormbase;
use Storable;
use Log_files;

use Bio::EnsEMBL::Registry;

my ($gff_source, 
    $gff_type,
    $gff3,
    $test,
    $debug,
    $species,
    $new_lines,
    $store,
    $ens_regconf,
    $wormbase,
    );

&GetOptions('gffsource=s' => \$gff_source,
            'gfftype=s'   => \$gff_type,
            'newlines=s'  => \$new_lines,
            'species=s'   => \$species,
            'test'        => \$test,
            'debug=s'     => \$debug,
            'ensreg=s'    => \$ens_regconf,
            );

$gff_type = 'CDS' if not defined $gff_type;
$gff_source = 'curated' if not defined $gff_source;

my $reg = "Bio::EnsEMBL::Registry";
$reg->load_all($ens_regconf);

my $member_adaptor = $reg->get_adaptor('compara_trees', 'compara', 'Member');
my $homology_adaptor = $reg->get_adaptor('compara_trees', 'compara', 'Homology');
my $genomedb_adaptor = $reg->get_adaptor('compara_trees', 'compara', 'GenomeDB');
my $mlss_adaptor = $reg->get_adaptor('compara_trees', 'compara', 'MethodLinkSpeciesSet');

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
                             );
}

my $log = Log_files->make_build_log($wormbase);

my %sa = $wormbase->species_accessors;

my %cds2wbgeneid; 
# the following will fail for tierIIIs, so stick in in an eval
eval { 
  %cds2wbgeneid = $wormbase->FetchData('cds2wbgene_id');
};
my %elegans_cds2cgc = $sa{elegans}->FetchData('cds2cgc');
my %elegans_cds2wbgeneid = $sa{elegans}->FetchData('cds2wbgene_id');

my $elegans_tax_id = $sa{elegans}->ncbi_tax_id;
my $this_tax_id = $wormbase->ncbi_tax_id;

my $elegans_genome_db = $genomedb_adaptor->fetch_by_taxon_id($elegans_tax_id);
my $this_genome_db = $genomedb_adaptor->fetch_by_taxon_id($this_tax_id);

my $mlss = $mlss_adaptor->fetch_by_method_link_type_GenomeDBs('ENSEMBL_ORTHOLOGUES',
                                                              [$elegans_genome_db,
                                                               $this_genome_db]);

my %elegans_aliases;

if ($elegans_tax_id ne $this_tax_id) {
  my @members = @{$member_adaptor->fetch_all_by_source_taxon("ENSEMBLGENE", $this_tax_id)};

  $log->write_to(sprintf("Fetched %d members for source=ENSEMBLGENE and taxid=$this_tax_id\n", scalar(@members)));;

  @members = @members[0..100];

  while( my $member = shift @members){
    
    my $name = $member->stable_id;
    if (exists $cds2wbgeneid{$name}) {
      $name = $cds2wbgeneid{$name};
    }

    my @homologies = @{$homology_adaptor->fetch_all_by_Member_MethodLinkSpeciesSet( $member, $mlss)};
    $log->write_to(sprintf("  Got %d homologies for %s\n", scalar(@homologies), $name));  

    foreach my $homology ( @homologies) {  

      my ($cgc_alias, $seq_name_alias, $wbgene_alias); 

      foreach my $ma ( @{ $homology->get_all_Member_Attribute } ) {
        my ( $me, $at ) = @{$ma};

        next if $me->taxon_id ne $elegans_tax_id;
        
        # Sequence_name
        $seq_name_alias = $me->stable_id;

        my ($pepm) = @{ $me->get_all_peptide_Members() };

        my $elegans_cds_id = $pepm->stable_id;

        $wbgene_alias = $elegans_cds2wbgeneid{ $elegans_cds_id };
        $cgc_alias = $elegans_cds2cgc{ $elegans_cds_id };
      }
      
      foreach my $alias ($cgc_alias, $seq_name_alias, $wbgene_alias) {
        if ($alias) {
          push @{$elegans_aliases{$name}}, $alias;
        }
      }

    }

  }
}


while(<>) {
  if (/^\#/) {
    print if not $new_lines;
    next;
  }

  chomp;

  my ( $ref, $source, $type, $start, $stop, $score, $strand, $phase, $group ) = split /\t/;
  
  if ( $source eq $gff_source && $type eq $gff_type){
    my $gene_id;

    if ($group =~ /ID=gene:([^;]+)/) {
      $gff3 = 1;
      $gene_id = $1;
    } elsif ($group =~ /Gene\s+\"(\S+)\"/) {
      $gene_id = $1;
    } else {
      die("Could not extract gene id from group field '$group'\n");
    }

    if (exists $elegans_aliases{$gene_id}) {
      foreach my $alias (@{$elegans_aliases{$gene_id}}) {
        my $alias_string = ($gff3) ? ";Alias=$alias" : " ; Alias \"$alias\"";
     
        if ($group !~ /$alias_string/) {
          $group .= "$alias_string";
        }
      }

    
      my ($nsrc, $ntype);
      
      if ($new_lines) {
        ($nsrc, $ntype) = split(/:/, $new_lines);
      } else {
        ($nsrc, $ntype) = ($source, $type);
      }

      print join("\t", $ref, $nsrc, $ntype, $start, $stop, $score, $strand, $phase, $group), "\n";
    } else {
      if (not $new_lines) {
        print join("\t", $ref, $source, $type, $start, $stop, $score, $strand, $phase, $group), "\n";
      }
    }

    next;
  }
  
  if (not $new_lines) {
    print join("\t", $ref, $source, $type, $start, $stop, $score, $strand, $phase, $group), "\n";
  }

}

$log->mail();

exit(0);
