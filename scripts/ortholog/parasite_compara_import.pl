#!/usr/bin/env perl

# small script to get the ParaSite Orthologs into Ace format (including stub genes)
# usage: perl parasite_compara_import.pl [-debug <username>] [-outfile <filename>]


use lib $ENV{CDS_DIR};

use Bio::EnsEMBL::Registry;
use Getopt::Long;
use Wormbase;
use Log_files;

my ($outfile,$debug);
GetOptions(
   'outfile=s' => \$outfile,
   'debug=s'   => \$debug,
)||die(@!);

my $wormbase=Wormbase->new(-debug => $debug);

$outfile||=$wormbase->acefiles . '/parasiteCompara.ace';
my $ofh = new IO::File $outfile,'w';

my %core_accessors = $wormbase->species_accessors;
$core_accessors{$wormbase->species} = $wormbase;
my @core = map {$_->full_name} values %core_accessors;
my %non_core_accessors = $wormbase->tier3_species_accessors;
my @noncore = map {$_->full_name} values %non_core_accessors;

my @productionNames;
foreach my $accessor (values (%core_accessors),values(%non_core_accessors)){
   push @productionNames, $accessor->full_name(-g_species=>1);
   my $long_name = lc $accessor->full_name;
   $long_name =~ s/\s/_/;
   push @productionNames, $long_name;
   $long_name.='_'. lc $accessor->ncbi_bioproject;
}

my @others = ('Danio rerio','Drosophila melanogaster','Mus musculus','Homo sapiens','Saccharomyces cerevisiae');
my @all = (@core,@noncore,@others);
my @wb = (@core,@noncore);

my $log = Log_files->make_build_log($wormbase);

# needs to pull it in from the WormBase module
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_multiple_dbs( 
      {
        -host    => 'mysql-ps-staging-2.ebi.ac.uk',
        -user    => 'ensro',
        -port    => 4467,
      },
      {
        -host    => 'mysql-ps-prod-1.ebi.ac.uk',
        -user    => 'ensro',
        -port    => 4450,
      }
);

my $cds2wbgene = &get_commondata_for_core_species('worm_gene2geneID_name');

my $homology_adaptor = $registry->get_adaptor('wbparasite','compara','Homology');
my $genomeDB_adaptor = $registry->get_adaptor('wbparasite','compara','GenomeDB');

my @genome_dbs = @{$genomeDB_adaptor->fetch_all};

my %species;
my %wbspecies;
my %gdb2taxonid;
foreach my $gdb(@genome_dbs){
    my $gdbtaxon = $gdb->taxon->binomial||($gdb->taxon->name=~/(\w+\s+\w+)/g)[0];
    $gdbtaxon='Caenorhabditis tropicalis' if $gdb->taxon_id == 886184; # csp11 hack
    if (grep {$gdbtaxon eq $_} @others){
      $species{$gdb->dbID} = $gdbtaxon;
    }elsif(grep {$gdb->name eq $_} @productionNames){
      $species{$gdb->dbID} = $gdbtaxon;
      $wbspecies{$gdb->dbID} = $gdbtaxon;
    }
    $gdb2taxonid{$gdb->dbID}=$gdb->taxon_id;
    $gdb2taxonid{$gdb->dbID}=1094322 if $gdb->taxon_id == 886184; # csp11 hack
}

foreach my $gdb1 (@genome_dbs) {
     next unless $wbspecies{$gdb1->dbID}; # skip non-wormbase species

     my %homols;

     foreach my $gdb2 (@genome_dbs) {
      next unless $species{$gdb2->dbID}; # skip all the ones we don't like

      my $mlss;
      if ($gdb1->dbID == $gdb2->dbID) {
       $mlss = $registry->get_adaptor('wbparasite','compara','MethodLinkSpeciesSet')->fetch_by_method_link_type_GenomeDBs('ENSEMBL_PARALOGUES', [$gdb1]);      
      } else {
       $mlss = $registry->get_adaptor('wbparasite','compara','MethodLinkSpeciesSet')->fetch_by_method_link_type_GenomeDBs('ENSEMBL_ORTHOLOGUES', [$gdb1, $gdb2]);
      }

      foreach my $homology (@{$homology_adaptor->fetch_all_by_MethodLinkSpeciesSet($mlss)}) {
       my ($m1, $m2) = @{ $homology->get_all_Members };
      
       if ($m1->genome_db->dbID != $gdb1->dbID) {
        ($m2, $m1) = ($m1, $m2); # members have been returned in the wrong order, so swap them
       }        
      
       my ($gid1, $gid2);
       if (grep {$species{$gdb1->dbID} eq $_} @core) {
         $gid1 =  $m1->gene_member->stable_id=~/WBGene\d+/ ? $m1->gene_member->stable_id:$cds2wbgene->{$m1->gene_member->stable_id};
        unless ($gid1){
         $log->write_to("skipping ${\$m1->gene_member->stable_id} (doesn't exist in current version)\n");
         next;
        }
       }else{
        $gid1 = sprintf('%s:%s',$gdb2taxonid{$gdb1->dbID},$m1->gene_member->stable_id);
       } 

       if (grep {$species{$gdb2->dbID} eq $_} @core) {
          $gid2 =  $m2->gene_member->stable_id=~/WBGene\d+/ ? $m2->gene_member->stable_id : $cds2wbgene->{$m2->gene_member->stable_id};
          unless ($gid2){
            $log->write_to("skipping ${\$m2->gene_member->stable_id} (doesn't exist in current version)\n");
            next;
          }
       } elsif(grep{ $species{$gdb2->dbID} eq $_ } @noncore){
          $gid2 = sprintf('%s:%s',$gdb2taxonid{$gdb2->dbID},$m2->gene_member->stable_id);
       } else{
          $gid2 = $m2->gene_member->stable_id;
       }
 
       if ($gdb1->dbID == $gdb2->dbID) {
        # we need to add the connection both ways, so that that evidence gets added to both
        $homols{$gid1}->{Paralog}->{$species{$gdb2->dbID}}->{$gid2} = 1;
        $homols{$gid2}->{Paralog}->{$species{$gdb1->dbID}}->{$gid1} = 1;
       } else {
        $homols{$gid1}->{Ortholog}->{$species{$gdb2->dbID}}->{$gid2} = 1; # save forward
        $homols{$gid2}->{Ortholog}->{$species{$gdb1->dbID}}->{$gid1} = 1; # and reverse
       }
    } # homologies
  } # gdb2
  print  $ofh "// Homologies for ${\$gdb1->name}\n";
  
  foreach my $g (keys %homols) {
    
     print $ofh "Gene : \"$g\"\n";

     foreach my $tag_group (keys %{$homols{$g}}) {
       foreach my $spe (sort keys %{$homols{$g}->{$tag_group}}) {
         foreach my $entity (sort keys %{$homols{$g}->{$tag_group}->{$spe}}) {
           print $ofh "$tag_group \"$entity\" \"$spe\" From_analysis ParaSite-Compara\n";
         }
       }
    }
    print $ofh "\n";
  }
}

$log->mail;

# slurps in all species commondata
sub get_commondata_for_core_species {
  my ($name)=@_;
  
  my %all_data;
  
  foreach my $wb (values %core_accessors) {
    my %hash;
    $wb->FetchData($name, \%hash);
    foreach my $k (keys %hash) {
      $all_data{$k} = $hash{$k};
    }
  }
  return \%all_data;
}

