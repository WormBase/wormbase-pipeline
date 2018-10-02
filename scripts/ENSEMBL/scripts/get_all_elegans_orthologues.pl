#!/usr/bin/env perl
#===============================================================================
#
#         FILE:  get_all_elegans_orthologues.pl
#
#      CREATED:  03/08/06 13:26:19 BST (mh6@sanger.ac.uk)
#===============================================================================

use strict;
use IO::File;
use Getopt::Long;

use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use lib $ENV{CVS_DIR};

use Wormbase;
use Log_files;

my ($verbose,$new,$debug,$test,$store,$outfile,$other, $wormbase,%species,@species);

my $comparadb = 'worm_compara';
my $dbhost    = $ENV{'WORM_DBHOST'};
my $dbuser    = 'wormro';
my $dbport    = $ENV{'WORM_DBPORT'};

GetOptions(
  'database=s'   => \$comparadb,
  'dbhost=s'     => \$dbhost,
  'dbuser=s'     => \$dbuser,
  'dbport=s'     => \$dbport,
  'verbose'      => \$verbose,
  'debug=s'      => \$debug,
  'test'         => \$test,
  'store=s'      => \$store,
  'outfile=s'    => \$outfile,
  'species=s'    => \@species,
) || die("cant parse the command line parameter\n");

if ($store){
 $wormbase = Storable::retrieve($store)
      or croak("cannot restore wormbase from $store"); 
}else{
 $wormbase = Wormbase->new(
    -test    => $test,
    -debug   => $debug,
 )
}

$outfile ||= $wormbase->acefiles . '/compara.ace';

my %core_accessors = $wormbase->species_accessors;
$core_accessors{$wormbase->species} = $wormbase;
my %non_core_accessors = $wormbase->tier3_species_accessors;
my %all_accessors = (%core_accessors, %non_core_accessors);

my $log = Log_files->make_build_log($wormbase);

my $cds2wbgene = &get_commondata_for_core_species('cds2wbgene_id');

my $compara_db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -host   => $dbhost,
    -user   => $dbuser,
    -port   => $dbport,
    -dbname => $comparadb
) or die(@!);

my $gdb_adaptor = $compara_db->get_GenomeDBAdaptor;
my $member_adaptor = $compara_db->get_GeneMemberAdaptor();
my $homology_adaptor = $compara_db->get_HomologyAdaptor();

my @genome_dbs = @{$gdb_adaptor->fetch_all};
foreach my $gdb (@genome_dbs) {
  if (not exists $all_accessors{$gdb->name}) {
    $log->log_and_die("Could not find genome_db name " . $gdb->name . " in WormBase\n");
  }
  $species{$gdb->dbID} = $all_accessors{$gdb->name}->full_name;
}

my $outfh = IO::File->new($outfile,'w')||die(@!);

foreach my $gdb1 (sort { $a->name cmp $b->name } @genome_dbs) {

  if (@species and not grep { $gdb1->name eq $_ } @species) {
    $log->write_to("Skipping " . $gdb1->name . "\n") if $verbose;
    next;
  }

  $log->write_to("Processing " . $gdb1->name . "...\n") if $verbose;
  
  my (%homols);
  
  foreach my $gdb2 (sort { $a->name cmp $b->name } @genome_dbs) {
    
    $log->write_to("   Comparing to " . $gdb2->name . "...\n") if $verbose;
    
    my $mlss;
    if ($gdb1->dbID == $gdb2->dbID) {
      $mlss = $compara_db->get_MethodLinkSpeciesSetAdaptor->fetch_by_method_link_type_GenomeDBs('ENSEMBL_PARALOGUES', [$gdb1]);      
    } else {
      $mlss = $compara_db->get_MethodLinkSpeciesSetAdaptor->fetch_by_method_link_type_GenomeDBs('ENSEMBL_ORTHOLOGUES', [$gdb1, $gdb2]);
    }
    
    my @homologies = @{$homology_adaptor->fetch_all_by_MethodLinkSpeciesSet( $mlss )};
    
    foreach my $homology (@homologies) {
      
      my ($m1, $m2) = sort { $a->stable_id cmp $b->stable_id } @{ $homology->get_all_Members };
      
      if ($m1->genome_db->dbID != $gdb1->dbID) {
        # members have been returned in the wrong order, so swap them
        ($m2, $m1) = ($m1, $m2);
      }        
      
      my ($gid1, $gid2);
      if (exists $core_accessors{$gdb1->name}) {
        $log->log_and_die("Could not find WBGene id for " . $m1->stable_id . "\n") if not exists $cds2wbgene->{$m1->stable_id};
        $gid1 =  $cds2wbgene->{$m1->stable_id};
      } else {
        $gid1 = sprintf("%s_%s", $all_accessors{$gdb1->name}->ncbi_bioproject, $m1->gene_member->stable_id);
      }
      if (exists $core_accessors{$gdb2->name}) {
        $log->log_and_die("Could not find WBGene id for " . $m2->stable_id . "\n") if not exists $cds2wbgene->{$m2->stable_id};
        $gid2 =  $cds2wbgene->{$m2->stable_id};
      } else {
        $gid2 = sprintf("%s_%s", $all_accessors{$gdb2->name}->ncbi_bioproject, $m2->gene_member->stable_id);
      }
      
      if ($gdb1->dbID == $gdb2->dbID) {
        # we need to add the connection both ways, so that that evidence gets added to both
        $homols{$gid1}->{Paralog}->{$species{$gdb2->dbID}}->{$gid2} = 1;
        $homols{$gid2}->{Paralog}->{$species{$gdb1->dbID}}->{$gid1} = 1;
      } else {
        $homols{$gid1}->{Ortholog}->{$species{$gdb2->dbID}}->{$gid2} = 1;
      }
    } 
  }
  
  print $outfh "// Homologies for " . $gdb1->name . "\n\n";
  
  foreach my $g (sort keys %homols) {
    
    print $outfh "\nGene : \"$g\"\n";
    
    foreach my $tag_group (keys %{$homols{$g}}) {
      foreach my $spe (sort keys %{$homols{$g}->{$tag_group}}) {
        foreach my $entity (sort keys %{$homols{$g}->{$tag_group}->{$spe}}) {
          print $outfh "$tag_group \"$entity\" \"$spe\" From_analysis WormBase-Compara\n";
        }
      }
    }
  }   
}

$log->mail;
exit(0);

##################################
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
