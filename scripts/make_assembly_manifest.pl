#!/usr/bin/env perl
#
# make_assembly_manifest.pl
#
# Makes the ASSEMBLIES.json file summarising assembly info for a particular release
# 
##########################################################################################

use strict;
use Getopt::Long;

use Storable;
use Ace;
use JSON;
use Wormbase;

my ($test,
    $debug,
    $store,
    $database,
    $wormbase,
    $manifile,
    );

GetOptions (
  "test"            => \$test,
  "debug=s"         => \$debug,
  "store:s"         => \$store,
  "database:s"      => \$database,
  "outfile:s"       => \$manifile,
    );


if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
      );
}

my $log = Log_files->make_build_log($wormbase);

$manifile = $wormbase->reports . "/ASSEMBLIES.json" if not defined $manifile;

open(my $fh, ">$manifile\n") or $log->log_and_die("Could not open $manifile for writing\n");
my $db = Ace->connect(-path => defined $database ? $database : $wormbase->autoace );

my (%accessors, %accessors_by_species, %json);

%accessors = ($wormbase->all_species_accessors);
$accessors{$wormbase->species} = $wormbase;

foreach my $acc (values %accessors) {
  my $species_name = $acc->full_name;
  push @{$accessors_by_species{$species_name}}, $acc;
}

foreach my $species (sort keys %accessors_by_species) {
  my (@accs) = @{$accessors_by_species{$species}};
  
  my $species_obj = $db->fetch(-class => 'Species', -name => "$species");
  my @seq_col = $species_obj->at('Assembly');
  
  my $g_species = $accs[0]->full_name(-g_species => 1); 
  
  my $obj = {
    full_name => $species,
    assemblies => [],
  };
  
  $json{$g_species} = $obj;
  
  foreach my $seq_col_name (@seq_col) {
    my $seq_col = $seq_col_name->fetch;
    
    my $strain = $seq_col->Strain;
    my $assembly_name = $seq_col->Name;
    unless (defined $assembly_name) {
      $log->write_to("ERROR: The $g_species assembly appears to be unnamed so the JSON file will reflect this....this DOES have implications for the website, so best to check and patch the build if there is an official assembly name and you have time (a patch can be supplied to the web team so you don't have to re-package from scratch).\n");
    }
    my $first_ws_rel = $seq_col->First_WS_release;
    my @labs;
    
    if ($seq_col->at('Origin.Laboratory')) {      
      my @laboratory = $seq_col->at('Origin.Laboratory');
      foreach my $lab (@laboratory) {
        push @labs, $lab->fetch->Mail->name;
      }
    }
    
    my ($bioproj, $gc_acc, $ucsc_name);
    
    my @db = $seq_col->at('Origin.DB_info.Database');
    foreach my $db (@db) {
      if ($db->name eq 'NCBI_BioProject') {
        $bioproj = $db->right->right->name;
      } elsif ($db->name eq 'NCBI_Genome_Assembly') {
        $gc_acc  = $db->right->right->name;
      } elsif ($db->name eq 'UCSC_Genome_Browser') {
        $ucsc_name = $db->right->right->name;
      }
    }
    
    # need to find the corresponding accessor, because only that
    # will tell us whether the bioproject is the canonical one
    my ($rel_acc) = grep { $_->ncbi_bioproject eq $bioproj } @accs;
    next unless $rel_acc;
    my $is_canonical = $rel_acc->is_canonical;
    my $bioproj_desc = $rel_acc->bioproject_description;
    
    push @{$obj->{assemblies}}, {
      bioproject => $bioproj,
      bioproject_description => $bioproj_desc,
      assembly_accession => $gc_acc,
      assembly_name => (defined $assembly_name) ? $assembly_name->name : "$g_species Assembly unnamed",
      ucsc_name => (defined $ucsc_name) ? $ucsc_name : undef, 
      appeared_in => 'WS'.$first_ws_rel->name,
      is_canonical => ($is_canonical) ? JSON::true : JSON::false,
      strain => (defined $strain) ? $strain->name : "Unknown strain",
      laboratory => \@labs,
    };
  }
}

$db->close;

my $json_obj = JSON->new;
my $string = $json_obj->allow_nonref->canonical->pretty->encode(\%json);

print $fh $string;
close($fh);

$log->mail();
exit(0);
