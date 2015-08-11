#!/usr/bin/env perl

use strict;
use Getopt::Long;
use LWP::UserAgent;

use lib $ENV{CVS_DIR};

use Wormbase;
use Log_files;

my ($debug, $test, $store, $species, $wb, $overwrite, $noload,
    $pid_table_file, $outace);


&GetOptions ("debug:s"        => \$debug,
             "test"           => \$test,
             "store:s"        => \$store,
             "species:s"      => \$species,
             "noload"         => \$noload,
             "pidtable=s"     => \$pid_table_file,
             'acefile=s'      => \$outace,
             'overwrite=s'    => \$overwrite,
    );


if ($store) { 
  $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n"); 
}
else { 
  $wb = Wormbase->new(
    -debug    => $debug, 
    -test     => $test,
    -organism => $species, 
      );
}

my $log = Log_files->make_build_log($wb);

$pid_table_file = $wb->acefiles . "/EBI_protein_ids.txt" if not defined $pid_table_file;
$outace = $wb->acefiles . "/EBI_Uniprot_brief_ids.ace" if not defined $outace;

my $cur_data = &fetch_current_data();

my (%cds_product); 

open(my $table_fh, $pid_table_file)
    or $log->log_and_die("Could not open $pid_table_file for reading\n");
while(<$table_fh>) {
  chomp;
  my @data = split("\t",$_);
  
  next unless scalar(@data) == 11;
  my( $cds, $uniprot_ac, $product_name) 
      = ($data[5],$data[6],$data[10]);  
  
  if (defined $product_name and $product_name ne '.' and defined $uniprot_ac and $uniprot_ac ne'.') {
    $cds_product{$cds} = {
      product => $product_name, 
      uniprot => $uniprot_ac,
    };
  }
}
close($table_fh) or $log->log_and_die("Could not close the protein_id command/file\n");


open(my $acefh, ">$outace")
    or $log->log_and_die("Could not open $outace for writing\n");

foreach my $gene (sort keys %$cur_data) {
  my %brief_ids;
  
  foreach my $cds (sort keys %{$cur_data->{$gene}}) {
    if (exists $cds_product{$cds}) {
      my $prod = $cds_product{$cds}->{product};
      my $uniprot = $cds_product{$cds}->{uniprot};
      
      my $update = 0;
      
      if (exists $cur_data->{$gene}->{$cds}->{brief_id}) {
        if ($cur_data->{$gene}->{$cds}->{brief_id} eq $prod) {
          if (exists $cur_data->{$gene}->{$cds}->{accession} and
              $cur_data->{$gene}->{$cds}->{accession} eq $uniprot) {
            $log->write_to("$cds : attached to UniProt $uniprot and product unchanged. Doing nothing\n") if $debug;
          } else {
            $log->write_to("$cds : UniProt product matches, but accession does not. Updating.\n");
            $update = 1;          
          }
        } else {
          $log->write_to("$cds : UniProt product is different from current. Updating.\n");
          $update = 1;
        }
      } else {
        $log->write_to("$cds : no brief_id. Updating.\n");
        $update = 1;
      }

      if ($update) {
        $prod =~ s/\//\\\//g; 
        print $acefh "\nCDS : \"$cds\"\n";
        print $acefh "-D Brief_identification\n";
        print $acefh "\nCDS : \"$cds\"\n";
        print $acefh "Brief_identification \"$prod\" Accession_evidence \"UniProt\" \"$uniprot\"\n";      
        $cur_data->{$gene}->{$cds} = {
          brief_id => $prod,
          evidence_tag => "Accession_evidence", 
          database => "UniProt",
          accession => $uniprot,
        };
      }
    } else {
      if (exists $cur_data->{$gene}->{$cds}->{brief_id}) {
        $log->write_to("$cds : No UniProt product, but existing product. Doing nothing\n") if $debug; 
      } else {
        $log->write_to("$cds : No UniProt product and no existing product. Doing nothing\n") if $debug; 
      }
    }
  }
}


close($acefh) or $log->log_and_die("Could not close $outace properly\n");

unless ($noload) {
  $wb->load_to_database($wb->autoace, $outace, 'ENA_protein_xrefs', $log);
}

$log->mail();
exit(0);


######################################################
sub fetch_current_data {
  my (%gene2brief);

  $log->write_to("Fetching current data...\n");
  my $db = Ace->connect(-path => $wb->autoace );
  
  my $full_name = $wb->full_name;
  
  my $iter = $db->fetch_many(-query => "FIND CDS Method = \"curated\" AND Species = \"$full_name\"");
  while (my $cds = $iter->next) {
    my $gene = $cds->Gene->name;

    my $data = { };

    my ($bid) = $cds->Brief_identification;
    if ($bid) {
      $data->{brief_id} = $bid->name;
      
      my ($evi_tag, $db, $acc);
      if ($bid->right) {
        my ($tag) = $bid->right;
        $data->{evidence_tag} = $tag->name;
        if ($tag eq 'Accession_evidence') {
          $data->{database} = $tag->right->name;
          $data->{accession} = $tag->right->right->name;
        }
      }
    }
    $gene2brief{$gene}->{$cds->name} = $data;
  }
  
  $db->close();
  
  return \%gene2brief;
}


1;
