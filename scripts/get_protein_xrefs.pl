#!/usr/bin/env perl

use strict;
use Getopt::Long;
use LWP::UserAgent;

use lib $ENV{CVS_DIR};

use Wormbase;
use Log_files;

my ($debug, $test, $store, $species, $wb, 
    $ncbi_tax_id, $bioproject_id, $pid_table_file);


&GetOptions ("debug:s"        => \$debug,
             "test"           => \$test,
             "store:s"        => \$store,
             "species:s"      => \$species,
             "pidtable=s"     => \$pid_table_file,
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


my ($ggenus, $gspecies) = $wb->full_name =~ /^(\S+)\s+(\S+)/;

$ncbi_tax_id = $wb->ncbi_tax_id;
$bioproject_id = $wb->ncbi_bioproject;
$pid_table_file = $wb->acefiles . "/EBI_protein_ids.txt" if not defined $pid_table_file;

&lookup_from_ebi_production_dbs($pid_table_file); 

exit(0);

#########################
sub lookup_from_ebi_production_dbs {
  my ($output_file, $type) = @_;

  my $ebi_prod_dir = $wb->wormpub . "/ebi_resources";
  my $ena_env      = "$ebi_prod_dir/ena_oracle_setup.sh";
  my $ena_cred     = "$ebi_prod_dir/ENAORACLE.s";
  my $uni_cred     = "$ebi_prod_dir/UNIPROTORACLE.s";

  my $cmd =  "source $ena_env &&"
      . " perl  $ENV{CVS_DIR}/get_protein_ids_ebiprod.pl"
      . "  -enacred $ena_cred" 
      . "  -uniprotcred $uni_cred"
      . "  -orgid $ncbi_tax_id"
      . "  -bioprojectid $bioproject_id";
  
  system("$cmd > $output_file") 
        and die("Could not successfully run '$cmd'\n");
  
}

