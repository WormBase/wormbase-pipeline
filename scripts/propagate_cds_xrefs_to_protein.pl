#!/usr/bin/env perl

use strict;
use Getopt::Long;

use lib $ENV{CVS_DIR};

use Wormbase;
use Log_files;

my ($debug, $test, $store, $species, $wb, $noload, $acefile, $common_data_dir);

&GetOptions ("debug:s"        => \$debug,
             "test"           => \$test,
             "store:s"        => \$store,
             "species:s"      => \$species,
             "acefile=s"      => \$acefile,
             "noload"         => \$noload,
             'commondata:s'   => \$common_data_dir,
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


my (%cds2wp, %pep_xrefs);

$acefile = $wb->acefiles . "/protein_object_xrefs.ace" if not defined $acefile;

$wb->FetchData('cds2wormpep', \%cds2wp, $common_data_dir);

my $tm_def = &write_CDS_xref_def();
my $tm_query = $wb->table_maker_query( $wb->autoace, $tm_def );

while(<$tm_query>) {
  chomp;  
  s/\"//g; 
  next if (/acedb/ or /\/\//);
  next if /^\s*$/;
  
  my ($cds_id, $db_name, $field, $acc) = split(/\t/, $_);

  if (exists $cds2wp{$cds_id}) {
    my $pep_acc = $cds2wp{$cds_id};
    $pep_xrefs{$pep_acc}->{$db_name}->{$field}->{$acc} = 1; 
  }
}
unlink $tm_def;

open(my $out_fh, ">$acefile") or $log->log_and_die("Could not open $acefile for writing\n");
foreach my $pep (sort keys %pep_xrefs) {
  print $out_fh "\nProtein : \"$pep\"\n";
  foreach my $db_name (sort keys %{$pep_xrefs{$pep}}) {
    foreach my $field (sort keys %{$pep_xrefs{$pep}->{$db_name}}) {
      foreach my $acc (sort keys %{$pep_xrefs{$pep}->{$db_name}->{$field}}) {
        print $out_fh "Database $db_name $field $acc\n";
      }
    }
  }
}
close($out_fh) or $log->log_and_die("Could not close $acefile after writing\n");

unless ($noload) {
  $wb->load_to_database($wb->autoace, $acefile, 'protein_object_xrefs', $log);
}


$log->mail();
exit(0);


sub write_CDS_xref_def {
  my $txt = <<"ENDE";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class CDS 
From 1 
Condition Method = "curated"
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class Database 
From 1 
Tag Database 
 
Colonne 3 
Width 12 
Optional 
Visible 
Class 
Class Database_field 
Right_of 2 
Tag  HERE  
 
Colonne 4 
Width 12 
Optional 
Visible 
Class 
Class Text 
Right_of 3 
Tag  HERE  

ENDE

  return &write_tm_def("CDS_xrefs", $txt);
}



###################################
sub write_tm_def {
  my ($fname, $string) = @_;
  
  my $file = "/tmp/$fname.def";

  open(my $fh, ">$file") or $log->log_and_die("Could not open $fname for writing\n");
  print $fh $string;
  close($fh) or $log->log_and_die("Could not close $fname after writing\n");

  return $file;
}
