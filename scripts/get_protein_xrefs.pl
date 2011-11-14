#!/usr/bin/env perl

use strict;
use Getopt::Long;
use LWP::UserAgent;

use lib $ENV{CVS_DIR};

use Wormbase;
use Log_files;

my ($debug, $test, $store, $species, $wb, $acefile, $load, $ncbi_tax_id, $table_file,
    %cds_xrefs, %pep_xrefs);


&GetOptions ("debug=s"    => \$debug,
             "test"       => \$test,
             "store:s"    => \$store,
             "species:s"  => \$species,
             "load"       => \$load,
             "table=s"    => \$table_file,
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

my %cds2wormpep = $wb->FetchData('cds2wormpep');
my %accession2clone   = $wb->FetchData('accession2clone');

my ($ggenus, $gspecies) = $wb->full_name =~ /^(\S+)\s+(\S+)/;

$ncbi_tax_id = $wb->ncbi_tax_id;
$acefile = $wb->acefiles . "/cds_embl_data.ace" if not defined $acefile;

my $table_fh;
if (defined $table_file) {
  open($table_fh, $table_file);
} else {
  $table_fh = &lookup_from_ebi_production_dbs();
}

open(my $acefh, ">$acefile")
    or $log->log_and_die("Could not open $acefile for writing\n");

while(<$table_fh>) {
  chomp;
  my @data = split("\t",$_);
  
  next unless scalar(@data) == 8;
  my($cloneacc, $pid, $version, $cds, $uniprot) = ($data[0],$data[2],$data[3],$data[-1],$data[-2]);  
  
  next unless (defined $pid);
  $log->write_to("Potential New Protein: $_\n") if $uniprot eq 'UNDEFINED';
  
  next unless $accession2clone{$cloneacc}; #data includes some mRNAs
  
  push @{$cds_xrefs{$cds}->{Protein_id}}, [$accession2clone{$cloneacc}, $pid, $version];

  if($cds2wormpep{$cds} and defined $uniprot and $uniprot ne 'UNDEFINED') {
    $cds_xrefs{$cds}->{UniProtAcc}->{$uniprot} = 1;
    $pep_xrefs{"WP:".$cds2wormpep{$cds}}->{UniProtAcc}->{$uniprot} = 1;
  }
}
close($table_fh) or $log->log_and_die("Could not close the protein_id command/file\n");

my $acc2idmap = &get_uniprot_acc2idmap();

foreach my $pair (["CDS",\%cds_xrefs], ["Protein", \%pep_xrefs]) {
  my ($class, $hash) = @$pair;

  foreach my $k (keys %$hash) {
    print $acefh "$class : \"$k\"\n";
    if (exists $hash->{$k}->{Protein_id}) {
      foreach my $pidl (@{$hash->{$k}->{Protein_id}}) {
        print $acefh "Protein_id\t@$pidl\n"; 
      }
    }

    if (exists $hash->{$k}->{UniProtAcc}) {
      foreach my $acc (keys %{$hash->{$k}->{UniProtAcc}}) {
        print $acefh "Database UniProt UniProtAcc $acc\n";

        if (exists $acc2idmap->{$acc}) {
          printf $acefh "Database UniProt UniProtID %s\n", $acc2idmap->{$acc};
        }
      }
    }
    print $acefh "\n";
  }
}

close($acefh) or $log->log_and_die("Could not close $acefile properly\n");

if ($load) {
  $wb->load_to_database($wb->autoace, $acefile, 'get_protein_xrefs', $log);
}

$log->mail();
exit(0);

#########################
sub lookup_from_ebi_production_dbs {

  my $ebi_prod_dir = $wb->wormpub . "/ebi_resources";
  my $ena_perl     = "$ebi_prod_dir/ena_perl";
  my $ena_env      = "$ebi_prod_dir/ena_oracle_setup.sh";
  my $swiss_idx    = "$ebi_prod_dir/swiss_idx.txt";
  my $trembl_idx   = "$ebi_prod_dir/trembl_idx.txt";
  
  my $cmd =  "source $ena_env &&"
      . " $ena_perl  $ENV{CVS_DIR}/get_protein_ids_ebiprod.pl"
      . "  -swissidx $swiss_idx"
      . "  -tremblidx $trembl_idx"
      . "  -enadb ENAPRO" 
      . "  -orgid $ncbi_tax_id";
  
  open(my $cmdfh, "$cmd |")
      or $log->log_and_die("Could not successfully invoke '$cmd'\n");

  return $cmdfh;
}


###########################
sub get_uniprot_acc2idmap {
  my (%acc2ids);

  my $ua       = LWP::UserAgent->new;

  my $base     = 'http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-noSession';
  my $query    = "+[uniprot-org:$ggenus]&[uniprot-org:$gspecies]";

  my $cResult  = '+-page+cResult+-ascii';
  my $fullview = '+-view+UniprotView+-ascii+-lv+';

  $log->write_to("Doing query: ${base}${query}${cResult}\n");


  my $qa1 = $ua->get($base.$query.$cResult);
  $log->log_and_die("Can't get URL -- " . $qa1->status_line) unless $qa1->is_success;

  if($qa1->content =~/^(\d+)/) {
    my $lv = $1;
    $log->write_to("EBI SRS server returned $lv entries; fetching...\n");
    
    my $tmp_file = "/tmp/srs_results.$$.txt";
    my $qa2 = $ua->get($base.$query.$fullview.$lv, ':content_file' => $tmp_file);
    $log->log_and_die("Could not fetch Uniprot entries using EBI SRS server") 
        if not $qa2->is_success;

    open(my $f, $tmp_file);
    while(<$f>) {
      /UNIPROT:(\S+)\s+(\S+)/ and do {
        $acc2ids{$2} = $1;
      }
    }
  } else {
    $log->log_and_die("Unexpected content from SRS query\n");
  }
   
  return \%acc2ids;
}
