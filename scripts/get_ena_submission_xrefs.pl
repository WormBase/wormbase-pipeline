#!/usr/bin/env perl

use strict;
use Getopt::Long;
use LWP::UserAgent;

use lib $ENV{CVS_DIR};

use Wormbase;
use Log_files;

my ($debug, $test, $store, $species, $wb, 
    $svacefile, $pidacefile, $noload, $ncbi_tax_id, $bioproject_id, $pid_table_file,$sv_table_file, $gid_table_file, $gidacefile,
    %cds_xrefs, %pep_xrefs, $generate_tables, $sequence_xrefs, $protein_xrefs, $gene_xrefs, $common_data_dir);


&GetOptions ("debug:s"        => \$debug,
             "test"           => \$test,
             "store:s"        => \$store,
             "species:s"      => \$species,
             "noload"         => \$noload,
             "svtable=s"      => \$sv_table_file,
             'svacefile=s'    => \$svacefile,
             "pidtable=s"     => \$pid_table_file,
             'pidacefile=s'   => \$pidacefile,
             'gidtable=s'     => \$gid_table_file,
             'gidacefile=s'   => \$gidacefile,
             'generatetables' => \$generate_tables,
             'sequencexrefs'  => \$sequence_xrefs,
             'proteinxrefs'   => \$protein_xrefs,
             'genexrefs'      => \$gene_xrefs,
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

my ($ggenus, $gspecies) = $wb->full_name =~ /^(\S+)\s+(\S+)/;

$ncbi_tax_id = $wb->ncbi_tax_id;
$bioproject_id = $wb->ncbi_bioproject;
$common_data_dir = $wb->common_data if not defined $common_data_dir;
$svacefile = $wb->acefiles . "/EBI_sequence_xrefs.ace" if not defined $svacefile;
$pidacefile = $wb->acefiles . "/EBI_pid_xrefs.ace" if not defined $pidacefile;
$gidacefile = $wb->acefiles . "/EBI_gene_xrefs.ace" if not defined $gidacefile;

$pid_table_file = $wb->acefiles . "/EBI_protein_ids.txt" if not defined $pid_table_file;
$sv_table_file = $wb->acefiles . "/EBI_sequence_versions.txt" if not defined $sv_table_file;
$gid_table_file = $wb->acefiles . "/EBI_gene_ids.txt" if not defined $gid_table_file;


if ($generate_tables) {
  &lookup_from_ebi_production_dbs($pid_table_file, 'proteinxrefs');
  &lookup_from_ebi_production_dbs($sv_table_file, 'seqversions');
  &lookup_from_ebi_production_dbs($gid_table_file, 'genexrefs');
}



if ($sequence_xrefs) {
  my %accession2clone;
  $wb->FetchData('accession2clone', \%accession2clone, $common_data_dir);
    
  open(my $vtable_fh, $sv_table_file) 
      or $log->log_and_die("Could not open $sv_table_file for reading\n");
  
  open(my $acefh, ">$svacefile")
      or $log->log_and_die("Could not open $svacefile for writing\n");
  
  while(<$vtable_fh>) {
    /^(\S+)\s+(\d+)/ and do {
      my ($cloneacc, $ver) = ($1, $2);
      
      my $clone =  $accession2clone{$cloneacc};
      next if not $clone;
      
      print $acefh "Sequence : \"$clone\"\n";
      print $acefh "-D Database EMBL  NDB_SV\n";
      print $acefh "\nSequence : \"$clone\"\n";
      print $acefh "Database EMBL NDB_SV $cloneacc.$ver\n\n";
    }
  }
  
  close($acefh) or $log->log_and_die("Could not close $svacefile properly\n");
  
  unless ($noload) {
    $wb->load_to_database($wb->autoace, $svacefile, 'ENA_sequence_xrefs', $log);
  }

}

if ($gene_xrefs) {

  open(my $gtable_fh, $gid_table_file) 
      or $log->log_and_die("Could not open $gid_table_file for reading\n");
  
  open(my $acefh, ">$gidacefile")
      or $log->log_and_die("Could not open $gidacefile for writing\n");

  my %tran2gene;
  $wb->FetchData('worm_gene2geneID_name', \%tran2gene, $common_data_dir);
  
  while(<$gtable_fh>) {
    my ($tran_name, $clone_acc, $locus_name, $locus_tag) = /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;

    if (exists $tran2gene{$tran_name} and $locus_tag ne '.' and $clone_acc ne '.') {
      my $gene = $tran2gene{$tran_name}; 
      print $acefh "\nGene : \"$gene\"\n";
      print $acefh "Other_name \"$locus_tag\" Accession_evidence \"NDB\" \"$clone_acc\"\n";
      print $acefh "Database \"NDB\" \"locus_tag\" \"$locus_tag\"\n";
    }
  }
  close($acefh) or $log->log_and_die("Could not close $gidacefile after writing (probably disk full)\n");

  unless ($noload) {
    $wb->load_to_database($wb->autoace, $gidacefile, 'ENA_gene_xrefs', $log);
  }
}



if ($protein_xrefs) {
  my ( %accession2clone, %cds2wormpep);
  $wb->FetchData('accession2clone', \%accession2clone, $common_data_dir);
  $wb->FetchData('cds2wormpep', \%cds2wormpep, $common_data_dir);

  open(my $table_fh, $pid_table_file)
      or $log->log_and_die("Could not open $pid_table_file for reading\n");

  open(my $acefh, ">$pidacefile")
      or $log->log_and_die("Could not open $pidacefile for writing\n");

  while(<$table_fh>) {
    chomp;
    my @data = split("\t",$_);
    
    next unless scalar(@data) == 10;
    my($cloneacc, $pid, $version, $cds, $uniprot_ac, $uniprot_id, $uniprot_iso_acc, $ec_num) 
        = ($data[0],$data[2],$data[3],$data[5],$data[6],$data[7],$data[8],$data[9]);  
    
    next unless (defined $pid);
    $log->write_to("Potential New Protein: $_\n") if $uniprot_ac eq 'UNDEFINED';
    
    next unless exists $cds2wormpep{$cds}; # just ENA is slightly out of date w.r.t. our latest annotation
    
    push @{$cds_xrefs{$cds}->{Protein_id}}, [$accession2clone{$cloneacc}, $pid, $version];
    
    if($cds2wormpep{$cds}) {
      if (defined $uniprot_ac and $uniprot_ac ne 'UNDEFINED') {
        $cds_xrefs{$cds}->{dblinks}->{UniProt}->{UniProtAcc}->{$uniprot_ac} = 1;
        $pep_xrefs{"WP:".$cds2wormpep{$cds}}->{dblinks}->{UniProt}->{UniProtAcc}->{$uniprot_ac} = 1;
      }
      if (defined $uniprot_id and $uniprot_id ne 'UNDEFINED') {
        $cds_xrefs{$cds}->{dblinks}->{UniProt}->{UniProtId}->{$uniprot_id} = 1;
        $pep_xrefs{"WP:".$cds2wormpep{$cds}}->{dblinks}->{UniProt}->{UniProtId}->{$uniprot_id} = 1;
      }
      if (defined $uniprot_iso_acc and $uniprot_iso_acc ne 'UNDEFINED') {
        $cds_xrefs{$cds}->{dblinks}->{UniProt}->{UniProtIsoformAcc}->{$uniprot_iso_acc} = 1;
        $pep_xrefs{"WP:".$cds2wormpep{$cds}}->{dblinks}->{UniProt}->{UniProtIsoformAcc}->{$uniprot_iso_acc} = 1;
      }
      if (defined $ec_num and $ec_num ne 'UNDEFINED') {
        $cds_xrefs{$cds}->{dblinks}->{KEGG}->{KEGG_id}->{$ec_num} = 1;
        $pep_xrefs{"WP:".$cds2wormpep{$cds}}->{dblinks}->{KEGG}->{KEGG_id}->{$ec_num} = 1;
      }
    }
  }
  close($table_fh) or $log->log_and_die("Could not close the protein_id command/file\n");
  
  foreach my $pair (["CDS",\%cds_xrefs], ["Protein", \%pep_xrefs]) {
    my ($class, $hash) = @$pair;
    
    foreach my $k (keys %$hash) {
      print $acefh "$class : \"$k\"\n";

      if (exists $hash->{$k}->{Protein_id}) {
        foreach my $pidl (@{$hash->{$k}->{Protein_id}}) {
          print $acefh "Protein_id\t@$pidl\n"; 
        }
      }
      
      if (exists $hash->{$k}->{dblinks}) {
        my $h = $hash->{$k}->{dblinks};
        foreach my $db (sort keys %$h) {
          foreach my $attr_k (sort keys %{$h->{$db}}) {
            foreach my $attr_v (sort keys %{$h->{$db}->{$attr_k}}) {
              print $acefh "Database $db $attr_k $attr_v\n";         
            }
          }
        }
      }
      print $acefh "\n";
    }
  }

  close($acefh) or $log->log_and_die("Could not close $pidacefile properly\n");
  
  unless ($noload) {
    $wb->load_to_database($wb->autoace, $pidacefile, 'ENA_protein_xrefs', $log);
  }

}

$log->mail();
exit(0);

#########################
sub lookup_from_ebi_production_dbs {
  my ($output_file, $type) = @_;

  my $ebi_prod_dir = $wb->wormpub . "/ebi_resources";
  my $ena_perl     = "$ebi_prod_dir/ena_perl";
  my $ena_env      = "$ebi_prod_dir/ena_oracle_setup.sh";
  my $ena_cred     = "$ebi_prod_dir/ENAORACLE.s";
  my $uni_cred     = "$ebi_prod_dir/UNIPROTORACLE.s";

  if ($type eq 'proteinxrefs') {  
    my $cmd =  "source $ena_env &&"
        . " $ena_perl  $ENV{CVS_DIR}/get_protein_ids_ebiprod.pl"
        . "  -enacred $ena_cred" 
        . "  -uniprotcred $uni_cred"
        . "  -orgid $ncbi_tax_id"
        . "  -bioprojectid $bioproject_id";
    
    system("$cmd > $output_file") 
        and $log->log_and_die("Could not successfully run '$cmd'\n");

  } elsif ($type eq 'seqversions') {

    my $cmd =  "source $ena_env &&"
        . " $ena_perl  $ENV{CVS_DIR}/get_sequence_versions_ebiprod.pl"
        . "  -enacred $ena_cred"
        . "  -orgid $ncbi_tax_id"
        . "  -bioprojectid $bioproject_id";
    
    system("$cmd > $output_file") 
        and $log->log_and_die("Could not successfully run '$cmd'\n");
  } elsif ($type eq 'genexrefs') {

    my $cmd =  "source $ena_env &&"
        . " $ena_perl  $ENV{CVS_DIR}/get_gene_ids_ebiprod.pl"
        . "  -enacred $ena_cred"
        . "  -orgid $ncbi_tax_id"
        . "  -bioprojectid $bioproject_id";
    
    system("$cmd > $output_file") 
        and $log->log_and_die("Could not successfully run '$cmd'\n");

  }

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
