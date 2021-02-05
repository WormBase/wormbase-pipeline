#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;


use lib $ENV{CVS_DIR};
use Wormbase;
use Log_files;

use lib "$ENV{CVS_DIR}/ONTOLOGY";
use GAF;

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($outfile, $acedbpath);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "database:s" => \$acedbpath,
	    "output:s"   => \$outfile,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
      );
}

my $tace = $wormbase->tace;
my $log  = Log_files->make_build_log($wormbase);
my $date = &get_GAF_date();
my $taxid = $wormbase->ncbi_tax_id;
my $full_name = $wormbase->full_name;

$acedbpath = $wormbase->autoace unless $acedbpath;
$outfile = $wormbase->ontology."/disease_association.".$wormbase->get_wormbase_version_name.".wb" unless $outfile;

$log->write_to("Connecting to database $acedbpath\n");
my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or $log->log_and_die("Connection failure: ". Ace->error);

my ($gene_info, $count, $it);

$gene_info = &get_gene_info( $acedbpath, $wormbase, $full_name );

$log->write_to( "Got name information for " . scalar(keys %$gene_info) . " genes\n");

open(my $out, ">$outfile") or $log->log_and_die("cannot open $outfile : $!\n");

&print_wormbase_GAF_header($out, $wormbase->get_wormbase_version_name);

$it = $db->fetch_many(-query=>'find Gene Disease_info');

while (my $obj=$it->next) {
  next unless $obj->isObject();
  next unless $obj->Species;
  next unless $obj->Species->name eq $full_name;

  my $g = $obj->name;

  foreach my $doterm ($obj->Potential_model) {

    my $meth = $doterm->right->right;
    if ($meth->name eq 'Inferred_automatically') {
      my $text = $meth->right->name;
      my ($with_from_list) = $text =~ /\((\S+)\)/;
      
      my @ens = map { "ENSEMBL:$_" } grep { $_ =~ /ENSG\d+/ } split(/,/, $with_from_list);
      my @omim = grep { $_ =~ /OMIM:/ } split(/,/, $with_from_list);

      &print_wormbase_GAF_line($out,  
                               $g, 
                               $gene_info->{$g}->{public_name}, 
                               "",  
                               $doterm,
                               "PMID:19029536",  # this is the reference for Ensembl Compara
                               "IEA", 
                               join("|", @ens,@omim),
                               "D",
                               $gene_info->{$g}->{sequence_name},
                               $taxid, 
                               $date);
    }
  }

  foreach my $doterm ($obj->Experimental_model) {
    my (@papers, @omims);
    foreach my $evi ($doterm->right->col) {
      if ($evi->name eq 'Paper_evidence') {
        foreach my $paper ($evi->col) {
          push @papers, $paper;
        }
      } elsif ($evi->name eq 'Accession_evidence') {
        foreach my $db ($evi->col) {
          if ($db->name eq 'OMIM') {
            foreach my $acc ($db->col) {
              push @omims, "OMIM:$acc";
            }
          }
        }
      }
    }
    foreach my $paper (@papers) {
      &print_wormbase_GAF_line($out,  
                               $g, 
                               $gene_info->{$g}->{public_name}, 
                               "",  
                               $doterm,
                               "WB_REF:$paper",
                               "IMP", 
                               join("|", @omims),
                               "D",
                               $gene_info->{$g}->{sequence_name},
                               $taxid, 
                               $date);
    }
  }
}

$db->close;
$log->mail;

exit(0);

	

