#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;

use lib $ENV{CVS_DIR};
use Wormbase;


my ($debug, $test, $verbose, $store, $wormbase);
my ($outfile, $acedbpath, $ws_version, $out_fh, $bgi_json);

GetOptions (
  "debug=s"     => \$debug,
  "test"        => \$test,
  "verbose"     => \$verbose,
  "store:s"     => \$store,
  "database:s"  => \$acedbpath,
  "outfile:s"   => \$outfile,
  "wsversion=s" => \$ws_version,
  "bgijson=s"   => \$bgi_json,
    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
      );
}

my $tace = $wormbase->tace;
my $date = &get_rfc_date();
my $alt_date = join("/", $date =~ /^(\d{4})(\d{2})(\d{2})/);
my $taxid = $wormbase->ncbi_tax_id;
my $full_name = $wormbase->full_name;

$acedbpath = $wormbase->autoace unless $acedbpath;
$ws_version = $wormbase->get_wormbase_version_name unless $ws_version;

if (not defined $outfile) {
  $outfile = "./wormbase.agr_phenotype.${ws_version}.json";
}

#
# restrict to genes in the BGI, if provded
#
my (%only_these_genes, @pheno_annots);

if (defined $bgi_json) {
  # read it so that we can restrict to genes that have basic info
  my %genes;

  my $json_string = "";
  open(my $json_fh, $bgi_json) or die "Could not open $bgi_json for reading\n";
  while(<$json_fh>) {
    $json_string .= $_;
  }

  my $json_reader = JSON->new;
  my $json = $json_reader->decode($json_string);
  
  foreach my $entry (@{$json->{data}}) {
    my $id = $entry->{primaryId};
    $id =~ s/WB://; 
    $only_these_genes{$id} = 1;
  }

}


my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or die("Connection failure: ". Ace->error);

my $it = $db->fetch_many(-query => 'find Variation WHERE Live AND COUNT(Gene) == 1 AND Phenotype AND NOT Natural_variant');

while (my $obj = $it->next) {
  next unless $obj->isObject();

  my ($gene) = $obj->Gene->name;
  next if defined $bgi_json and not exists $only_these_genes{$gene};
  
  foreach my $pt ($obj->Phenotype) {
    my $phen_id = $pt->name;
    my $phen_desc = $pt->Primary_name->name;
    my @paper;
    
    foreach my $evi ($pt->col()) {
      if ($evi->name eq 'Paper_evidence') {
        foreach my $wb_paper ($evi->col ) {
          push @paper, &get_paper( $wb_paper );
        }
      }
    }
    
    foreach my $pap (@paper) {
      my $json_obj = {
        objectId     => "WBVar:$obj", 
        phenotypeTermIdentifiers => [ { termId => $phen_id, termOrder => 1 } ],
        phenotypeStatement => $phen_desc,
        dateAssigned => $date,
      };
      foreach my $k (keys %$pap) {
        $json_obj->{$k} = $pap->{$k};
      }
      
      push @pheno_annots, $json_obj;
    }
  }
}


my $meta_data = {
  dateProduced => $date,
  dataProvider => [
    { 
      crossReference => {
        id => "WB",
        pages => ["homepage"],
      },
      type => "curated",
    },
      ],
  release      => (defined $ws_version) ? $ws_version : $wormbase->get_wormbase_version_name(),
};

my $data = {
  metaData => $meta_data,
  data     => \@pheno_annots,
};

if (defined $outfile) {
  open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";
} else {
  $out_fh = \*STDOUT;
}


my $json_obj = JSON->new;
my $string = $json_obj->allow_nonref->canonical->pretty->encode($data);
print $out_fh $string;

$db->close;

exit(0);


exit(0);


##############################################
sub get_rfc_date {

  my $date;
  
  open(my $date_fh, "date --rfc-3339=seconds |");
  while(<$date_fh>) {
    if (/^(\S+)\s+(\S+)/) {
      $date = "${1}T${2}";
    }
  }
  
  return $date;
}

##############################################
sub get_paper_json {
  my ($wb_paper) = @_;
  
  my $json_paper = {};
  
  my $pmid;
  foreach my $db ($wb_paper->Database) {
    if ($db->name eq 'MEDLINE') {
      $pmid = $db->right->right->name;
      last;
    }
  }
  $json_paper->{modPublicationId} = "WBPaper:$wb_paper";
  if ($pmid) {
    $json_paper->{pubMedId} = "PMID:$pmid";
  }
  
  return $json_paper;
}

sub get_paper {
  my ($wb_paper) = @_;
  
  my $pmid;
  foreach my $db ($wb_paper->Database) {
    if ($db->name eq 'MEDLINE') {
      $pmid = $db->right->right->name;
      last;
    }
  }
  
  my $obj = {
    pubModId => "WBPaper:$wb_paper",
  }; 
  if ($pmid) {
    $obj->{pubMedId} = "PMID:$pmid";
  }
  
  return $obj;
}
