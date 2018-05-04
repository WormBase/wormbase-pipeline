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
  $outfile = "./wormbase.agr_allele.${ws_version}.json";
}

#
# restrict to genes in the BGI, if provded
#
my (%only_these_genes, @alleles);

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

my ( $it, @annots);


$it = $db->fetch_many(-query => 'find Variation WHERE Live AND COUNT(Gene) == 1 AND (Description OR Interactor OR Disease_info) AND NOT Natural_variant');

while (my $obj = $it->next) {
  next unless $obj->isObject();
  #next unless $obj->Species;
  #next unless $obj->Species->name eq $full_name;

  my ($has_disease, $has_interaction, $has_phenotype) = (0,0,0);

  $has_disease = 1 if $obj->Disease_info;
  $has_phenotype = 1 if $obj->Description;
  if ($obj->Interactor) {
    foreach my $item ($obj->Interactor->col) {
      $has_interaction = 1 if $item eq 'Regulatory';
    }
  }
  
  next unless $has_disease or $has_phenotype or $has_interaction;

  my ($gene) = $obj->Gene->name;
  next if not exists $only_these_genes{$gene};

  my $symbol = $obj->Public_name->name;
  my %synonyms;
  foreach my $on ($obj->Other_name) {
    $synonyms{$on->name} = 1;
  }

  my $json_obj = {
    primaryId     => "WBVar:$obj", 
    symbol        => $symbol,
    symbolText    => $symbol,
    synonyms      => [keys \%synonyms],
    secondaryIds => [],
    taxonId       => "NCBITaxon:" . $taxid,
    gene          => "WB:$gene",
    crossReferences => [ { id => "WBVar:$obj", pages => ["allele"] }],
  };


  push @alleles, $json_obj;
}


my $meta_data = {
  dateProduced => &get_rfc_date(),
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
  data     => \@alleles,
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

