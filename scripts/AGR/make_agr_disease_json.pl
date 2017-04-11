#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;

use lib $ENV{CVS_DIR};
use Wormbase;


my ($debug, $test, $verbose, $store, $wormbase);
my ($outfile, $acedbpath, $ws_version, $outfh);

GetOptions (
  "debug=s"     => \$debug,
  "test"        => \$test,
  "verbose"     => \$verbose,
  "store:s"     => \$store,
  "database:s"  => \$acedbpath,
  "outfile:s"   => \$outfile,
  "wsversion=s" => \$ws_version,
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

my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or die("Connection failure: ". Ace->error);

my ( $it, @annots);

$it = $db->fetch_many(-query=>'find Gene Disease_info');

while (my $obj=$it->next) {
  next unless $obj->isObject();
  next unless $obj->Species;
  next unless $obj->Species->name eq $full_name;

  my $g = $obj->name;


  foreach my $doterm ($obj->Experimental_model) {
    my (@papers, $evi_date);
    foreach my $evi ($doterm->right->col) {
      if ($evi->name eq 'Paper_evidence') {
        foreach my $paper ($evi->col) {
          my $pmid;
          foreach my $db ($paper->Database) {
            if ($db->name eq 'MEDLINE') {
              $pmid = $db->right->right->name;
              last;
            }
          }
          push @papers, {
            publicationModId => "$paper",
          };
          if ($pmid) {
            $papers[-1]->{pubMedId} = $pmid;
          }
        }
      } elsif ($evi->name eq 'Date_last_updated') {
        $evi_date = $evi->right;
        $evi_date->date_style('ace');
        my ($y, $m, $d) = split(/\-/, $evi_date);
        $evi_date = sprintf("%4d-%02d-%02dT00:00:00+00:00", $y, $m, $d);
      }
    }
    
    push @annots, {
      objectId => "WB:$g",
      DOid     => $doterm->name,
      taxonId  => $taxid,
      dataProvider => "WB", 
      dateAssigned => defined $evi_date ? $evi_date : $date,
      geneticSex   => "hermaphrodite",
      evidence     => [
        {
          evidenceCode => "ECO:0000015",  # inferred from mutant phenotype; hard-coded for now
          publications => \@papers,
        },
          ],
      objectRelation => {
        associationType => "causes_condition",
        objectType      => "gene",
      },
    };
  }
}

$db->close;


my $meta_data = {
  dateProduced => $date,
  dataProvider => "WB", 
  release      => (defined $ws_version) ? $ws_version : $wormbase->get_wormbase_version_name(),
};

my $data = {
  metaData => $meta_data,
  data => \@annots,
};

if ($outfile) {
  open($outfh, ">$outfile") or die("cannot open $outfile : $!\n");  
} else {
  $outfh = \*STDOUT;
}

my $json_obj = JSON->new;
my $string = $json_obj->allow_nonref->canonical->pretty->encode($data);
print $outfh $string;


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

