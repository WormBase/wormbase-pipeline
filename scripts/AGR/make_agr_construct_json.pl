#!/usr/bin/env perl
# create a JSON dump of the construct class for AGR
#
# WARNING: the JSON schema will change and the script will need rewriting
# Schema : https://github.com/alliance-genome/agr_schemas/tree/data_dictionary/ingest/transgenicConstruct

use strict;
use Ace;
use Modules::AGR;
use JSON;
use Getopt::Long;
use Wormbase;

my ($acedbpath,$ws_version,$outfile,$debug,$test,$verbose,$store,$wormbase);

GetOptions (
  "debug=s"     => \$debug,
  "test"        => \$test,
  "verbose"     => \$verbose,
  "store:s"     => \$store,
  "database:s"  => \$acedbpath,
  "outfile:s"   => \$outfile,
  "wsversion=s" => \$ws_version,
)||die(@!);

if ($store) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug => $debug, -test => $test,);
}

$acedbpath||=$wormbase->autoace;
my $date = AGR::get_rfc_date();

my $db = Ace->connect(-path => $acedbpath) or die("Connection failure: ". Ace->error);
my $it = $db->fetch_many(-query => 'find Construct WHERE (Public_name OR Summary OR Driven_by_gene OR Gene)');

my @constructs;
while (my $obj = $it->next) {
     my %json_obj;

     $json_obj{primaryID} = "WB:$obj";
     $json_obj{symbol} = "${\$obj->Public_name}" if $obj->Public_name;
     $json_obj{synonyms} = ["${\$obj->Summary}"] if $obj->Summary;

     # the two below are actually lists in the database, so it needs changing the JSON schema
     $json_obj{constructComponents}=[];
     push @{$json_obj{constructComponents}},{ componentRelation => 'is_regulated_by',componentID => "WB:${\$obj->Driven_by_gene}", componentSymbol => "${\$obj->Driven_by_gene->Public_name}"} if $obj->Driven_by_gene;
     push @{$json_obj{constructComponents}},{ componentRelation => 'expresses',componentID => "WB:${\$obj->Gene}", componentSymbol => "${\$obj->Gene->Public_name}"} if $obj->Gene;

     push @constructs, \%json_obj;
}

my $data = {
  metaData => AGR::get_file_metadata_json( (defined $ws_version) ? $ws_version : $wormbase->get_wormbase_version_name(), $date ),
  data     => \@constructs,
};

my $out_fh;
if ($outfile) {
  open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";
} else {
  $out_fh = \*STDOUT;
}


my $json_obj = JSON->new;
my $string = $json_obj->allow_nonref->canonical->pretty->encode($data);
print $out_fh $string;

$db->close;
exit(0);

