#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Getopt::Long;

use WormBase2Ensembl;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Utils::IO::GTFSerializer;

my ($dbname, $dbuser, $dbpass, $dbport, $dbhost, $debug,
    @dump_slice, $out_file, $out_fh);

GetOptions(
  'host=s'       => \$dbhost,
  'port=s'       => \$dbport,
  'user=s'       => \$dbuser,
  'pass=s'       => \$dbpass,
  'dbname=s'     => \$dbname,
  'slice=s@'     => \@dump_slice,
  'outfile:s'    => \$out_file,
  'debug'        => \$debug,
          )or die ("Couldn't get options");

my $ensdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
	-dbname  => $dbname,
	-host    => $dbhost,
	-port    => $dbport,
	-user    => $dbuser,
	-pass    => $dbpass,
);


my $sa = $ensdb->get_SliceAdaptor();
my $mc = $ensdb->get_MetaContainer();

my $gb_version = $mc->single_value_by_key('genebuild.version');

my (@slices, %genes_by_slice);
@dump_slice = split(/,/,join(',',@dump_slice));

if (@dump_slice) {
  foreach (@dump_slice) {
    my $sl = $sa->fetch_by_region('toplevel',$_);
    if (not defined $sl) {
      die "Could not fetch slice for $_\n";
    }
    push @slices, $sl;
  }
}
else {
  @slices = sort { $b->length <=> $a->length } @{$sa->fetch_all('toplevel')};
}


if ($out_file) {
  open( $out_fh, ">$out_file") or die "Could not open $out_file for writing\n";
} else {
  $out_fh = \*STDOUT;
}

print $out_fh "#!genebuild-version $gb_version\n";
my $serializer = Bio::EnsEMBL::Utils::IO::GTFSerializer->new($out_fh);

foreach my $slice (@slices) {
  foreach my $g (sort { $a->start <=> $b->start } @{$slice->get_all_Genes}) {
    $serializer->print_Gene($g);
  }
}

exit(0);

