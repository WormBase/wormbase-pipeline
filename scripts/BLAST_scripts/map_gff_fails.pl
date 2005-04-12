#!/usr/local/perl

use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Coords_converter;
use Wormbase;

my $database = shift;
$database = "/wormsrv2/autoace" unless $database;

#CHROMOSOME_I   waba_coding     similarity      7504218 7504560 378.0000        -       .       Target "Sequence:cb25.fpc2695" 196866 197208 197208,7504218,0

my $coords = Coords_converter->invoke("$database");
#my $coords = Coords_converter->invoke();

while(<>) {
  chomp;
  my @data = split(/\t/,$_);
  my @clone_coords = $coords->LocateSpan($data[0], $data[3], $data[4]);

  $data[0] = $clone_coords[0];
  $data[3] = $clone_coords[1];
  $data[4] = $clone_coords[2];

  print join("\t",@data),"\n";
}
