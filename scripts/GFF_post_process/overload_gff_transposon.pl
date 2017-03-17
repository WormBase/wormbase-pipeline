#!/usr/bin/env perl
#
# overload_gff_transposon.pl
#
# Overloads Transpososn lines with Family info
#

use lib $ENV{CVS_DIR};
use Wormbase;
use Storable;
use IO::File;
use Getopt::Long;
use strict;
use Ace;
 

my ($debug,$test,$species,$store,$wormbase,$database);
my ($infile,$outfile, $gff3, $changed_lines);

GetOptions(
  'debug=s'       => \$debug,
  'test'          => \$test,
  'species:s'     => \$species,
  'store:s'       => \$store,
  'infile:s'      => \$infile,
  'outfile:s'     => \$outfile,
  'gff3'          => \$gff3,
  'database:s'    => \$database,
    );

if ($store) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
			     );
}

$database = $wormbase->autoace if not defined $database;
my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

$changed_lines = 0;

# fill some reference hashes
my ($family_names) = &get_family_name_data();
while (<$gff_in_fh>){
  unless(/WormBase_transposon\s+transposable_element|Transposon\s+transposable_element/){
    print $gff_out_fh $_;
    next;
  }
  chomp;

  print $gff_out_fh $_;
  my ($transposond) = /(WBTransposon\d+|Predicted_\w+\d+)/;

  if (exists $family_names->{$transposond}) {
    my $gfamily=$family_names->{$transposond}[0];
    if ($gff3) {
      print $gff_out_fh ";Family=$gfamily\n";
    } else {
      print $gff_out_fh " ; Family \"$gfamily\"\n";
    }
    $changed_lines++;
  }
  else {
    print "Warning: $transposond does not have an assigned Transposon family\n";
    print $gff_out_fh "\n";
  }
}
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
$log->write_to("Finished processing : $changed_lines lines modified\n");
$log->mail();
exit(0);


# populates Operon::Gene hash
sub get_family_name_data {
  my %family_names;
  my $db = Ace->connect(-path => $database);
  my $cursor = $db->fetch_many(-query => 'Find Transposon where Member_of');
  while (my $Transposon = $cursor->next){
    foreach my $family ($Transposon->Member_of) {
      push @{$family_names{"$Transposon"}}, $family->name;
    }
  }
  $db->close;
  return \%family_names;
}

1;
