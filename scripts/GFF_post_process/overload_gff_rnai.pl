#!/usr/bin/env perl
#
# overload_gff_rnai.pl
#
# Overloads RNAi mapping lines with extra info (Lab of clone etc)
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2014-08-27 21:50:10 $      

use lib $ENV{CVS_DIR};
use Wormbase;
use Storable;
use IO::File;
use Getopt::Long;
use strict;
use Ace;

my ($debug,$test,$species,$store,$wormbase,$database);
my ($infile,$outfile,$gff3, $changed_lines);

GetOptions(
  'debug=s'    => \$debug,
  'test'       => \$test,
  'species:s'  => \$species,
  'store:s'    => \$store,
  'infile:s'   => \$infile,
  'outfile:s'  => \$outfile,
  'database:s' => \$database,
  'gff3'       => \$gff3,
)||die(@!);

if ($store) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
			     );
}

$database = $wormbase->autoace if not defined $database;

# establish log file.
my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

my ($r2lab,$r2hist) = &get_rnai2lab();

while ( <$gff_in_fh> ) {
  unless(/RNAi_(primary|secondary)\s+RNAi_reagent/){
    print $gff_out_fh $_;
    next;
  }
  chomp;
  my @l = split(/\t/, $_);

  my ($rnaid) = /(WBRNAi\d+)/;
  my $old_attr = $l[8];
  if ($gff3) {
    $l[8] .= ";laboratory=$$r2lab{$rnaid}" if $$r2lab{$rnaid};
    $l[8] .= ";history_name=$$r2hist{$rnaid}" if $$r2hist{$rnaid};
  } else {
    $l[8] .= " ; Laboratory \"$$r2lab{$rnaid}\"" if $$r2lab{$rnaid};
    $l[8] .= " ; History_name \"$$r2hist{$rnaid}\"" if $$r2hist{$rnaid};
  }
  $changed_lines++ if $l[8] ne $old_attr;
  print $gff_out_fh join("\t", @l), "\n";
}
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
$log->write_to("Finished processing : $changed_lines lines modified\n");
$log->mail();
exit(0);


#################################################
sub get_rnai2lab {
  my %rnai2lab;
  my %rnai2history;
  my $db = Ace->connect(-path => $database);
  my $cursor = $db->fetch_many(RNAi => '*');
  while (my $rnai = $cursor->next){
    $rnai2lab{"$rnai"}="${\$rnai->Laboratory}" if $rnai->Laboratory;
    $rnai2history{"$rnai"}="${\$rnai->History_name}" if $rnai->History_name;
  }
  $db->close;
  return \%rnai2lab,\%rnai2history;
}

1;
