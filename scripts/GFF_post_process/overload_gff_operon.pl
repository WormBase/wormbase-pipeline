#!/usr/bin/env perl
#
# overload_gff_operon.pl
#
# Overloads operon lines with Gene info
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-07-21 11:07:59 $

use lib $ENV{CVS_DIR};
use Wormbase;
use Storable;
use IO::File;
use Getopt::Long;
use strict;
use Ace;

my ($debug,$test,$species,$store,$wormbase);
my ($gff3,$infile,$outfile, $changed_lines);

GetOptions(
  'debug=s'       => \$debug,
  'test'          => \$test,
  'species:s'     => \$species,
  'store:s'       => \$store,
  'infile:s'      => \$infile,
  'outfile:s'     => \$outfile,
  'gff3'          => \$gff3,
    );

if ($store) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

# fill some reference hashes
my ($operon_genes) = &get_operon_data();
my $prefix_name = $wormbase->pep_prefix; # 'CE' for elegans, 'BM' for brugia

while (<$gff_in_fh>){
  unless(/operon\s+operon/){
    print $gff_out_fh $_;
    next;
  }
  chomp;
  print $gff_out_fh $_;
  my ($operond) = /(${prefix_name}OP\d+|${prefix_name}OP\w+\d+)/;
  if (/Gene/) {
    $log->write_to("Error:$_\n");
    $log->log_and_die("\nIt appears that you have already overloaded the ${prefix_name}OP lines for $species\n\n");
  }
  foreach my $genes (@{$$operon_genes{$operond}}) { 
    if (defined $$operon_genes{$operond}) {
      print $gff_out_fh " ; Gene \"$genes\"";
      $changed_lines++;
    }
  }
  print $gff_out_fh "\n";
}
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
$log->write_to("Finished processing : $changed_lines lines modified\n");
$log->mail();
exit(0);


# populates Operon::Gene hash
sub get_operon_data {
  
  my %operon_genes;
  my $db = Ace->connect(-path => $wormbase->autoace);
  my $cursor = $db->fetch_many(Operon => '*');
  while (my $operon = $cursor->next){
    my @tmp_genes = $operon->Contains_gene;
    my @tmp_genes2 = map{"$_"} @tmp_genes;
    if ($debug) {print "$operon\n" unless (defined $operon->Contains_gene->name);}
    push @{$operon_genes{"$operon"}}, @tmp_genes2 if (defined $operon->Contains_gene);
  }
  $db->close;
  return \%operon_genes;
}

1;
