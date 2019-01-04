#!/usr/bin/env perl
#
# overload_gff_operon.pl
#
# Overloads operon lines with Gene info
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
my ($gff3,$infile,$outfile, $changed_lines);

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

# fill some reference hashes
my ($operon_genes) = &get_operon_data();
my $prefix_name = $wormbase->pep_prefix; # 'CE' for elegans, 'BM' for brugia

while (<$gff_in_fh>){
  unless((/operon\s+operon/) || (/dicistronic_mRNA\s+operon/)){
    print $gff_out_fh $_;
    next;
  }
  chomp;

  print $gff_out_fh $_;
  my ($operond) = /(${prefix_name}OP\d+|${prefix_name}OP\w+\d+)/;

  if (exists $operon_genes->{$operond}) {
    my @g = sort grep { defined $_ } @{$operon_genes->{$operond}};

    if ($gff3) {
      my $g_str = join(",", @g);
      print $gff_out_fh ";genes=$g_str";
    } else {
      foreach my $genes (@g) {
        print $gff_out_fh " ; Gene \"$genes\"";
      }
    }
    $changed_lines++;
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
  my $db = Ace->connect(-path => $database);
  my $cursor = $db->fetch_many(-query => 'Find Operon where Canonical_parent AND Contains_gene');
  while (my $operon = $cursor->next){
    foreach my $gene ($operon->Contains_gene) {
      push @{$operon_genes{"$operon"}}, $gene->name;
    }
  }
  $db->close;
  return \%operon_genes;
}

1;
