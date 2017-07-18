#!/usr/bin/env perl

use strict;

use IO::File;
use Storable;
use Getopt::Long;

use lib $ENV{CVS_DIR};

use Wormbase;
use Log_files;

my ($store,$debug,$test,$database,$species,$outfile);
GetOptions(
  'store=s' => \$store,
  'debug=s' => \$debug,
  'test'    => \$test,
  'species=s'  => \$species,
  'database=s' => \$database,
  'outfile=s'  => \$outfile,
    )||die(@!);

my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n");
} 
else {
  $wormbase = Wormbase->new( -debug => $debug, -test => $test,-autoace=> $database,-organism => $species);
}

my $log = Log_files->make_build_log($wormbase);

$database = $wormbase->autoace if not defined $database;
    
$outfile = $wormbase->reports . '/pcr_product2gene.txt'
    if not defined $outfile;
my $of = IO::File->new($outfile,'w');
$log->write_to("writing to $outfile\n");

my $tace    = $wormbase->tace;
my $command = "Table-maker -p $database/wquery/pcr_product2gene.def\nquit\n";

# hashes needed because one pcr product may hit two or more genes
my %pcr2gene;
my %gene2sequence;
my %gene2cgc;

open (my $tace_fh, "echo '$command' | $tace $database | ") or $log->log_and_die("Could nt tace that $database\n");
while (<$tace_fh>){
  if (m/^\"/){
    s/\"//g;
    chomp;
    # split the line into various fields
    my ($pcr,$gene,$cgc,$sequence) = split(/\t/, $_) ;
    
    # fill hashes and arrays, allow for multiple genes per PCR product
    $pcr2gene{$pcr}      .= "$gene,";
    ($gene2cgc{$gene} = $cgc) if($cgc);
    $gene2sequence{$gene} = "$sequence";
  }
}
close($tace_fh);

foreach my $pcr (keys %pcr2gene){
  my @genes = split(/,/,$pcr2gene{$pcr});
  my $counter =0;
  print $of "$pcr";
  foreach my $gene (@genes){
    
    # remove next element if it is the same gene ID and start loop again
    if (defined($genes[$counter+1]) && $genes[$counter+1] eq $genes[$counter]){
      splice(@genes, $counter+1,1);
      redo;
    }
    $counter++;
    print $of "\t$gene";
    print $of "($gene2cgc{$gene})" if (exists($gene2cgc{$gene}));
    
    # now print sequence name
    print $of ",$gene2sequence{$gene}";      
  }
  print $of "\n";
}
close($of);

$of->close;
$log->mail;
exit(0);
