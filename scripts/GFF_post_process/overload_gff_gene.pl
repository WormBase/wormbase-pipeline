#!/usr/bin/env perl
#
# overload_gff_gene.pl
#
# Adds interpolated map positions and other information to gene and allele lines
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-09-25 14:09:20 $


use strict;                                      
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Ace;
use File::Copy;

my ($help, $debug, $test, $store, $wormbase, $species);
my ( $gff3, $infile, $outfile, $changed_lines);

GetOptions (
  "debug=s"   => \$debug,
  "test"      => \$test,
  "store:s"   => \$store,
  "species:s" => \$species,
  "gff3"      => \$gff3,
  "infile:s"  => \$infile,
  "outfile:s" => \$outfile,
    );

if ( $store ) {
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
my $species_full = $wormbase->full_name;

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

my ($locus_h, $sequence_name_h, $biotype_h) = &get_data();

while (<$gff_in_fh>) {
  if (/^\#/) {
    print $gff_out_fh $_;
    next;
  }
  chomp;
  my @f = split(/\t/, $_);

  my $orig_attr = $f[8];

  if ($f[2] eq 'gene' and ($f[1] eq 'gene' or $f[1] eq 'WormBase')) {
    my $changed = 0;
    if ($gff3) {
      my ($first) = split(/;/, $f[8]);
      my ($gene) = $first =~ /ID=Gene:(\S+)/;

      $f[8] .= ";locus=$locus_h->{$gene}" 
          if exists $locus_h->{$gene} and $f[8] !~ /locus/;

      $f[8] .= ";sequence_name=$sequence_name_h->{$gene}" 
          if exists $sequence_name_h->{$gene} and $f[8] !~ /sequence_name/;

      $f[8] .= ";biotype=$biotype_h->{$gene}" 
          if exists $biotype_h->{$gene} and $f[8] !~ /biotype/;
    } else {
      my ($gene) = $f[8] =~ /Gene\s+\"(\S+)\"/;

      $f[8] .= " ; Locus \"$locus_h->{$gene}\"" 
          if exists $locus_h->{$gene} and $f[8] !~ /Locus/;

      $f[8] .= " ; Sequence_name \"$sequence_name_h->{$gene}\"" 
          if exists $sequence_name_h->{$gene} and $f[8] !~ /Sequence_name/;

      $f[8] .= " ; Biotype\"$biotype_h->{$gene}\"" 
          if exists $biotype_h->{$gene} and $f[8] !~ /Biotype/;
    }	
  }

  $changed_lines++ if $orig_attr ne $f[8];
  print $gff_out_fh join("\t", @f), "\n";
}
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
$log->write_to("Finished processing : $changed_lines lines modified\n");
$log->mail();
exit(0);

#######################
sub get_data {
  my $db = Ace->connect('-path' => $wormbase->autoace) 
      or $log->log_and_die("cant open Ace connection to db\n".Ace->error."\n");

  my (%locus, %sequence_name, %biotype);
  
  # get the CGC/WGN name of the Genes
  $log->write_to("Reading WGN names of genes\n");
  
  my $query = "find Gene where Species = \"$species_full\"";
  my $genes = $db->fetch_many('-query' => $query);
  while (my $gene = $genes->next){
    if ($gene->CGC_name) {
      $locus{$gene->name} = $gene->CGC_name;
    }
    if ($gene->Sequence_name) {
      $sequence_name{$gene->name} = $gene->Sequence_name;
    }
  }

  my (%gene_biotypes);

  $log->write_to("Fetching transcript info\n");

  my $query = "find Transcript where Gene";
  my $trans = $db->fetch_many('-query' => $query);
  while(my $tran = $trans->next) {
    my $gene = $tran->Gene;
    next if not $gene;
    if ($tran->Corresponding_CDS) {
      $gene_biotypes{$gene}->{protein_coding}++;
    } elsif ($tran->Transcript) {
      $gene_biotypes{$gene}->{$tran->Transcript}++;
    }
  }
  
  $query = "find Pseudogene where Gene";
  my $pseudos = $db->fetch_many('-query' => $query);
  while(my $pseudo = $pseudos->next) {
    my $gene = $pseudo->Gene;
    $gene_biotypes{$gene}->{pseudogene}++;
  }

  $db->close();
  $log->write_to("Closed connection to DB\n");

  # a gene may have more than one kind of biotype, so rule as follows:
  # at least one protein-coding => protein_coding
  # at least one kind of RNA => 
  #   different kinds => ncRNA
  # only pseudogene => pseudogene
  foreach my $g (keys %gene_biotypes) {
    if(exists $gene_biotypes{$g}->{protein_coding}) {
      $biotype{$g} = 'protein_coding';
    } else {
      my @types = sort keys %{$gene_biotypes{$g}};
      my @rna_types = grep { /RNA/ } @types;

      if (@rna_types) {
        if (scalar(@rna_types) == 1) {
          $biotype{$g} = $rna_types[0];
        } else {
          $biotype{$g} = 'ncRNA';
        }
      } else {
        # has to be a pseudogene, right? 
        if (not exists $gene_biotypes{$g}->{pseudogene}) {
          $log->log_and_die("Gene $g could not be assigned a biotype (@types); exiting\n");
        } else {
          $biotype{$g} = 'pseudogene';
        }
      }
    }
  }
  
  return (\%locus, \%sequence_name, \%biotype);
}
