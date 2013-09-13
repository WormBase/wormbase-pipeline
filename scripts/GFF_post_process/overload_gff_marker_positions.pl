#!/usr/bin/env perl
#
# overload_gff_marker_positions.pl
#
# Adds interpolated map positions and other information to gene and allele lines
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-09-13 12:18:29 $


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

my (%variation, %gene, %gene_exact, %locus);
&get_map_data();

while (<$gff_in_fh>) {
  if (/^\#/) {
    print $gff_out_fh $_;
    next;
  }
  chomp;
  my @f = split(/\t/, $_);

  my $orig_attr = $f[8];

  if (defined $f[8] && $f[1] eq 'Allele') {
    if ($gff3) {
      my ($first) = split(/;/, $f[8]);
      my ($allele) = $first =~ /Variation:(\S+)/;
      $f[8] .= ";interpolated_map_position=$variation{$allele}"
          if exists $variation{$allele} and $f[8] !~ /Interpolated_map_position/;
    } else {
      my ($allele) = $f[8] =~ /Variation\s+\"(\S+)\"/;
      $f[8] .= " ; Interpolated_map_position \"$variation{$allele}\""
          if exists $variation{$allele} and $f[8] !~ /Interpolated_map_position/;
    }	    
  } elsif ($f[2] eq 'gene' and ($f[1] eq 'gene' or $f[1] eq 'WormBase')) {
    my $changed = 0;
    if ($gff3) {
      my ($first) = split(/;/, $f[8]);
      my ($gene) = $first =~ /ID=Gene:(\S+)/;
      $f[8] .= ";interpolated_map_position=$gene{$gene}"
          if exists $gene{$gene} and $f[8] !~ /interpolated_map_position/;
      $f[8] .= ";position=$gene_exact{$gene}" 
          if exists $gene_exact{$gene} and $f[8] !~ /position/;
      $f[8] .= ";locus=$locus{$gene}" 
          if exists $locus{$gene} and $f[8] !~ /locus/;
    } else {
      my ($gene) = $f[8] =~ /Gene\s+\"(\S+)\"/;
      $f[8] .= " ; Interpolated_map_position \"$gene{$gene}\"" 
          if exists $gene{$gene} and $f[8] !~ /Interpolated_map_position/;
      $f[8] .= " ; Position \"$gene_exact{$gene}\"" 
          if exists $gene_exact{$gene} and $f[8] !~ /Position/;
      $f[8] .= " ; Locus \"$locus{$gene}\"" 
          if exists $locus{$gene} and $f[8] !~ /Locus/;
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
sub get_map_data {
  my $db = Ace->connect('-path' => $wormbase->autoace) 
      or $log->log_and_die("cant open Ace connection to db\n".Ace->error."\n");
  
  $log->write_to("Getting interpolated map position for variations\n");
  
  my $tmfh = $wormbase->table_maker_query($wormbase->autoace, &write_interpolated_map_position_query());
  while (<$tmfh>) {
    chomp;
    s/\"//g; 
    next if (/acedb/ or /\/\//);
    next if /^\s*$/;
    
    my @data = split(/\t/, $_);
    $variation{$data[0]} = $data[2];
  }
  
  # get the interpolated physical mapping position of the Genes
  $log->write_to("Reading interpolated Gene mapping data\n");
  
  my $query = "find Gene where Interpolated_map_position";
  my $genes = $db->fetch_many('-query' => $query);
  while (my $gene = $genes->next){
    $gene{$gene->name} = $gene->Interpolated_map_position(2);
  }
  
# get the exact physical mapping position of the Genes
  $log->write_to("Reading exact Gene mapping data\n");
  
  $query = "find Gene where Species = \"'$species_full\"' AND Map AND NEXT AND NEXT = Position";
  $genes = $db->fetch_many('-query' => $query);
  while (my $gene = $genes->next){
    $gene_exact{$gene->name} = $gene->Map(3);
  }
  
  # get the CGC/WGN name of the Genes
  $log->write_to("Reading WGN names of genes\n");
  
  $query = "find Gene where Species = \"$species_full\"";
  $genes = $db->fetch_many('-query' => $query);
  while (my $gene = $genes->next){
    if ($gene->CGC_name) {
      $locus{$gene->name} = $gene->CGC_name;
    } elsif ($gene->Sequence_name) {
      $locus{$gene->name} = $gene->Sequence_name;
    }
  }
  
  $db->close();
  $log->write_to("Closed connection to DB\n");
}


#######################
sub write_interpolated_map_position_query {

  my $tmp_file = "/tmp/interpolated_map.tmquery.txt";
  open(my $qfh, ">$tmp_file") or $log->log_and_die("Could not open TM query file $tmp_file\n");

  my $query = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Variation 
From 1 
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class Map 
From 1 
Tag Interpolated_map_position 
 
Colonne 3 
Width 12 
Mandatory 
Visible 
Float 
Right_of 2 
Tag  HERE  

EOF

  print $qfh $query;
  close($qfh);

  return $tmp_file;
}
