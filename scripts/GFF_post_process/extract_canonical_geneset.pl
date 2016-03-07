#!/usr/bin/env perl
#
# extract_canonical_geneset.pl
#
# Extracts the canonical gene set from the final GFF3 files, and writes to another
# GFF3 file by default, or optionally a GTF file.
# (a flavout of GFF2 required by many software packages)
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2015-02-23 17:20:49 $

use strict;
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;


######################################
# variables and command-line options #
######################################

my ( $infile, $outfile, %trans2gene, $gtf );

GetOptions(
  'infile:s'  => \$infile,
  'outfile:s' => \$outfile,
  'gtf'       => \$gtf,
);


if (not defined $infile or not defined $outfile) { 
  die("You must define -infile and -outfile\n");
}

my ($gff_in_fh, $gff_out_fh);
if ($infile =~ /\.gz$/) {
  open( $gff_in_fh, "gunzip -c $infile |" ) 
      or die("Could not open gzip stream to $infile\n");
} else {
  open($gff_in_fh, $infile) or die("Could not open $infile for reading\n");
}

if ($outfile =~ /\.gz/) {
  open($gff_out_fh, "| gzip -n -c > $outfile") 
      or die("Could not open gzip write stream to $outfile\n");
} else {
  open($gff_out_fh, ">$outfile") 
      or die("Could not open $outfile for writing\n");
}

if ($gtf) {
  print $gff_out_fh "\#\#gff-version 2\n"; 
}

while(<$gff_in_fh>) {
  /^\#/ and do {
    print $gff_out_fh $_ unless $gtf;    
    next;
  };

  chomp; 

  my @l = split(/\t/, $_);

  next if ($l[1] ne 'WormBase' and $l[1] ne 'WormBase_imported');

  if (not $gtf) {
    print $gff_out_fh join("\t", @l), "\n";
    next;
  }

  my ($id)   = $l[8] =~ /ID=([^;]+)/;
  $id =~ s/\S+://;

  if ($l[2] eq 'gene') {
    #
    # gene
    #
    my ($biotype) = $l[8] =~ /biotype=([^;]+)/;
    my ($locus)   = $l[8] =~ /locus=([^;]+)/;

    $l[8] = "gene_id \"$id\"";

    if (defined $biotype) {
      $l[8] .= "; gene_biotype \"$biotype\""; 
    }
    if (defined $locus) {
      $l[8] .= " ; locus \"$locus\"";
    }

    print $gff_out_fh join("\t", @l), "\n";

  } elsif ($l[2] =~ /RNA/ or $l[2] =~ /pseudo/i or $l[2] =~ /nc_primary_transcript/) {
    #
    # transcript
    #
    my ($gene) = $l[8] =~ /Parent=([^;]+)/;
    $gene =~ s/\S+://; 
    $trans2gene{$id} = $gene;
    my $biotype = $l[2];
    $biotype = "protein_coding" if $biotype eq 'mRNA';

    $l[2] = "transcript";
    $l[8] = "gene_id \"$gene\" ; transcript_id \"$id\" ; transcript_biotype \"$biotype\""; 

    print $gff_out_fh join("\t", @l), "\n";    

  } elsif ($l[2] eq 'exon') {
    #
    # exon
    #    
    my ($transcript) = $l[8] =~ /Parent=([^;]+)/;
    $transcript =~ s/\S+://; 
    my $gene = $trans2gene{$transcript};

    $l[8] = "gene_id \"$gene\" ; transcript_id \"$transcript\"" ;

    print $gff_out_fh join("\t", @l), "\n";
  } elsif ($l[2] eq 'CDS') {
    #
    # CDS
    #
    my ($tran_list) = $l[8] =~ /Parent=([^;]+)/;
    my @transcripts = map { s/\S+://; $_ } split(/,/, $tran_list);

    foreach my $tran (@transcripts) {
      my $gene = $trans2gene{$tran};
      $l[8] = "gene_id \"$gene\" ; transcript_id \"$tran\"";

      print $gff_out_fh join("\t", @l), "\n";
    }
  }
} 
close($gff_out_fh) or die("Could not close $outfile after writing\n");
exit(0);
