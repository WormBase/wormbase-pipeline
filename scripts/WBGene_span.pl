#!/usr/local/bin/perl5.8.0 -w
#
# WBGene_span.pl
#
# by Anthony Rogers
#
# Creates SMapped Gene spans for Gene objects
#
# Last edited by: $Author: krb $
# Last edited on: $Date: 2004-08-06 09:05:40 $

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Coords_converter;
use Getopt::Long;
use Log_files;

my $database = "/wormsrv2/autoace";
my ($test, $gff, $no_ace, $debug);
GetOptions (
	    'database=s' => \$database,
	    'test'       => \$test,
	    'gff'        => \$gff,
	    'no_ace'     => \$no_ace,
	    'debug=s'    => \$debug
	   );

my $log = Log_files->make_build_log($debug);

$log->write_to("Generating WBGene spans from database $database\n");


my %CDS_data;
my $cds;
my %gene2seq;

my $def_file = "$database/wquery/WBGene2seq.def";
my $tace = &tace;

my $coords = Coords_converter->invoke($database);

my $command="Table-maker -p $def_file\nquit\n";

open (FH, "echo '$command'| $tace $database | ");

while (<FH>) {
  s/{\t}*/\t/g;
  s/\"//g;
  my @data = split;
  next unless $data[1];
  my $gene = shift @data;
  $gene2seq{$gene} = [@data];
}

close FH;



my @chromosomes = qw( I II III IV V X MtDNA);
@chromosomes = qw(III) if $test;
foreach my $chrom ( @chromosomes ) {
  my %gene_coords;
  open (GFF,"<$database/CHROMOSOMES/CHROMOSOME_${chrom}.gff") or die "cant open GFF $chrom\n";
  while ( <GFF> ) {
    my @data = split;
    if ( ($data[1] eq "Coding_transcript") or
	 ($data[1] eq "Non_coding_transcript") or
	 ($data[1] eq "Pseudogene") or
	 ($data[1] eq "curated") or
	 ($data[1] eq "tRNAscan-SE-1.23") or
	 ($data[1] eq "tRNA") or
	 ($data[1] eq "snRNA") or
	 ($data[1] eq "miRNA") or
	 ($data[1] eq "rRNA") or
	 ($data[1] eq "scRNA") or
	 ($data[1] eq "snoRNA") or
	 ($data[1] eq "tRNA") or
	 ($data[1] eq "stRNA") or
	 ($data[1] eq "snRNA")
       ) {
      next if ( $data[2] eq "exon" or $data[2] eq "coding_exon" or $data[2] eq "intron" );
      my ($gene) = $data[9] =~ /\"(\w+\.\w+?)/;
      $gene_coords{$gene} = [($data[3], $data[4], $data[6]) ];
    }
  }
  close GFF;
  my %gene_span;
 WBG:foreach my $WBgene ( keys %gene2seq ){
    foreach my $CDS ( @{$gene2seq{$WBgene}} ) {
      next WBG unless $gene_coords{$CDS};
      if ( !(defined $gene_span{$WBgene}) ) {
	$gene_span{$WBgene}->{'min'}    = $gene_coords{$CDS}->[0];
	$gene_span{$WBgene}->{'max'}    = $gene_coords{$CDS}->[1];
	$gene_span{$WBgene}->{'strand'} = $gene_coords{$CDS}->[2];
      } else {
	if ( $gene_coords{$CDS}->[0] < $gene_span{$WBgene}->{'min'} ) {
	  $gene_span{$WBgene}->{'min'}    = $gene_coords{$CDS}->[0];
	}
	if ( $gene_span{$WBgene}->{'max'} < $gene_coords{$CDS}->[1] ) {
	  $gene_span{$WBgene}->{'max'}    = $gene_coords{$CDS}->[1];
	}
      }
    }
  }

  if( $gff ) {
    open (OUTGFF,">$database/CHROMOSOMES/CHROMOSOME_${chrom}_WBgene.gff") or do{ $log->write_to("cant open output\n"); die "cant open output\n"; }
  }
  my $acefile = "$database/acefiles/WBgene_spans_${chrom}.ace";
  unless ( $no_ace){ 
    open (ACE,">$acefile") or do{ $log->write_to("cant open output $acefile:\t$!\n"); die "cant open output $acefile:\t$!\n";}
  }
  foreach my $gene ( keys %gene_span ) {

    if( $gff ) {
      print OUTGFF "CHROMOSOME_${chrom}\tgene\tgene\t";
      print OUTGFF $gene_span{$gene}->{'min'},"\t",$gene_span{$gene}->{'max'},"\t";
      print OUTGFF ".\t",$gene_span{$gene}->{'strand'},"\t.\tGene \"$gene\"\n";
    }

    # write S-map details back to database
    unless ($no_ace) {
      my @coords;
      @coords = $gene_span{$gene}->{'strand'} eq "+" ? $coords->LocateSpan($chrom,$gene_span{$gene}->{'min'},$gene_span{$gene}->{'max'}) : $coords->LocateSpan($chrom,$gene_span{$gene}->{'max'},$gene_span{$gene}->{'min'});
      print ACE "\nSequence : \"$coords[0]\"\n";
      print ACE "Gene_child $gene $coords[1] $coords[2]\n";
    } 
  }
  close ACE;
  close OUTGFF;
  &load_to_database("$database","$acefile","WBGene_span") unless ( $no_ace ); # I know this would be better done all together but . . 
}




exit(0); 


=pod

=Title

  WBGene_span.pl

=head1 Overview

  parses GFF files to create a gene span for each WBGene from the start of the 5' most point of all transcripts to the 3'.

=head1 Output

  writes acefile output to $database/acefiles and gff to $database/CHROMOSOMES
  if acefile is created it will be loaded in to the database

=head1 Options

  -database  default is /wormsrv2/autoace
  -test      to use the test database
  -gff       to write a GFF file of the genespans as they are created rather than dumping them after loading to database
  -no_ace    dont write an acefile

=head1 Example

 perl WBGene_span.pl -database ~wormpub/DATABASES/current_DB -gff

=cut
