#!/usr/bin/perl

# Generate a GFF format of the orthologs
# This file is used for generating plots and such

#  Todd Harris (harris@cshl.org)
#  DATE : 05 May 2003


use Parse;
use strict;

use constant VER => 'orthologs-2.00';

$|++;
my $file = shift;
my $parse = Parse->new();

# Fetch and order all elegans genes and relevant information
my $elegans    = $parse->elegans('protein');
my $briggsae   = $parse->briggsae('gene');
my $orthologs  = $parse->orthologs_temp($file);
my $alignments = $parse->sw();

open OUT,">orthologs-" . VER;
print OUT "////////////////////////////////////////////\n";
print OUT "//     C. elegans/C. briggsae orthologs\n";
print OUT "// Generated : ",`date`,"\n";
print OUT "// Version   : ",VER,"\n";
print OUT '// Todd Harris (harris@cshl.org)',"\n";
print OUT "// elegans briggsae evalue SW_percent_identity confidence method\n";
print OUT "////////////////////////////////////////////\n";




foreach my $ortho (sort keys %$orthologs) {
  my $ce = $orthologs->{$ortho}->{ce};
  my $cb = $orthologs->{$ortho}->{cb};
  
  # Only want to print out ce->cb (orthos hash contains both directions)
  next if ($briggsae->{$ortho});
  
  print STDERR $ce,"\t",$cb,"\n";

  # Fetch out the percent identity for this pair
  my @als = eval { @{$alignments->{$cb}} };
  my $per_id;
  foreach my $al (@als) {
    if ($al->{ce} eq $ce) {
      $per_id = $al->{per_id};
      last;
    }
  }
  
  # Fetch out the start and stop positions of both genes.
  my $chrom      = $elegans->{$ce}->{chrom};
  my $ce_start   = $elegans->{$ce}->{chrom_start};
  my $ce_stop    = $elegans->{$ce}->{chrom_stop};
  my $ce_strand  = $elegans->{$ce}->{strand};
  $ce_strand = ($ce_strand eq '+1') ? '+' : '-';

  my $cb_start = $briggsae->{$cb}->{start};
  my $cb_stop = $briggsae->{$cb}->{stop};
  my $super   = $briggsae->{$cb}->{supercontig};
  my $cb_strand  = $briggsae->{$cb}->{strand};

  ($ce_start,$ce_stop) = ($ce_stop,$ce_start) if ($ce_stop < $ce_start);
  ($cb_start,$cb_stop) = ($cb_stop,$cb_start) if ($cb_stop < $cb_start);

  print join("\t",
	     $chrom,
	     $ce,
	     $ce_start,
	     $ce_stop,
	     $ce_strand,
	     $super,
	     $cb,
	     $cb_start,
	     $cb_stop,
	     $cb_strand,
	     $orthologs->{$ortho}->{eval},
	     $per_id,
	     $orthologs->{$ortho}->{conf},
	     $orthologs->{$ortho}->{meth}),"\n";
  
  my $meth = $orthologs->{$ortho}->{meth};
  my $conf   = $orthologs->{$ortho}->{conf};
  # Rest the confidence values if the method is by synteny

  unless ($meth eq 'seg-on' || $meth eq 'seg-off') {
    $conf = '-';
  }
  
  print OUT join("\t",
		 $ce,
		 $cb,
		 $orthologs->{$ortho}->{eval},
		 $per_id,
		 $conf,
		 $meth),"\n";
}

close OUT;
