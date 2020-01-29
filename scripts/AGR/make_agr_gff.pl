#!/usr/bin/env perl

use strict;
use Getopt::Long;

use lib $ENV{CVS_DIR};
use Modules::AGR;

my ($bgi_json, $gff_in, $gff_out, $ws_version);

GetOptions (
	    "gffout:s"    => \$gff_out,
            "gffin=s"     => \$gff_in,
            "bgijson=s"   => \$bgi_json,
            "wsversion=s" => \$ws_version,
	   )||die();

$gff_in   || die "You must supply an input GFF file\n";
$gff_out  || die "You must supply an output file name\n";
$bgi_json || die "You must supply an BGI JSON file\n";

my %bgi_genes = %{AGR::get_bgi_genes($bgi_json)};
my $date      = AGR::get_rfc_date();

my $in_fh = &open_gff_file($gff_in);
open(my $out_fh, ">$gff_out") or die "Could not open $gff_out for writing\n";
print $out_fh "# WormBase release $ws_version\n" if defined $ws_version;

my %headers; # do not print redundant headers

while(<$in_fh>) {
  if (/^\#/){
    print $out_fh $_ unless $headers{$_};
    $headers{$_}=1;
    if (/gff-version 3/){
	    print $outfh<<HERE;
#!date-prodiced $date
#!data-source WormBase
#!assembly WBcel235
HERE
    }
    next;
  }

  if (/^\S+\s+WormBase\s+gene\s+/) {
    chomp;
    my @l = split(/\t/, $_);
    my %attr;
    foreach my $kv (split(/;/, $l[8])) {
      my ($k, $v) = split(/\=/, $kv);
      $attr{$k} = $v;
    }
    my $gid = $attr{curie};

    # use symbol as Name
    $attr{Name} = $bgi_genes{$gid}->{symbol};
    if (exists $bgi_genes{$gid}->{name}) {
      $attr{long_name} = $bgi_genes{$gid}->{name};
    }
    if (exists $bgi_genes{$gid}->{geneSynopsis}) {
      $attr{description} = $bgi_genes{$gid}->{geneSynopsis};
    }
    $attr{Ontology_term}='SO:0000704';

    $l[8] = join(";", map { $_ . "=" . $attr{$_} } keys %attr);
    print $out_fh join("\t", @l), "\n";

  } elsif(/^\S+\sWormBase\smRNA/){
          my @l = split(/\t/, $_);
	  my $ontologyID = 'SO:0000234';
	  /Name=([^;]+)/;
	  my $transcriptID="$1";
	  /(WBGene\d+)/;
          my $geneID="$1";
	  s/;Parent/curie=WB:$transcriptID;curie=WB:$geneID;Ontology_term=$ontologyID;Parent/;
	  print $out_fh $_;

  } elsif (/^\S+\s+WormBase\s+/) {
    print $out_fh $_;
  }
}
close($out_fh) or die "Could not close $gff_out after writing (probably file system full)\n";

exit(0);

sub open_gff_file {
  my ($file) = @_;

  my $fh;
  if ($file =~ /\.gz$/) {
    open( $fh, "gunzip -c $file |") or die "Could not open gunzip stream to $file\n";
  } else {
    open( $fh, $file) or die "Could not open $file for reading\n";
  }

  return $fh;
}

