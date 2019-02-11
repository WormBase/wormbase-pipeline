#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;

use lib $ENV{CVS_DIR};
use Wormbase;
use Modules::AGR;

my ($help, $debug, $test, $verbose, $store, $wormbase, $bgi_json);
my ($gff_in, $gff_out, $ws_version);

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "test"        => \$test,
	    "verbose"     => \$verbose,
	    "store:s"     => \$store,
	    "gffout:s"    => \$gff_out,
            "gffin=s"     => \$gff_in,
            "bgijson=s"   => \$bgi_json,
            "wsversion=s" => \$ws_version,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
      );
}

die "You must supply an input GFF file\n" if not defined $gff_in or not -e $gff_in;
die "You must supply an output file name\n" if not defined $gff_out;
die "You must supply an BGI JSON file\n" if not defined $bgi_json;

my %bgi_genes = %{&get_bgi_genes($bgi_json)};

my $in_fh = &open_gff_file($gff_in);
open(my $out_fh, ">$gff_out") or die "Could not open $gff_out for writing\n";
print $out_fh "# WormBase release $ws_version\n" if defined $ws_version;

while(<$in_fh>) {
  /^\#/ and do {
    print $out_fh $_;
    next;
  };

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
    $l[8] = join(";", map { $_ . "=" . $attr{$_} } keys %attr);
    print $out_fh join("\t", @l), "\n";

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

