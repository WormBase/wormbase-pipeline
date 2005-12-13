#!/usr/local/bin/perl -w

use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;

my $simple;
GetOptions ( 'simple' => \$simple );

my %CE2gene;
my %brigpep2gene;

FetchData('wormpep2cds',\%CE2gene);
FetchData('brigpep2cds',\%brigpep2gene);

while( <> ) {
#  if( /(CBP\d+).*(CE\d{5})/ ) {
#    my $brigpep = $1;
#    my $wormpep = $2;
#    my $bgenes = $brigpep2gene{ $brigpep };
#    my $cegenes = $CE2gene{ $wormpep };
#    if( $bgenes and $cegenes ) {
#      foreach $cbgene (split(/\s/,$bgenes) ){
#	my $line = $_;
#	$line =~ s/$brigpep/$cbgene/;
#	foreach $cegene (split(/\s/,$cegenes) ){
#	  $line =~ s/$wormpep/$cegene/;
#	  print $line;
#	}
#      }
#    }
#  }

  if( $simple ) {
    if (/(CBP\d{5})/) {
      my $gene = $brigpep2gene{"$1"};
      s/(CBP\d{5})/$gene/ if $gene;
      print ;
    }
    else {print}
  }
  else {
    
    if ( /^\d+\s+(CBP\d+)/ ) {
      my $brigpep = $1;
      my $bgenes = $brigpep2gene{ $brigpep };
      if ( $bgenes ) {
	foreach $cbgene (split(/\s/,$bgenes) ) {
	  my $line = $_;
	  $line =~ s/$brigpep/$cbgene/;
	  print $line;
	}
      }
    } elsif ( /(CE\d{5}).*(CBP\d{5})/ ) {
      my $brigpep = $2;
      my $wormpep = $1;
      my $bgenes = $brigpep2gene{ $brigpep };
      my $cegenes = $CE2gene{ $wormpep };
      if ( $bgenes and $cegenes ) {
	foreach $cbgene (split(/\s/,$bgenes) ) {
	  my $line = $_;
	  $line =~ s/$brigpep/$cbgene/;
	  foreach $cegene (split(/\s/,$cegenes) ) {
	    $line =~ s/$wormpep/$cegene/;
	    print $line;
	  }
	}
      }
    }
  }
}
