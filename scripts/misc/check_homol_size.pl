#!/usr/local/bin/perl -w

use strict;
use Ace;
use lib $ENV{'CVS_DIR'};

my $database = shift;
my $acefile = shift;
open (ACE,">$acefile") or die "cant open $acefile $!\n" if $acefile;

my $db = Ace->connect( -path => $database) or die Ace->error();

my @badones;
my $clones = $db->fetch_many(-query=>'find Sequence; Genomic_canonical; species = "Caenorhabditis elegans"');
while ( my $clone = $clones->next ) {
  my $clonesize = length($clone->DNA->right->name);

  # check Homol_data
  my @homols = $clone->Homol_data;
  &check( \@homols, $clonesize, $clone->name  );

  # . . and Feature_data
  @homols = $clone->Feature_data;
  &check(\@homols, $clonesize, $clone->name );
}
$db->close;

print join("\n",@badones);
close ACE if $acefile;

sub check
  {
    my $homols= shift;
    my $clonesize = shift;
    my $clone = shift;
    foreach ( @$homols ) {
      if( $_->right and $_->right->right ){ 
	my $length = $_->right->right->name;
	if ($length != $clonesize ) {
	  print "\nERROR : ",$_->name," discrepant length\n\n";
	  print ACE "//discrepant length\nSequence : $clone\nHomol_data ",$_->name," 1\t$clonesize\n" if $acefile;
	  push(@badones,$_->name);
	}
      }
      else {
	print "\nERROR : ",$_->name," missing coordinates\n\n";
	  print ACE "//missing coordinates\nSequence : $clone\nHomol_data ",$_->name," 1\t$clonesize\n" if $acefile;
	push(@badones,$_->name);
      }
    }
  }
