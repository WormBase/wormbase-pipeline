#!/usr/local/bin/perl5.6.0 -w                    # perl5.6.0 and -w flag

use strict;                                      # always use
use Wormbase;


open (HISTORY, "/wormsrv2/WORMPEP/wormpep81/wormpep.history81");
my %newpep_oldpep;
my %store;#hash of arrays GENE => protein, in, out

#array index variables to make script readable
my $protein = 0;
my $in = 1;
my $out = 2;

#read file in 
while(<HISTORY>)
      {
	my @data = split(/\s+/,$_);
	my $key = shift(@data);
	if( defined($store{$key}) )
	  {
	    if ($data[$in] >= $store{$key}[$out])
	      {
		#note in replacement hash
		my $new_protein = $data[$protein];
		my $old_protein = $store{$key}[$protein];
		$newpep_oldpep{$new_protein} = $old_protein;
	      }
	  }
  	$store{$key} = [ @data ];
      }

close HISTORY;

foreach my $ce (keys %newpep_oldpep)
  {
    print "$ce replaces $newpep_oldpep{$ce}\n"
  }

foreach my $geen (keys %store)
  {
    print "$geen $store{$geen}[0]\n";
  }

open (PEPACE,">/nfs/disk56/ar2/pepace.ace") || die "cant write pepace.ace";
open (NEWLIVE,">/nfs/disk56/ar2/newlive.ace") || die "cant write newlive.ace";
print NEWLIVE "KeySet : \"new_pepace\"\n";
foreach my $gene(keys %store)
  {
    my $peptide = $store{$gene}[$protein];
    print PEPACE "Protein : \"WP:$peptide\"\n";
    print NEWLIVE "Protein : \"WP:$peptide\"\n";
    
    if (defined($newpep_oldpep{$peptide})) {
      if ($peptide ne $newpep_oldpep{$peptide}) {
	print PEPACE "Replaces \"WP:$newpep_oldpep{$peptide}\"\n";
	undef $newpep_oldpep{$peptide};
      }
    }
    unless( defined($store{$gene}[$out]) )#was ever removed
      {
	print PEPACE "Corresponding_DNA \"$gene\"\n";
	print PEPACE "Live\n";
      }

    print PEPACE "\n";
  }

#now throw out those proteins whose replacements were replaced.
foreach my $CE (keys %newpep_oldpep)
  {
    if( defined ($newpep_oldpep{$CE}) ) {
      print PEPACE "Protein : \"WP:$newpep_oldpep{$CE}\"\n";
      print PEPACE "Replaced_by \"WP:$CE\"\n\n";
      print NEWLIVE "Protein : \"WP:$newpep_oldpep{$CE}\"\n";
      print NEWLIVE "Protein : \"WP:$CE\"\n";
    }
  }
exit(0);
