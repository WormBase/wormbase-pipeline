#!/usr/local/bin/perl5.6.0 -w                    # perl5.6.0 and -w flag

use strict;                                      # always use
use Wormbase;


open (HISTORY, "/wormsrv2/WORMPEP/wormpep82/wormpep.history82");
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
		#note in replacement hash (if not being replaced by itself (reappeared)
		my $new_protein = $data[$protein];
		my $old_protein = $store{$key}[$protein];
		if( "$new_protein" ne "$old_protein" ){
		  $newpep_oldpep{$new_protein} = $old_protein;
		}
	      }
	  }
  	$store{$key} = [ @data ];
      }$ARGV[0]

close HISTORY;

foreach my $ce (keys %newpep_oldpep)
  {
    print "$ce replaces $newpep_oldpep{$ce}\n"
  }

foreach my $geen (keys %store)
  {
    print "$geen $store{$geen}[0]\n";
  }
my $kill_pep;

open (PEPACE,">/nfs/disk56/ar2/pepace.ace") || die "cant write pepace.ace";
foreach my $gene(keys %store)
  {
    my $peptide = $store{$gene}[$protein];
    print PEPACE "Protein : \"WP:$peptide\"\n";
    
    if (defined($newpep_oldpep{$peptide})) {
      if ($peptide ne $newpep_oldpep{$peptide}) {
	print PEPACE "Replaces \"WP:$newpep_oldpep{$peptide}\"\n";
	$kill_pep = "WP:$newpep_oldpep{$peptide}";
	undef $newpep_oldpep{$peptide};
      }
    }
    unless( defined($store{$gene}[$out]) )#was ever removed
      {
	print PEPACE "Corresponding_DNA \"$gene\"\n";
	print PEPACE "Live\n";
      }
    else {
      print PEPACE "-D Live\n-D Corresponding_DNA\n";
    }
	
    print PEPACE "\n";

    if( defined($kill_pep) ){
      print PEPACE "Protein : \"$kill_pep\"\n-D Live\n-D Corresponding_DNA\n";
      print PEPACE "\n";
      undef $kill_pep;
    }
  }

#now throw out those proteins whose replacements were replaced.
foreach my $CE (keys %newpep_oldpep)
  {
    if( defined ($newpep_oldpep{$CE}) ) {
      print PEPACE "Protein : \"WP:$newpep_oldpep{$CE}\"\n";
      print PEPACE "Replaced_by \"WP:$CE\"\n"; 
      print PEPACE "-D Live\n";
      print PEPACE "-D Corresponding_DNA\n";
      print PEPACE "\n";
    }
  }
exit(0);
