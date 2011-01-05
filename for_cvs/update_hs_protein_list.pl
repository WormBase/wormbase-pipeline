#!/usr/bin/perl
# from Norie de la Cruz
# Version 2nd February 2009

# usage:
#   generates/updates a list of the all the proteins as well as a list of human proteins.
#   You need to start it off with an empty file.

use strict;
use Ace;

my $list_dir = "./"; ## orthology_data/integration_test
my $all_protein_list = "$list_dir/all_proteins.txt";
# my $protein_test_out = "$list_dir/all_proteins_test.txt";
my $hs_protein_list = "$list_dir/hs_proteins.txt";

open ALL_PROTEIN_LIST, "$all_protein_list" or die "Cannot open all protein list\n";


my $DB = Ace->connect(-host=>'aceserver.cshl.org',-port=>2005);


## build protein hash

my %all_proteins;

foreach my $protein_id (<ALL_PROTEIN_LIST>) {
	
	chomp $protein_id;
	$all_proteins{$protein_id} = 1;
#	print "$protein_id\n";

}

close ALL_PROTEIN_LIST;

open ALL_PROTEIN_TEST_OUT, ">> $all_protein_list" or die "Cannot open test protein list\n";
open HS_PROTEIN_LIST, ">> $hs_protein_list" or die "Cannot open hs protein list\n";

### get and check protein data 


my @acedb_proteins = $DB->fetch(-class=>'Protein');

foreach my $ace_protein (@acedb_proteins) {

	if ($all_proteins{$ace_protein}) {
	
		next;
	
	} else {
	
			#print "$ace_protein\n";
		 	my $sp = $ace_protein->Species;
                
        	if($sp =~ m/sapien/){
             
             	print ALL_PROTEIN_TEST_OUT "$ace_protein\n";
             	print HS_PROTEIN_LIST "$ace_protein\n";
             
             } else {
               
               	print ALL_PROTEIN_TEST_OUT "$ace_protein\n";
                
           	}


	}
}

print "done\n";
