#!/usr/local/bin/perl -w
use strict;
use lib '../blib/lib';
use lib '/nfs/WWWdev/SANGER_docs/lib/Projects/C_elegans';
use NameDB_handler;
use Getopt::Long;

=pod

=head batch_addname.pl

=item Options:

  -user      username
  -password  password
  -file 	 file containing list of GeneIDs and CGC name eg WBGene00008040 ttr-5
  -species   what species these are for - default = elegans
  -test      use the test nameserver

e.g. perl batch_addname.pl -u fred -p secret -file genenames.txt -species briggsae

=cut

my ($USER,$PASS, $test, $file, $species);
GetOptions(
	   'user:s'     => \$USER,
	   'password:s' => \$PASS,
	   'test'       => \$test,
	   'file:s'     => \$file,
	   'species:s'  => \$species
	  ) or die;

$species = 'elegans' unless $species;

my $DB;
if ($test) {
    $DB = 'test_wbgene_id;mcs2a';
  } else {
    $DB = 'wbgene_id;mcs2a';
}

my $db = NameDB_handler->new($DB,$USER,$PASS);

$db->setDomain('Gene');

#open file and read

open (FILE,"<$file") or die "can't open $file : $!\n";
while(<FILE>) {
	my($id,$name) = split;
	$db->add_name($id,$name,'CGC',$species);
}

















