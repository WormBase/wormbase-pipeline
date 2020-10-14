#!/bin/env perl
# create a JSON version of the WormBase Molecule class
# spec is here: ....

# Some Thoughts
# * maybe skip CHEBI ?
# * maybe change the prop names to be in line with OBO ?

use Getopt::Long;
use Ace;
use JSON;
use lib $ENV{CVS_DIR};
use Modules::AGR;
use strict;

my ($database,$alleles,$outfile,$ws_version,$diseases,$skipchebi);
GetOptions (
        'database:s' => \$database,
        'outfile:s'  => \$outfile,
        'wsversion:s'=> \$ws_version,
	'skipchebi'  => \$skipchebi, # in case we have to skip entries with CHEBI ids
)||die @!;

my $db = Ace->connect(-path => $database) or die(Ace->error);

my $it = $db->fetch_many(Molecule => '*');
my @molecules;
while (my $mol = $it->next){

	my $chebi = get_chebi($mol);
	next if ($chebi && $skipchebi);
	
	my %m = (
	 id => 'WB:'.$mol->name,
	 name => "${\$mol->Public_name}",
         crossReferences => [{id => 'WB:'.$mol->name,pages => ['molecule']}] # do we need to add CHEBI xrefs?, Do we need an array?
        );

	if ($chebi){
		$m{id} = $chebi;
		push @{$m{crossReferences}},{id => $chebi,pages =>['entry']};
	}

	### that block is potentially not needed ###
	$m{iupac}     = "${\$mol->IUPAC}"    if $mol->IUPAC;
	$m{smiles}    = "${\$mol->SMILES}"   if $mol->SMILES;
	$m{inchi}     = "${\$mol->InChi}"    if $mol->InChi;
	$m{inchikey}  = "${\$mol->InChiKey}" if $mol->InChiKey;
	$m{formula}   = "${\$mol->Formula}"  if $mol->Formula;
	############################################

	if ($mol->Synonym) {
		my @syns = map{"$_"} $mol->Synonym;
	        $m{synonyms}  = \@syns;
        }

	push @molecules,\%m;
}

$db->close;

my $date = AGR::get_rfc_date();
my $data = {
  metaData => AGR::get_file_metadata_json($ws_version, $date),
  data     => \@molecules,
};

my $out_fh;
if ($outfile) {
  open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";
} else {
  $out_fh = \*STDOUT;
}

my $json_obj = JSON->new;
my $string   = $json_obj->allow_nonref->canonical->pretty->encode($data);
print $out_fh $string;

sub get_chebi{
	my ($obj)=@_;
	for my $d ($obj->Database) {
        	next unless $d->name eq 'ChEBI';
        	return 'CHEBI:' . $d->right->right->name;
    	}
}
