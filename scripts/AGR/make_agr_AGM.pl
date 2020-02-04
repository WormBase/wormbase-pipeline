#!/usr/bin/env perl
# create the AGM , which is basically strains
# https://github.com/alliance-genome/agr_schemas/tree/master/ingest/affectedGenomicModel

use lib $ENV{CVS_DIR};
use strict;
use Getopt::Long;
use JSON;
use Ace;
use lib $ENV{CVS_DIR};
use Modules::AGR;

my ($database,$alleles,$outfile,$ws_version);
GetOptions (
	'database:s' => \$database,
	'alleles:s'  => \$alleles,
	'outfile:s'  => \$outfile,
	'wsversion:s'=> \$ws_version,
)||die @!;

my $db = Ace->connect(-path => $database) or die(Ace->error);

my $jsondoc = read_json($alleles);
my @ids = map {$_->{primaryId}} @{$jsondoc->{data}};
my %strains;

while (my $id = shift @ids){
	$id =~s/WB://;
	my $aceVar;
	my @strains;
	if ($id=~/WBVar/){
      	    $aceVar = $db->fetch(Variation => $id);
	}elsif($id=~/WBTransgene/){
      	    $aceVar = $db->fetch(Transgene => $id);
	}
        @strains = $aceVar->Strain;

	next unless @strains; # for the alleles without curated strain data
	foreach my $strain(@strains){
		# skip the ones without disease annotation
                # next unless $strain->Disease_info;
	        $strains{"$strain"}||=[];
		push @{$strains{"$strain"}},"$id";
	}
}

my @annotation;
while (my($k,$v)= each %strains){
	my $s = $db->fetch(Strain => $k);
	my $strain={
		primaryID => "WB:$k",
		name      => ($s->Public_name ?"${\$s->Public_name}" :$k),
		taxonId   => 'NCBITaxon:6239',
		crossReference => {id => "WB:$k", pages => ['strain']},
		affectedGenomicModelComponents => [map {{alleleID => "WB:$_",zygosity => 'GENO:0000137'}} @$v],
	};
	push @annotation,$strain;
	# last;
}

my $date = AGR::get_rfc_date();
my $data = {
  metaData => AGR::get_file_metadata_json($ws_version, $date),
  data     => \@annotation,
};

my $json_obj = JSON->new;
my $string   = $json_obj->allow_nonref->canonical->pretty->encode($data);
my $out_fh;
if ($outfile) {
  open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";
} else {
  $out_fh = \*STDOUT;
}

print $out_fh $string;
$db->close;

###########################
sub read_json{
	my ($file) = @_;

	my ($fh,$buf);
	open($fh,$file);
	while (<$fh>){
		$buf.= $_
	}

	my $json_reader = JSON->new;
	return $json_reader->decode($buf);
}
