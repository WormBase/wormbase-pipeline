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

my ($database,$alleles,$outfile,$ws_version,$diseases);
GetOptions (
	'database:s' => \$database,
	'alleles:s'  => \$alleles,
	'diseases:s' => \$diseases,
	'outfile:s'  => \$outfile,
	'wsversion:s'=> \$ws_version,
)||die @!;

my $db = Ace->connect(-path => $database) or die(Ace->error);

my @ids = read_ids($alleles);
my @disease_ids = read_disease_ids($diseases);

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
	warn("cannot find $id\n") unless $aceVar;
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
}

foreach my $dId (@disease_ids){
    if ($dId =~ /^WBGenotype/) {
	my $d = $db->fetch(Genotype => $dId);
	my $disease_annotation = {
		primaryID => "WB:$dId",
		name      => ($d->Genotype_name ? "${\$d->Genotype_name}" : $dId),
		taxonId   => 'NCBITaxon:6239',
		crossReference => {id => "WB:$dId",pages => ['genotype']} # needs to be added to the resources file
	};
	push @annotation, $disease_annotation;
    }
    else {
	if (exists $strains{$dId}) {
            warn "$dId already retrieved from alleles file, ignoring entry in disease file\n";
            next;
        }
	my $s = $db->fetch(Strain => $dId);
	my $strain = {
	    primaryID => "WB:$dId",
            name      => ($s->Public_name ?"${\$s->Public_name}" :$dId),
            taxonId   => 'NCBITaxon:6239',
            crossReference => {id => "WB:$dId", pages => ['strain']},
        };
        push @annotation, $strain;
    }
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

sub read_disease_ids{
	my ($file) = @_;

	my %seen;
	open IN, "<$file";
	my @ids;
	while (<IN>){
	    if ((/(WBGenotype\d+)/ or /(WBStrain\d+)/) and !exists $seen{$1}) {
		push @ids, "$1" ;
		$seen{$1}++;
	    }
	}
	return @ids;
}

sub read_ids{
	my ($file)=@_;
	my $jsondoc = read_json($file);
	my @ids = map {$_->{primaryId}} @{$jsondoc->{data}};
	return @ids;
}

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
