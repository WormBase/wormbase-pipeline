#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;

use lib $ENV{CVS_DIR};
use Wormbase;
use Modules::AGR;

my ($debug, $test, $verbose, $store, $wormbase);
my ($outfile, $acedbpath, $ws_version, $out_fh, $bgi_json,$disease_file);

GetOptions (
  'debug=s'       => \$debug,
  'test'          => \$test,
  'verbose'       => \$verbose,
  'store:s'       => \$store,
  'database:s'    => \$acedbpath,
  'outfile:s'     => \$outfile,
  'wsversion=s'   => \$ws_version,
  'bgijson=s'     => \$bgi_json,
  'diseasefile=s' => \$disease_file,
)||die();

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
      );
}

my $tace = $wormbase->tace;
my $date = AGR::get_rfc_date();
my $alt_date = join("/", $date =~ /^(\d{4})(\d{2})(\d{2})/);
my $taxid = $wormbase->ncbi_tax_id;
my $full_name = $wormbase->full_name;

$acedbpath  ||= $wormbase->autoace;
$ws_version ||= $wormbase->get_wormbase_version_name;
$outfile    ||= "./wormbase.agr_allele.${ws_version}.json";

my (@alleles, @annots);
my %objectsIds;
my $bgi_genes = AGR::get_bgi_genes( $bgi_json ) if $bgi_json;

my $db = Ace->connect(-path => $acedbpath, -program => $tace) or die("Connection failure: ". Ace->error);

my $it = $db->fetch_many(-query => 'find Variation WHERE Live AND COUNT(Gene) == 1 AND (Phenotype OR Disease_info OR Interactor) AND NOT Natural_variant');
process($it);

$it = $db->fetch_many(-query => 'find Transgene');
process_transgenes($it);

process_disease_variations(grab_disease_variations($disease_file)) if $disease_file;

sub process_disease_variations{
	my @vars = @_;
	foreach my $v(@vars){
	        $it = $db->fetch_many(Variation => $v);
		process($it);
	}
}

sub grab_disease_variations{
	my ($file) = @_;
	my @variations;

	my $in = IO::File->new($file,'r');
	while(<$in>){
		push(@variations,"$1") if /(WBVariation\d+)/
	}
	return @variations;
}

sub process{
  my ($it) = @_;

  while (my $obj = $it->next) {
    next unless $obj->isObject();
    
    if ($objectsIds{"$obj"}){
	    next
    }else{
	    $objectsIds{"$obj"}=1;
    }

    my $has_interaction;
    my $has_disease = 1 if $obj->Disease_info;
    my $has_phenotype = 1 if $obj->Description;
    if ($obj->Interactor) {
      foreach my $item ($obj->Interactor->col) {
        $has_interaction = 1 if $item eq 'Regulatory';
      }
    }
  
    next unless $has_disease or $has_phenotype or $has_interaction or $obj->Corresponding_transgene;

    my ($gene) = $obj->Gene->name if $obj->Gene;
    $gene ||= $obj->Corresponding_transgene->Construct->Gene if $obj->Corresponding_transgene && $obj->Corresponding_transgene->Construct;

    next if defined $bgi_genes and not exists $bgi_genes->{"WB:$gene"};

    my $symbol = $obj->Public_name->name;
    my %synonyms = map {$_->name => 1}$obj->Other_name;

    my $json_obj = {
      primaryId     => "WB:$obj", 
      symbol        => $symbol,
      symbolText    => $symbol,
      synonyms      => [keys \%synonyms],
      secondaryIds => [],
      taxonId       => "NCBITaxon:" . $taxid,
      crossReferences => [ { id => "WB:$obj", pages => ['allele','allele/references']}],
    };
    map { push @{$json_obj->{crossReferences}}, {id => "WB:$_", pages => ['reference']} } $obj->Reference;
    $$json_obj{alleleObjectRelations}=[{objectRelation => {associationType => 'allele_of', gene => "WB:$gene"}}];

    if ($obj->Corresponding_transgene){
	    my $transgene = $obj->Corresponding_transgene;
	    my $construct = $transgene->Construct;
	    next unless $construct;

	    push @$json_obj{alleleObjectRelations} , {objectRelation => {associationType => 'contains',construct => "WB:$construct"}};
    }

    push @alleles, $json_obj;
  }
}

sub process_transgenes{
  my ($it) = @_;

  while (my $obj = $it->next) {
    next unless $obj->isObject();

    my $interaction;
    my $disease = $obj->Disease_info;
    my $phenotype = $obj->Phenotype;
    if ($obj->Interactor) {
      foreach my $item ($obj->Interactor->col) {
        $interaction = 1 if $item eq 'Regulatory';
      }
    }
  
    next unless $disease or $phenotype or $interaction;

    my $symbol = $obj->Public_name ? $obj->Public_name->name : "$obj";
    my %synonyms = map {$_->name => 1} $obj->Synonym;

    my $json_obj = {
      primaryId     => "WB:$obj", 
      symbol        => $symbol,
      symbolText    => $symbol,
      synonyms      => [keys \%synonyms],
      secondaryIds => [],
      taxonId       => "NCBITaxon:" . $taxid,
      crossReferences => [ { id => "WB:$obj", pages => ['transgene','transgene/references'] }],
    };
    $json_obj->{description} = "${\$obj->Summary}" if $obj->Summary;
    map { push @{$json_obj->{crossReferences}}, {id => "WB:$_", pages => ['reference']} } $obj->Reference;

    # $json_obj{alleleObjectRelations}=[{objectRelation => {associationType => 'contains',construct => "WB:$construct"}}];
    $json_obj->{alleleObjectRelations}=[] if $obj->Construct;
    map {push @{$json_obj->{alleleObjectRelations}},{objectRelation => {associationType =>'contains',construct => "WB:$_"}}} $obj->Construct;

    push @alleles, $json_obj;
  }
}

my $data = {
  metaData => AGR::get_file_metadata_json( (defined $ws_version) ? $ws_version : $wormbase->get_wormbase_version_name() ),
  data     => \@alleles,
};

if (defined $outfile) {
  open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";
} else {
  $out_fh = \*STDOUT;
}

my $json_obj = JSON->new;
my $string = $json_obj->allow_nonref->canonical->pretty->encode($data);
print $out_fh $string;

$db->close;
