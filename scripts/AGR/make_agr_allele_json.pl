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
my ($outfile, $acedbpath, $ws_version, $out_fh, $bgi_json);

GetOptions (
  "debug=s"     => \$debug,
  "test"        => \$test,
  "verbose"     => \$verbose,
  "store:s"     => \$store,
  "database:s"  => \$acedbpath,
  "outfile:s"   => \$outfile,
  "wsversion=s" => \$ws_version,
  "bgijson=s"   => \$bgi_json,
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

my $bgi_genes = AGR::get_bgi_genes( $bgi_json ) if $bgi_json;

my $db = Ace->connect(-path => $acedbpath, -program => $tace) or die("Connection failure: ". Ace->error);

my $it = $db->fetch_many(-query => 'find Variation WHERE Live AND COUNT(Gene) == 1 AND (Phenotype OR Disease_info OR Interactor) AND NOT Natural_variant');
process($it);

# $it = $db->fetch_many(-query => 'find Variation WHERE Live AND Corresponding_transgene');
# process($it);

$it = $db->fetch_many(-query => 'find Transgene');
process_transgenes($it);

sub process{
  my ($it) = @_;

  while (my $obj = $it->next) {
    next unless $obj->isObject();

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
      gene          => "WB:$gene",
      crossReferences => [ { id => "WB:$obj", pages => ["allele"] }],
    };
    map { push @$json_obj->{crossReferences}, {id => "WB:$_" => ['reference']} } $obj->Reference;

    if ($obj->Corresponding_transgene){
	    my $transgene = $obj->Corresponding_transgene;
	    my $construct = $transgene->Construct;
	    next unless $construct;
	    $$json_obj{construct}="WB:$construct";
	    next unless $transgene->Genetic_information; # skip the transgenes where the type hasn't been curated
            if(grep {$_ eq 'Integrated'} $transgene->Genetic_information){
	        $$json_obj{constructInsertionType}='Transgenic Insertion'
	    }elsif(grep {$_ eq 'Extrachromosomal'} $transgene->Genetic_information){
                $$json_obj{constructInsertionType}='Extrachromosomal Array'
	    }else{next}
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

#   my ($gene) = $obj->Construct->Gene->name if $obj->Construct && $obj->Construct->Gene;
#   next if defined $bgi_genes and not exists $bgi_genes->{"WB:$gene"};

    my $symbol = $obj->Public_name ? $obj->Public_name->name : "$obj";
    my %synonyms = map {$_->name => 1} $obj->Synonym;

    my $json_obj = {
      primaryId     => "WB:$obj", 
      symbol        => $symbol,
      symbolText    => $symbol,
      synonyms      => [keys \%synonyms],
      secondaryIds => [],
      taxonId       => "NCBITaxon:" . $taxid,
#     gene          => "WB:$gene",
      crossReferences => [ { id => "WB:$obj", pages => ['transgene'] }],
    };
    $json_obj->{description} = "$\{$obj->summary}" if $obj->summary;
    map { push @$json_obj->{crossReferences}, {id => "WB:$_" => ['reference']} }$obj->Reference;

    my $construct = $obj->Construct;
    next unless $construct;

    $$json_obj{construct}="WB:$construct";
    
    next unless $obj->Genetic_information; # skip the transgenes where the type hasn't been curated
    if(grep {$_ eq 'Integrated'} $obj->Genetic_information){
        $$json_obj{constructInsertionType}='Transgenic Insertion'
    }elsif(grep {$_ eq 'Extrachromosomal'} $obj->Genetic_information){
        $$json_obj{constructInsertionType}='Extrachromosomal Array'
    }else{next} # skip the ones without a AGR type

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
