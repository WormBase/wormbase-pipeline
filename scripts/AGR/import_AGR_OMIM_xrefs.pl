#!/bin/env perl
# trawl throught he rest endpoint and import from the AGR OMIM Xrefs into WormBase
# potential issues:
#  * the REST endpoint will change regularly (or its payload)
#  * AGM schema will also change regulary (https://github.com/alliance-genome/agr_schemas/tree/master/ingest/disease)
#
# URLs:
# * https://fms.alliancegenome.org/api/datafile/by/DAF?latest=true
# * https://download.alliancegenome.org


use LWP::Simple;
use JSON;
use Wormbase;
use Log_files;
use Storable;
use IO::File;
use Getopt::Long;
use Compress::Zlib;
use strict;

my ($debug,$test,$store,$wormbase,$load,$outfile);
GetOptions( 'debug=s'      => \$debug,
            'test'         => \$test,
            'store=s'      => \$store,
	    'load'         => \$load,
            'outfile=s'    => \$outfile,
    )||die("invalid commandline options\n");

if ( $store ) {
    $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
    $wormbase = Wormbase->new(
	-debug   => $debug,
	-test    => $test,
	);
}

my $log = Log_files->make_build_log($wormbase); # establish log file.

$outfile||=$wormbase->acefiles . '/AGR_OMIM_xrefs.ace';
my $outfh = IO::File->new($outfile,'w');
$outfh->binmode(':utf8'); # that might bite us when importing into acedb
$log->write_to("writing to $outfile\n\n");

my $content = get('https://fms.alliancegenome.org/api/datafile/by/DAF?latest=true');
die "Couldn't query the REST endpoint\n" unless $content;


my $latest = decode_json($content);
foreach my $db (@$latest){
    if ($db->{dataSubType}->{name} eq 'HUMAN'){
	my $downloadPath='https://download.alliancegenome.org/'.$db->{s3Path};
	$log->write_to("$1 => $downloadPath\n");
	get_agm($downloadPath);
    }
}

$wormbase->load_to_database($wormbase->autoace,$outfile,'AGR_import',$log) if $load;
$log->mail();


sub get_agm{
    my ($url)=@_;

    my $payload = get($url);
    $log->log_and_die("couldn't get human disease data\n") unless $payload;
    my $unzipped =Compress::Zlib::memGunzip($payload);

    my %genes_processed;
    my $json = decode_json($unzipped);
    foreach my $obj (@{$json->{data}}){
	my $gene_id = $obj->{objectId};
	next if exists $genes_processed{$gene_id};
	for my $prov(@{$obj->{dataProvider}}) {
	    for my $x ($prov->{crossReference}) {
		if ($x->{id} =~ /^OMIM:(\w+)$/) {
		    print $outfh "Gene: \"$gene_id\"\n";
		    print $outfh "Database \"OMIM\" \"OMIM_id\" \"$1\"\n\n";
		    $genes_processed{$gene_id} = 1;
		    last;
		}
	    }
	}
    }
}

