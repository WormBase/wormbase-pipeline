#!/bin/env perl
# import OMIM Xrefs and DOids from AGR DAF file into WormBase
# potential issues:
#  * the REST endpoint will change regularly (or its payload)
#  * AGR schema will also change regulary (https://github.com/alliance-genome/agr_schemas/tree/master/ingest/disease)
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

my ($debug,$test,$store,$wormbase,$load,$outfile, $database);
GetOptions( 'debug=s'      => \$debug,
            'test'         => \$test,
	    'database=s'   => \$database,
            'store=s'      => \$store,
	    'load'         => \$load,
            'outfile=s'    => \$outfile,
    )||die("invalid commandline options\n");

if ( $store ) {
    $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
    $wormbase = Wormbase->new(
	-autoace  => $database,
	-debug    => $debug,
	-test     => $test,
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
    my $unzipped = Compress::Zlib::memGunzip($payload);

    my %xrefs;
    my $json = decode_json($unzipped);
    foreach my $obj (@{$json->{data}}) {
	my $gene_id = $obj->{objectId};
	$xrefs{$gene_id}{'DO'}{$obj->{DOid}} = 1 if $obj->{DOid};
	for my $prov(@{$obj->{dataProvider}}) {
	    for my $x ($prov->{crossReference}) {
		if ($x->{id} =~ /^O?MIM:(\w+)$/) {
		    $xrefs{$gene_id}{'OMIM'} = $1;
		    last;
		}
	    }
	}
    }

    for my $gene_id (keys %xrefs) {
	print $outfh "Gene : \"$gene_id\"\n";
	printf $outfh "Database \"OMIM\" \"gene\" \"%s\"\n", $xrefs{$gene_id}{'OMIM'} if exists $xrefs{$gene_id}{'OMIM'};
	for my $do_id (keys %{$xrefs{$gene_id}{'DO'}}) {
	    print $outfh "Database \"DO\" \"id\" \"$do_id\"\n";
	}
	print $outfh "\n";
    }

    return;
}

