#!/bin/env perl
# import C.elegans AGR orthology data
#  from REST endpoint: https://fms.alliancegenome.org/api/datafile/by/ORTHO?latest=true


use lib $ENV{CVS_DIR};

use LWP::Simple;
use JSON;
use IO::File;
use Log_files;
use Storable;
use Wormbase;
use IO::File;
use Getopt::Long;
use strict;
use Compress::Zlib;

my ($debug,$test,$store,$wormbase,$load,$outfile);
GetOptions( 'debug=s'    => \$debug,
            'test'       => \$test,
            'store=s'    => \$store,
            'load'       => \$load,
            'outfile=s'  => \$outfile,
            )||die("invalid commandline options\n");

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             );
}

my $log = Log_files->make_build_log($wormbase); # establish log file.

$outfile||=$wormbase->acefiles . '/AGR_orthologs.ace';
my $outfh = IO::File->new($outfile,'w');
$log->write_to("writing to $outfile\n\n");

my %taxon2name = (
        9606   => 'Homo sapiens',            # human
        7227   => 'Drosophila melanogaster', # fly
        10090  => 'Mus musculus',            # mouse
        10116  => 'Rattus norvegicus',       # rat
        7955   => 'Danio rerio',             # fish
        559292 => 'Saccharomyces cerevisiae',# mushroom
        4932   => 'Saccharomyces cerevisiae',# also mushroom
	6239   => 'Caenorhabditis elegans',  # worm
);

my %taxon2prefix = (
        9606   => '',      # human
        7227   => 'FB:',   # fly
        10090  => '',      # mouse
        10116  => '',      # rat
        7955   => 'ZFIN:', # fish
        559292 => 'SGD:',  # mushroom
	4932   => 'SGD:',  # also mushroom
	6239   => '',      # worm
);

my %gene_ids = reverse $wormbase->FetchData('worm_gene2geneID_name');

my $content = get('https://fms.alliancegenome.org/api/datafile/by/ORTHO?latest=true');
$log->log_and_die("Couldn't query the REST endpoint\n") unless $content;

my $latest = decode_json($content);

foreach my $db (@$latest){
	if ($db->{dataSubType}->{name} eq 'WB'){
		get_orthos('https://download.alliancegenome.org/'.$db->{s3Path});
	}
}


$wormbase->load_to_database($wormbase->autoace,$outfile,'AGR_orthologs',$log) if $load;
$log->mail;

sub get_orthos{
	my ($url)=@_;
	my $payload = get($url);
	$log->log_and_die("cannot download orthologs from $url\n") unless $payload;
	$payload = Compress::Zlib::memGunzip($payload) if $url =~ /\.gz$/;
	my $json = decode_json($payload);
        foreach my $data(@{$json->{data}}){
		next unless $data->{strictFilter};
		$data->{gene1}=~s/DRSC://;
		$data->{gene2}=~s/DRSC://;
		if ($data->{gene1}||$data->{gene2}){ # skip pairs that don't exist in WB
		    next unless exists($taxon2name{$data->{gene1Species}}) && exists($taxon2name{$data->{gene2Species}}); # skip unrecognised species
		    print_ace($data->{gene1},$data->{gene1Species},$data->{gene2},$data->{gene2Species},$data->{predictionMethodsMatched});
		    print_ace($data->{gene2},$data->{gene2Species},$data->{gene1},$data->{gene1Species},$data->{predictionMethodsMatched});
		} else {
		    $log->write_to("cannot find ${\$data->{gene1}} / ${\$data->{gene2}}\n");
		}
	}
}

sub print_ace{
	my ($g1,$s1,$g2,$s2,$evidence)=@_;
	print $outfh "Gene : $taxon2prefix{$s1}$g1\n";
        map {
	 s/\s/\-/g;
	 print $outfh "Ortholog $taxon2prefix{$s2}$g2 \"$taxon2name{$s2}\" From_Analysis $_\n";
	} @$evidence;
	print $outfh "\n";
}
