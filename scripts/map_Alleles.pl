#!/usr/bin/env perl
#
#  Created by Michael Han on 2006-08-10.
#  Copyright (c) 2006. All rights reserved.

use strict;
use lib $ENV{'CVS_DIR'};
use lib "$ENV{CVS_DIR}/Modules";
use map_Alleles;
use Wormbase;
use Getopt::Long;
use IO::File;    

my ( $debug, $store, $outdir,$allele ,$noload,$database);

GetOptions(
    'debug=s'  => \$debug,
    'store=s'  => \$store,
	'outdir=s' => \$outdir,
	'allele=s' => \$allele,
	'noload'   => \$noload,
	'noupdate' => \$noload,
	'database' => \$database,
);

# WormBase template
my $maintainer = 'All';
my $wb;
$debug=$ENV{LOGNAME} if $debug==1;

if ($store) {
    $wb = Storable::retrieve($store)
      or croak("cannot restore wormbase from $store");
}
else { $wb = Wormbase->new( -debug => $debug, -test => $debug, -autoace => $database ) }

my $log = Log_files->make_build_log($wb);
MapAlleles::setup($log,$wb);

my $release=$wb->get_wormbase_version;
my $acefile=( $outdir ? $outdir : $wb->acefiles ) . "/allele_mapping.WS$release.ace";

# DEBUG mode
if ($debug) {
    $maintainer = "$debug\@sanger.ac.uk";
    print "DEBUG \"$debug\"\n\n";
}

# get filtered arrayref of the alleles
my $alleles = $allele? MapAlleles::get_allele($allele) : MapAlleles::get_all_alleles();

# map them
my $mapped_alleles = MapAlleles::map($alleles);
undef $alleles;# could theoretically undef the alleles here 

# for other databases don't run through the GFF_SPLITs
&finish() if $database;

my $fh = new IO::File ">$acefile";
# create mapping Ace file
while( my($key,$allele)=each %$mapped_alleles){
	print $fh "Sequence : \"$allele->{clone}\"\nAllele $key $allele->{clone_start} $allele->{clone_stop}\n\n";
}

# get overlaps with genes
# gene_name->[allele_names,...]
my $genes=MapAlleles::get_genes($mapped_alleles);

# create the gene Ace file
my $inversegenes=MapAlleles::print_genes($genes,$fh);

# compare old<->new genes
MapAlleles::compare($mapped_alleles,$inversegenes);

# get overlaps with CDSs (intron,exon,coding_exon,cds)
# cds_name->allele_name->type
my $cds=MapAlleles::get_cds($mapped_alleles);
MapAlleles::print_cds($cds,$fh);

# get overlaps with Transcripts                        
# transcript_name->allele_name->type
my $utrs=MapAlleles::load_utr;
my $hit_utrs=MapAlleles::search_utr($mapped_alleles,$utrs);
MapAlleles::print_utr($hit_utrs,$fh);
 
# load to ace and close filehandle
MapAlleles::load_ace($fh,$acefile) unless $noload;

&finish();

# send the report
sub finish {
	if (MapAlleles::get_errors()){$log->mail( $maintainer, "BUILD REPORT: map_Alleles.pl ${\MapAlleles::get_errors()} ERRORS" )}
	else {$log->mail( $maintainer, 'BUILD REPORT: map_Alleles.pl')}
	exit 0;
}
