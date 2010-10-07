#!/usr/bin/env perl
# Copyright (c) 2007, Genome Research Ltd (GRL).
# # All rights reserved.
# # Author: Michael Han <mh6@sanger.ac.uks>
# #
# # Redistribution and use in source and binary forms, with or without
# # modification, are permitted provided that the following conditions
# # are met:
# #     * Redistributions of source code must retain the above copyright
# #       notice, this list of conditions and the following disclaimer.
# #     * Redistributions in binary form must reproduce the above
# #       copyright notice, this list of conditions and the following
# #       disclaimer in the documentation and/or other materials
# #       provided with the distribution.
# #     * Neither the name of the <organization> nor the
# #       names of its contributors may be used to endorse or promote
# #       products derived from this software without specific prior
# #       written permission.
# #
# # THIS SOFTWARE IS PROVIDED BY GRL ``AS IS'' AND ANY EXPRESS OR
# # IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# # WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# # ARE DISCLAIMED.  IN NO EVENT SHALL GRL BE LIABLE FOR ANY DIRECT,
# # INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# # (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# # SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# # HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# # STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# # ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
# # OF THE POSSIBILITY OF SUCH DAMAGE.

use strict;
use lib $ENV{'CVS_DIR'};
use lib "$ENV{CVS_DIR}/Modules";
use map_Alleles;
use Wormbase;
use Getopt::Long;
use IO::File;              

sub print_usage{
print  <<USAGE;
map_Allele.pl options:
	-debug USER_NAME    sets email address and debug mode
	-store FILE_NAME    use a Storable wormbase configuration file
	-outdir DIR_NAME    print allele_mapping_VERSION.ace to DIR_NAME
	-allele ALLELE_NAME check only ALLELE_NAME instead of all
	-noload             dont update AceDB
	-noupdate           same as -noload
	-database			DATABASE_DIRECTORY use a different AceDB
	-weak_checks        relax sequence sanity checks
	-help               print this message
	-test               use the test database
	-idfile             use ids from an input file (one id per line)
	-species SPECIES_NAME use species as reference
USAGE

exit 1;	
}

my ( $debug, $species, $store, $outdir,$allele ,$noload,$database,$weak_checks,$help,$test,$idfile);

GetOptions(
    'species=s'=> \$species,
    'debug=s'  => \$debug,
    'store=s'  => \$store,
    'outdir=s' => \$outdir,
    'allele=s' => \$allele,
    'noload'   => \$noload,
    'noupdate' => \$noload,
    'database=s'  => \$database,
    'weak_checks' => \$weak_checks,
    'help'        => \$help,
    'test'        => \$test,
    'idfile=s'    => \$idfile,
) or &print_usage();

&print_usage if $help;

my $maintainer = 'All';
my $wb;

if ($store) {
    $wb = Storable::retrieve($store)
      or croak("cannot restore wormbase from $store");
}
else { $wb = Wormbase->new( -debug => $debug, -test => $test, -organism => $species, -autoace => $database ) }

my $log = Log_files->make_build_log($wb);
MapAlleles::setup($log,$wb) unless $database;
MapAlleles::set_wb_log($log,$wb,$weak_checks) if $database;

my $release=$wb->get_wormbase_version;
my $acefile=( $outdir ? $outdir : $wb->acefiles ) . "/allele_mapping.WS$release.$$.ace";

if ($debug) {
    print "DEBUG \"$debug\"\n\n";
}

# get filtered arrayref of the alleles
my $alleles;
if ($allele){
       $alleles= MapAlleles::get_allele($allele);
}elsif ($idfile){
       $alleles= MapAlleles::get_alleles_fromFile($idfile);
}else{
       $alleles= MapAlleles::get_all_alleles();
}

# map them
my $mapped_alleles = MapAlleles::map($alleles);
undef $alleles;# could theoretically undef the alleles here 


# for other databases don't run through the GFF_SPLITs
&finish() if $database;

my $fh = new IO::File ">$acefile" || die($!);
# create mapping Ace file
while( my($key,$allele)=each %$mapped_alleles){
	print $fh "Sequence : \"$allele->{clone}\"\nAllele $key $allele->{clone_start} $allele->{clone_stop}\n\n";
	print $fh "Variation : \"$key\"\nSequence  $allele->{clone}\n\n";
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
$utrs = $hit_utrs = undef; # cleanup memory

# get overlaps with Pseudogenes                        
# pseudogene_name->allele_name->1
my $pgenes=MapAlleles::load_pseudogenes;
my $hit_pgenes=MapAlleles::search_pseudogenes($mapped_alleles,$pgenes);
MapAlleles::print_pseudogenes($hit_pgenes,$fh);
$pgenes = $hit_pgenes = undef; # cleanup memory

# get overlaps with non-coding transcripts                        
# transcript_name->allele_name->1
my $nc_rnas=MapAlleles::load_ncrnas;
my $hit_ncrnas=MapAlleles::search_ncrnas($mapped_alleles,$nc_rnas);
MapAlleles::print_ncrnas($hit_ncrnas,$fh);
$hit_ncrnas = $hit_ncrnas = undef; # cleanup memory
 
# load to ace and close filehandle
MapAlleles::load_ace($fh,$acefile) unless $noload;

&finish();

# send the report
sub finish {
	if (MapAlleles::get_errors()){$log->mail( $maintainer, "BUILD REPORT: map_Alleles.pl ${\MapAlleles::get_errors()} ERRORS" )}
	else {$log->mail( $maintainer, 'BUILD REPORT: map_Alleles.pl')}
	exit 0;
}
