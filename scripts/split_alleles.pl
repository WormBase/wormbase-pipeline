#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  split_alleles.pl
#
#        USAGE:  ./split_alleles.pl 
#
#  DESCRIPTION: creates bins of alleles for mapping and submits them to LSF
#
#       AUTHOR:   (Michael Paulini), <mh6@sanger.ac.uk>
#      COMPANY:  
#      CREATED:  13/05/10 12:14:18 BST
#===============================================================================

use Ace;
use IO::File;
use Getopt::Long;
use FindBin qw($Bin);

use Modules::map_Alleles;

use lib '/software/worm/lib/site_perl/5.8.8/'; # that is where the Wormbaseified LSF module is
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;

use strict;

sub print_usage{
print  <<USAGE;
split_alleles.pl options:
	-debug USER_NAME             sets email address and debug mode
	-store FILE_NAME             use a Storable wormbase configuration file
	-outdir DIR_NAME             print allele_mapping_VERSION.ace to DIR_NAME
	-database DATABASE_DIRECTORY use a different AceDB
	-noload                      don't write back to AceDB
	-help                        print this message
	-test                        use the test database
	-species SPECIES_NAME        specify a non-elegans species
USAGE

exit 1;	
}

my $outdir = '/lustre/cbi4/scratch1/worm/tmp/map_allele_test';
my ( $debug, $store,$database,$help,$test,$species,$wb,$noload);

GetOptions(
    'species=s'=> \$species,
    'debug=s'  => \$debug,
    'store=s'  => \$store,
    'outdir=s' => \$outdir,
    'database=s'  => \$database,
    'help'        => \$help,
    'test'        => \$test,
    'noload'      => \$noload,
) or &print_usage();

&print_usage if $help;

my $maintainer = 'All';
if ($store) {
    $wb = Storable::retrieve($store) 
	    or croak("cannot restore wormbase from $store");
}
else { $wb = Wormbase->new( -debug => $debug, -test => $test, -organism => $species, -autoace => $database ) }

my $log = Log_files->make_build_log($wb);

if ($debug) {
    print "DEBUG \"$debug\"\n\n";
}

$database||=$wb->autoace();
MapAlleles::set_wb_log($log,$wb); # that is a bit crude, but makes $log available to the MapAlleles funtions

my $lsf = LSF::JobManager->new();
my @bsub_options =(-e => '/dev/null', -o => '/dev/null',-M => "3500000", -R => "\"select[mem>3500] rusage[mem=3500]\"");

my $variations = MapAlleles::get_all_alleles();

my $binsize = int(@$variations / 10 );
my $counter = 0;
my $bucket=1;

my $of = new IO::File "$outdir/map_alleles.$bucket",'w';

my @outfiles = ("$outdir/map_alleles.$bucket"); # to collect the names for later

while (my $a = shift @$variations){
	if ($counter++ > $binsize){
		$bucket++;
		$counter=1;
		$of->close;
		&mapAlleles($bucket-1);
		$of->open("> $outdir/map_alleles.$bucket");
		push(@outfiles,"$outdir/map_alleles.$bucket");
	}
	print $of "$a\n";
}
$of->close;
&mapAlleles($bucket);

$lsf->wait_all_children( history => 1 );

unless($noload){
    map {$wb->load_to_database($wb->autoace,$_,'map_alleles.pl',$log)} glob("$outdir/*.ace") ;
    map {unlink $_} @outfiles;
    map {unlink "$_.ace"} glob("$outdir/*.ace");
}

sub mapAlleles {
	my ($lastBin) = @_;
	my $binfile="$outdir/map_alleles.$lastBin";
	my $submitstring="/software/worm/perl_512/bin/perl $Bin/map_Alleles.pl -idfile $binfile -noload -outdir $outdir";
	$submitstring.=" -debug $debug" if $debug;
	$lsf->submit($submitstring);
}
