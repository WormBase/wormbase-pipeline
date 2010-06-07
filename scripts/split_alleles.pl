#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  split_alleles.pl
#
#        USAGE:  ./split_alleles.pl 
#
#  DESCRIPTION: iterates over 10 bins of alleles for mapping 
#
#       AUTHOR:   (Michael Paulini), <mh6@sanger.ac.uk>
#      COMPANY:  
#      VERSION:  $version:$
#      CREATED:  13/05/10 12:14:18 BST
#===============================================================================

use Ace;
use IO::File;
use Getopt::Long;
use lib  "$ENV{CVS_DIR}/Modules";
use map_Alleles;
use strict;

sub print_usage{
print  <<USAGE;
split_alleles.pl options:
	-debug USER_NAME             sets email address and debug mode
	-store FILE_NAME             use a Storable wormbase configuration file
	-outdir DIR_NAME             print allele_mapping_VERSION.ace to DIR_NAME
	-database DATABASE_DIRECTORY use a different AceDB
	-help                        print this message
	-test                        use the test database
	-species SPECIES_NAME        specify a non-elegans species
USAGE

exit 1;	
}

my ( $debug, $store, $outdir,$database,$help,$test,$species,$wb);

GetOptions(
    'species=s'=> \$species,
    'debug=s'  => \$debug,
    'store=s'  => \$store,
    'outdir=s' => \$outdir,
    'database=s'  => \$database,
    'help'        => \$help,
    'test'        => \$test,
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
MapAlleles::set_wb_log($log,$wb);

my $variations = MapAlleles::get_all_alleles();
my $binsize = int(@$variations / 10);
my $counter = 0;
my $bin=1;
my $of = new IO::File "$outdir/map_alleles.$bin",'w';
while (my $a = shift @$variations){
	if ($counter++ > $binsize){
		$bin++;
		$counter=1;
		$of->close;
                &mapAlleles($bin-1);
		$of->open("> $outdir/map_alleles.$bin");
	}
	print $of "$a\n";
}
$of->close;
&mapAlleles($bin-1);

sub mapAlleles {
	my ($lastBin) = @_;
	my $binfile="$outdir/map_alleles.$lastBin";
	$wb->run_script("map_alleles.pl -idfile $binfile",$log);
	unlink $binfile;
}
