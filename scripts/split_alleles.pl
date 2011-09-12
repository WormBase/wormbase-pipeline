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
use lib "$Bin";

use Modules::map_Alleles;

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
	-noload                      do not write back to AceDB
	-help                        print this message
	-test                        use the test database
	-species SPECIES_NAME        specify a non-elegans species
USAGE

exit 1;	
}

my ( $debug, $store,$database,$help,$test,$species,$wb,$noload, $outdir);

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

if ($store) {
  $wb = Storable::retrieve($store) 
      or croak("cannot restore wormbase from $store");
}
else { 
  $wb = Wormbase->new( -debug => $debug, 
                       -test => $test, 
                       -organism => $species, 
                       -autoace => $database ); 
}

my $log = Log_files->make_build_log($wb);

if ($debug) {
    print "DEBUG \"$debug\"\n\n";
}

$outdir = $wormbase->autoace . "/TMP" if not defined $outdir;
mkdir $outdir if not -d $outdir;

$database = $wb->autoace() if not defined $database;
$species = $wb->species if not defined $species;

MapAlleles::set_wb_log($log,$wb); # that is a bit crude, but makes $log available to the MapAlleles funtions

my $lsf = LSF::JobManager->new();
my $variations = MapAlleles::get_all_allele_ids();
my $binsize = int(@$variations / 10 );


my (@bins, @id_files, @out_files);

while (my $a = shift @$variations){
  if (not @bins or scalar(@{$bins[-1]}) == $binsize) {
    push @bins, [];
  }

  push @{$bins[-1]}, $a;
}

for(my $bn=1; $bn <= @bins; $bn++) {
  my $list = $bins[$bn-1];

  my $id_file = "$outdir/map_alleles.$$.${species}.${bn}.txt";
  my $out_file = "$outdir/allele_mapping.$$.${species}.${bn}.ace";

  open(my $fh, ">$id_file") or 
      $log->log_and_die("Could not open IDfile $id_file for writing\n");  
  foreach my $nm (@$list) {
    print $fh "$nm\n";
  }
  close($fh);

  &mapAlleles($id_file, $out_file);

  push @id_files, $id_file;
  push @out_files, $out_file;
}

$lsf->wait_all_children( history => 1 );

if (not $noload) {
  foreach my $ace (@out_files) {
    if (-e $ace) {
      $wb->load_to_database($wb->autoace, $ace, 'map_alleles.pl',$log);
    } else {
      croak("Expected to find $ace, but its not there. Strange.\n");
    }
  }
  my $outfile = $wb->acefiles . "/map_alleles4geneace.ace";
  $wormbase->run_command("cat @out_files > $outfile", $log);

  map { unlink $_ } (@out_files, @id_files);
}


$log->mail();

sub mapAlleles {
  my ($id_file, $out_file) = @_;
  
  my $submitstring="/software/bin/perl $Bin/map_Alleles.pl -idfile $id_file -outfile $out_file -noload -species $species";
  $submitstring .= " -debug $debug" if $debug;
  $submitstring .= " -test"  if $test;

  my $job_name = "worm_".$species."_splitalleles";
  my @bsub_options =(-e => '/dev/null', 
                     -o => '/dev/null',
                     -M => "3500000", 
                     -R => "\"select[mem>3500] rusage[mem=3500]\"",
                     -J => $job_name);
  $lsf->submit(@bsub_options, $submitstring);
}
