#!/usr/local/bin/perl5.8.0 -w
#
# autoace_builder.pl
# 
# based on the original autoace_minder.pl
#
# Usage : autoace_builder.pl [-options]
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2005-12-19 17:04:58 $

my $script_dir =  $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'};

use strict;
use Wormbase;
use Getopt::Long;
use File::Copy;
use Coords_converter;
use Log_files;

my ($debug, $test, $database);
my ($initiate, $prepare_databases, $acefile, $build, $first_dumps);
my ($make_wormpep, $finish_wormpep);
my ($run_blat, $finish_blat);
my ($gff_dump, $processGFF, $gff_split);
my $gene_span;
my $load;

GetOptions (
	    'debug:s'     => \$debug,
	    'test'        => \$test,
	    'database:s'  => \$database,
	    'initiate:s'  => \$initiate,
	    'prepare'     => \$prepare_databases,
	    'acefiles'    => \$acefile,
	    'build'       => \$build,
	    'first_dumps' => \$first_dumps,
	    'make_wormpep'=> \$make_wormpep,
	    'finish_wormpep'=> \$finish_wormpep,
	    'gff_dump:s'  => \$gff_dump,
	    'processGFF:s'=> \$processGFF,
	    'gff_split'   => \$gff_split,
	    'gene_span'   => \$gene_span,
	    'load'        => \$load,
	    'run_blat'    => \$run_blat,
	    'finish_blat' => \$finish_blat
	   );

my $wormbase = Wormbase->new(
			     -test    => $test,
			     -debug   => $debug,
			     -version => $initiate
			    );

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$wormbase->run_script("initiate_build.pl", $log)            if defined($initiate);
$wormbase->run_script('prepare_primary_databases.pl', $log) if $prepare_databases;
$wormbase->run_script('make_acefiles.pl', $log)             if $acefile;
$wormbase->run_script('make_autoace.pl', $log)              if $build;
$wormbase->run_script("build_dumpGFF.pl -stage $gff_dump",$log) if $gff_dump;

#//--------------------------- batch job submission -------------------------//
$wormbase->run_script("processGFF.pl -$processGFF",$log)        if $processGFF;  #clone_acc
&first_dumps if $first_dumps;
$wormbase->run_script('make_wormpep.pl -initial', $log)         if $make_wormpep;

#########   BLAT  ############
$wormbase->run_script('BLAT_controller.pl -mask -dump -run',$log) if $run_blat;

#//--------------------------- batch job submission -------------------------//

$wormbase->run_script('BLAT_controller.pl -process -postprocess -ace -load', $log) if $finish_blat;

######## WBGene spans ########
$wormbase->run_script('WBGene_span.pl -prepare', $log) if $gene_span;

if( $load ) {
  my $file = shift;
  my $tsuser = shift;

  $log->("loading $file to $database\n");
  $log->("\ttsuser = $tsuser\n\n");
  $wormbase->load_to_database($database, $file, $tsuser) if( $file );
}

$log->mail;
exit(0);

############################
#       SUBROUTINES        #
############################

sub first_dumps
  {
    $wormbase->run_script("chromosome_dump.pl --dna --composition", $log);
    $wormbase->run_script("make_agp_file.pl", $log);
    $wormbase->run_script("agp2dna.pl", $log);

    my $agp_errors = 0;

    my @chrom = qw( I II III IV V X);
    foreach my $chrom (@chrom) {
      open (AGP, ">".$wormbase->autoace."/yellow_brick_road/CHROMOSOME_${chrom}.agp_seq.log") or die "Couldn't open agp file : $!";
      while (<AGP>) {
	$agp_errors++ if (/ERROR/);
      }
      close(AGP);
    }
    $log->write_to("ERRORS ( $agp_errors ) in agp file\n");
  }
