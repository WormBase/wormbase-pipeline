#!/usr/local/ensembl/bin/perl -w
#
# copy_files_to_acari.pl
#
# dl
#
# copy the latest dna and wormpep files from wormsrv2 to acari

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use vars qw($opt_c $opt_w);
# $opt_w copy wormpep data
# $opt_c copy chromosomes
use Log_files;

my ($wormpep, $chroms, $test, $debug, $store);
GetOptions (
	    'wormpep:i' => \$wormpep,
	    'chrom'     => \$chroms,
	    "debug=s"   => \$debug,
	    "test"      => \$test,
	    "store:s"   => \$store
	   );

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( 'debug'   => $debug,
                             'test'    => $test,
			     );
}

my $dir = glob("~wormpipe/BlastDB");
my $log = Log_files->make_build_log($wormbase);

my @CHROMOSOME = ('I','II','III','IV','V','X');
if ( $chroms )
  {
    print "-------------------------------------------------------\nCopying new Chromosomal DNA\n\n";
    $log->write_to("Removing old DNA files . . . ");
    $wormbase->run_command( "rm -f $dir/*dna") ;
    
    foreach my $chrom (@CHROMOSOME) {
      print "about to get Chromosome $chrom\n";
      $wormbase->run_command( "scp ".$wormbase->chromosomes."/CHROMOSOME_$chrom.dna $dir/" );
      $log->write_to("chromosome $chrom complete\n");
    }

    $log->write_to("Finished copying DNA to acari\n\n-------------------------------------------------------\n");
  }

if( $wormpep )
  {
    our $WS_version = $wormbase->get_wormbase_version;
    
    my $wp_file = "wormpep" . $WS_version . ".pep";
    my $wp_old  = "wormpep" . ($WS_version - 1) . ".pep";
    
    $wormbase->run_command("rcp ".$wormbase->wormpep."/wormpep$WS_version $dir/wormpep${WS_version}.pep");

    
    $log->write_to("Removed old version of wormpep . . ");
    $wormbase->run_command("rm -f $dir/${wp_old}*") ;

    $log->write_to("Wormpep copied and ready for action\n\n-------------------------------------------------------\n");
    print "finished\n";
  }

$log->mail;
exit(0);
