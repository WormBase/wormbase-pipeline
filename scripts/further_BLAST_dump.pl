#!/usr/local/bin/perl5.8.0 -w
#
# further_BLAST_dump.pl
#
# script to copy files resulting from blast pipeline on acari to build directories
#
# Author: Chao-Kung CHen
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2006-03-02 17:50:55 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Log_files;
use Getopt::Long;
use File::Path;
use Storable;

my ($help, $debug, $test, $verbose, $store, $wormbase);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"      => \$store,
	  );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

my $log = Log_files->make_build_log($wormbase);

my $source_dir    = "/acari/work2a/wormpipe/dumps";
my $target_dir = $wormbase->acefiles;

my $farm_ace = glob("~wormpipe/ace_files") ;  # this is the only place where a path is specified outside of Wormbase.pm as cant access wormpipe and wormpub acefiles at same time
unlink("$source_dir/ensembl_protein_info.ace");
$wormbase->run_command("cat $farm_ace/flybase.ace $farm_ace/yeast.ace $source_dir/ipi_hits.ace $source_dir/swissproteins.ace $source_dir/tremblproteins.ace > $source_dir/ensembl_protein_info.ace");

my @files        = (
		    "worm_pep_best_blastp_hits",
		    "worm_brigpep_best_blastp_hits",
		    "TRF.ace",
		    "repeat_homologies.ace",
		    "worm_pep_blastp.ace",
		    "worm_brigpep_blastp.ace",
		    "worm_dna_blastx.ace",
		    "worm_pep_motif_info.ace",
		    "worm_brigpep_motif_info.ace",
		    "worm_pep_interpro_motif_info.ace",
		    "worm_brigpep_interpro_motif_info.ace",
		    "ensembl_protein_info.ace",
		    "waba.ace",
		   );

foreach my $file (@files){
  if ( -e $file ) {
    $log->write_to("scping new version of $file\n");
    $wormbase->run_command("scp acari:${source_dir}/${file} ${target_dir}/${file}");
  }
  else {
    $log->write_to($file." does not exist\n");
  }
}

$log->mail;
exit(0);
