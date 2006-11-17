#!/usr/local/bin/perl5.8.0 -w
#
# further_BLAST_dump.pl
#
# script to copy files resulting from blast pipeline on farm-login to build directories
#
# Author: Chao-Kung CHen
#
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2006-11-17 15:10:58 $

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

my $source_dir    = "/lustre/work1/ensembl/wormpipe/dumps";
my $target_dir = $wormbase->acefiles;

my $farm_ace = glob("~wormpipe/ace_files") ;  # this is the only place where a path is specified outside of Wormbase.pm as cant access wormpipe and wormpub acefiles at same time
unlink("$source_dir/ensembl_protein_info.ace");
$wormbase->run_command("cat $farm_ace/flybase.ace $farm_ace/yeast.ace $source_dir/ipi_hits.ace $source_dir/swissproteins.ace $source_dir/tremblproteins.ace > $source_dir/ensembl_protein_info.ace", $log);

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
  if (( $file eq "waba.ace") && (-e "/nfs/acari/wormpipe/ace_files/waba.ace")){
    $log->write_to("scping new version of $file\n");
    $wormbase->run_command("scp acari:/nfs/acari/wormpipe/ace_files/waba.ace ${target_dir}/waba.ace", $log);
  }
  elsif ( -e "$source_dir/$file" ) {
    $log->write_to("scping new version of $file\n");
    $wormbase->run_command("scp farm-login:${source_dir}/${file} ${target_dir}/${file}", $log);
  }
  else {
    $log->write_to($file." does not exist\n");
    $log->error;
  }
}


$log->mail;
exit(0);
