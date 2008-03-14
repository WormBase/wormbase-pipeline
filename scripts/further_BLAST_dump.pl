#!/usr/local/bin/perl5.8.0 -w
#
# further_BLAST_dump.pl
#
# script to copy files resulting from blast pipeline on farm-login to build directories
#
# Author: Chao-Kung CHen
#
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2008-03-14 13:33:48 $

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
my $backup_dir = "$source_dir/BACKUP";

my $farm_ace = glob("~wormpipe/ace_files") ;  # this is the only place where a path is specified outside of Wormbase.pm as cant access wormpipe and wormpub acefiles at same time
unlink("$source_dir/ensembl_protein_info.ace");
$wormbase->run_command("cat $farm_ace/flybase.ace $farm_ace/yeast.ace $source_dir/ipi_hits.ace $source_dir/swissproteins.ace $source_dir/tremblproteins.ace > $source_dir/ensembl_protein_info.ace", $log);

my @files        = (
		    "worm_pep_best_blastp_hits",
		    "repeat_homologies.ace",
		    "worm_pep_blastp.ace",
		    "worm_dna_blastx.ace",
		    "worm_ensembl_elegans_motif_info.ace",
		    "worm_ensembl_elegans_interpro_motif_info.ace",
		    "ensembl_protein_info.ace",
		    "waba.ace",
		   );

foreach my $file (@files){
  if (( $file eq "waba.ace") && (-e "/nfs/acari/wormpipe/ace_files/waba.ace")){
    $log->write_to("scping new version of $file\n");
    $wormbase->run_command("scp farm-login:/nfs/acari/wormpipe/ace_files/waba.ace ${target_dir}/waba.ace", $log);
    $wormbase->run_command("cp /nfs/acari/wormpipe/ace_files/waba.ace $backup_dir", $log);
  }
  elsif ( -e "$source_dir/$file" ) {
    $log->write_to("scping new version of $file\n");
    $wormbase->run_command("scp farm-login:${source_dir}/${file} ${target_dir}/${file}", $log);
    $wormbase->run_command("cp ${source_dir}/${file} $backup_dir", $log);
  }
  else {
    $log->write_to($file." does not exist\n");
    $log->error;
  }
}


$log->mail;
exit(0);
