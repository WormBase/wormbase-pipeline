#!/usr/local/bin/perl5.8.0 -w
#
# further_BLAST_dump.pl
#
# script to copy files resulting from blast pipeline on farm-login to build directories
#
# Author: Chao-Kung CHen
#
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2008-03-14 14:03:49 $

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
my $farm_base = '/lustre/work1/ensembl/wormpipe';
my $farm_ace = "$farm_base/ace_files"; 
my $source_dir    = "$farm_base/dumps";
my $target_dir = $wormbase->acefiles;
my $backup_dir = "$source_dir/BACKUP";

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
  if (( $file eq "waba.ace") && (-e "$farm_ace/waba.ace")){
    $log->write_to("copying new version of $file\n");
    $wormbase->run_command("cp $farm_ace/waba.ace $target_dir/waba.ace", $log);
    $wormbase->run_command("cp $farm_ace/waba.ace $backup_dir", $log);
  }
  elsif ( -e "$source_dir/$file" ) {
    $log->write_to("copying new version of $file\n");
    $wormbase->run_command("cp $source_dir/$file $target_dir/$file", $log);
    $wormbase->run_command("cp $source_dir/$file $backup_dir", $log);
  }
  else {
    $log->write_to($file." does not exist\n");
    $log->error unless $file eq 'waba.ace';
  }
}


$log->mail;
exit(0);
