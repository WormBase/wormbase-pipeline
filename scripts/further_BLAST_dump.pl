#!/usr/local/bin/perl5.8.0 -w
#
# further_BLAST_dump.pl
#
# script to copy files resulting from blast pipeline on acari to build directories
#
# Author: Chao-Kung CHen
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2006-02-17 11:32:47 $

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

my $source_dir    = "acari:/acari/work2a/wormpipe/dumps";
my $target_dir = $wormbase->acefiles;


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
		   );

foreach my $file (@files){
  if ( -e $file ) {
    $log->write_to("scping new version of $file\n");
    &run_command("scp ${source_dir}/${file} ${target_dir}/${file}");
  }
  else {
    $log->write_to($file." does not exist\n");
  }
}

#will do this as part of blast pipeline and only copy single file
#$log->write_to("creating new ensembl_protein_info.ace\n");
#&run_command("cat ipi_hits.ace flybase.ace yeast.ace swissproteins.ace tremblproteins.ace brigpep.ace > ensembl_protein_info.ace");

my $farm_ace = glob("~wormpipe/acefiles") ;  # this is the only place where a path is specified outside of Wormbase.pm as cant access wormpipe and wormpub acefiles at same time

$wormbase->run_command("cat $farm_ace/flybase.ace $farm_ace/yeast.ace $farm_ace/ipi_human.ace $farm_ace/swissproteins.ace $farm_ace/tremblproteins.ace >! $farm_ace/ensembl_protein_info.ace");
&run_command("scp $farm_ace/ensembl_protein_info.ace  ${target_dir}/ensembl_protein_info.ace");

#will do this as part of blast pipeline and only copy single file
#$log->write_to("creating new ensembl_protein_info.ace\n");
#&run_command("cat ipi_hits.ace flybase.ace yeast.ace swissproteins.ace tremblproteins.ace brigpep.ace > ensembl_protein_info.ace");

$log->mail;
exit(0);
