#!/usr/local/bin/perl5.8.0 -w
#
# further_BLAST_dump.pl
#
# script to copy files resulting from blast pipeline on acari to wormsrv2
#
# Author: Chao-Kung CHen
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2005-12-16 11:18:55 $

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Log_files;
use Getopt::Long;
use File::Path;

my $debug;
my $test;

GetOptions(
	   "debug=s"     => \$debug,
	   "test"        => \$test,
	  );

my $log = Log_files->make_build_log($debug);

my $source_dir    = "/acari/work2a/wormpipe/dumps";
my $target_dir = "/wormsrv2/wormbase/ensembl_dumps";
$target_dir = glob("~wormpub/TEST_BUILD/wormbase/ensembl_dumps") if $test;

mkpath($target_dir,0) unless ( -e $target_dir);
mkpath("$target_dir/BACKUP",0) unless (-e "${target_dir}/BACKUP");

chdir "$target_dir";

my @files        = (
		    "ipi_hits.ace",
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
		    "swissproteins.ace",
		    "tremblproteins.ace",
		   );

foreach my $file (@files){
  if( ( $file eq "TRF.ace" or $file eq "repeat_homologies.ace") and !(-e $file) ){
    next;
  }
  if ( -e $file ) {
    $log->write_to("gzipping and backing up $file\n");
    &run_command("gzip -f ${target_dir}/$file > ${target_dir}/BACKUP/$file.gz");
  }
  else {
    $log->write_to($file." does not exist\n");
  }

  $log->write_to("scping new version of $file\n");
  &run_command("scp acari:${source_dir}/${file} ${target_dir}/${file}");
}


# Also move big ensembl_protein file to BACKUP dir before making new one. 
if (-e "$target_dir/ensembl_protein_info.ace") {
  $log->write_to("backing up protein info\n");
  &run_command("gzip -f $target_dir/ensembl_protein_info.ace > $target_dir/BACKUP/ensembl_protein_info.ace.gz");
}
$log->write_to("creating new ensembl_protein_info.ace\n");
&run_command("cat ipi_hits.ace flybase.ace yeast.ace swissproteins.ace tremblproteins.ace brigpep.ace > ensembl_protein_info.ace");

$log->mail;
exit(0);


##################################################################################################################


sub run_command{
  my $command = shift;
  my $status = system($command);
  if($status != 0){
    $log->write_to("ERROR: $command failed\n");
  }

  # for optional further testing by calling subroutine
  return($status);
}
