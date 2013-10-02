#!/usr/local/bin/perl5.8.0 -w
#
# further_BLAST_dump.pl
#
# script to copy files resulting from blast pipeline on farm-login to build directories
#
# Author: Chao-Kung CHen
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-10-02 12:16:21 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Log_files;
use Getopt::Long;
use File::Path;
use Storable;
use Carp;

my ($help, $debug, $test, $verbose, $store, $wormbase, $source_dir);
my $species;
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "sourcedir:s"=> \$source_dir,
	  );

$source_dir ||= $ENV{'PIPELINE'}.'/dumps';

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species
			     );
}

$species = $wormbase->species;
my $log = Log_files->make_build_log($wormbase);
my $farm_ace = "$source_dir/../ace_files"; 
my $target_dir = $wormbase->acefiles;
my $backup_dir = "$source_dir/BACKUP/";



my @files        = (
		    "SPECIES_best_blastp_hits",
		    "SPECIES_blastp.ace",
		    "SPECIES_blastx.ace",
		    "worm_ensembl_SPECIES_motif_info.ace",
		    "worm_ensembl_SPECIES_interpro_motif_info.ace",
		   );


# each separate species build will deal with copying its own files.
# the common files that are needed for the final build database will
# only be copied over by the elegans build

if ($species eq 'elegans'){
	# copy the non-species-specific files
	unlink("$source_dir/ensembl_protein_info.ace");

	$wormbase->run_command("cat $farm_ace/flybase.ace $farm_ace/yeast.ace $source_dir/ipi_hits.ace $source_dir/swissproteins.ace  > $source_dir/ensembl_protein_info.ace", $log);


	if (-e "$farm_ace/waba.ace") {
	$log->write_to("copying new version of waba.ace\n");
	$wormbase->run_command("cp -f $farm_ace/waba.ace $target_dir/waba.ace", $log);
	$wormbase->run_command("cp -f $farm_ace/waba.ace $backup_dir", $log);
	}
	if ( -e "$source_dir/ensembl_protein_info.ace" ) {
	$log->write_to("copying new version of ensembl_protein_info.ace\n");
	$wormbase->run_command("cp $source_dir/ensembl_protein_info.ace $target_dir/ensembl_protein_info.ace", $log);
	$wormbase->run_command("cp $source_dir/ensembl_protein_info.ace $backup_dir", $log);
	}
}

foreach my $f (@files){
	my $file = $f;		# don't want to change $f as it is a reference to the element in @files
	$file =~ s/SPECIES/$species/;
	if ( -e "$source_dir/$file" ) {
		$log->write_to("copying new version of $file\n");
		$wormbase->run_command("cp $source_dir/$file $target_dir/$file", $log);
		$wormbase->run_command("cp $source_dir/$file $backup_dir", $log);
	} else {
		# can give error as if this script is being run for this species it should be present!
		$log->error("$source_dir/$file does not exist\n");
	}
}



$log->mail;
exit(0);
