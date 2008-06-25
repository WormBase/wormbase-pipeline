#!/usr/local/bin/perl5.8.0 -w
#
# further_BLAST_dump.pl
#
# script to copy files resulting from blast pipeline on farm-login to build directories
#
# Author: Chao-Kung CHen
#
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2008-06-25 16:06:02 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Log_files;
use Getopt::Long;
use File::Path;
use Storable;
use Carp;

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

my $species = $wormbase->species;
my $log = Log_files->make_build_log($wormbase);
my $farm_base = '/lustre/work1/ensembl/wormpipe';
my $farm_ace = "$farm_base/ace_files"; 
my $source_dir    = "$farm_base/dumps";
my $target_dir = $wormbase->acefiles;
my $backup_dir = "$source_dir/BACKUP";

unlink("$source_dir/ensembl_protein_info.ace");

$wormbase->run_command("cat $farm_ace/flybase.ace $farm_ace/yeast.ace $source_dir/ipi_hits.ace $source_dir/swissproteins.ace $source_dir/tremblproteins.ace > $source_dir/ensembl_protein_info.ace", $log);

my @files        = (
		    "SPECIES_best_blastp_hits",
		    "SPECIES_blastp.ace",
		    "SPECIES_blastx.ace",
		    "worm_ensembl_SPECIES_motif_info.ace",
		    "worm_ensembl_SPECIES_interpro_motif_info.ace",
		   );


# copy the non-species-specific files
if (-e "$farm_ace/waba.ace") {
  $log->write_to("copying new version of waba.ace\n");
  $wormbase->run_command("cp $farm_ace/waba.ace $target_dir/waba.ace", $log);
  $wormbase->run_command("cp $farm_ace/waba.ace $backup_dir", $log);
}
if ( -e "$source_dir/ensembl_protein_info.ace" ) {
  $log->write_to("copying new version of ensembl_protein_info.ace\n");
  $wormbase->run_command("cp $source_dir/ensembl_protein_info.ace $target_dir/ensembl_protein_info.ace", $log);
  $wormbase->run_command("cp $source_dir/ensembl_protein_info.ace $backup_dir", $log);
}


# copy these file for each tier 2 species
my %accessors = ($wormbase->species_accessors);
$accessors{elegans} = $wormbase;
foreach my $species (keys %accessors) {

  foreach my $f (@files){
    my $file = $f;		# don't want to change $f as it is a reference to the element in @files
    $file =~ s/SPECIES/$species/;
    if ( -e "$source_dir/$file" ) {
      $log->write_to("copying new version of $file\n");
      $wormbase->run_command("cp $source_dir/$file $target_dir/$file", $log);
      $wormbase->run_command("cp $source_dir/$file $backup_dir", $log);
    } else {
      warn "$source_dir/$file does not exist\n";
    }
  }
}


$log->mail;
exit(0);
