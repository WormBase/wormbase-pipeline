#!/software/bin/perl -w
#
#
# dump_protein_features.pl
# 
# Dumps protein features from an Ensembl database as Ace. 
# Calls on a number of other scripts
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2013-05-08 09:49:52 $      

use strict;                                      

use FindBin qw($Bin);
use lib "${Bin}/..";

use Getopt::Long;
use Carp;
use Storable;

use Wormbase;
use Log_files;


######################################
# variables and command-line options # 
######################################

# 
# some scripts spawned by this one do a sort, and these
# are the core parameters required for those sorts
#
$ENV{SORT_OPTS} = "-k2,2 -k8,8n -k10,10nr";


my $Blast_scripts = "BLAST_scripts";
my $WORMPIPE_DIR = $ENV{'PIPELINE'};
my $SORT_DUMP_DIR = "$WORMPIPE_DIR/sort_dump";

croak("The target directory $WORMPIPE_DIR must exist") if not -d $WORMPIPE_DIR;

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($checkonly, $species);

GetOptions("help"       => \$help,
           "debug=s"    => \$debug,
           "test"       => \$test,
           "store:s"    => \$store,
           "check_only" => \$checkonly, #allows you to run just the checking part of the script
	   "species:s"  => \$species,
           );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species,
			   );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

##########################
# MAIN BODY OF SCRIPT
##########################

my $version = $wormbase->get_wormbase_version;
my $dumpdir = "$WORMPIPE_DIR/dumps";
my $acedir  = $wormbase->acefiles;
$species = $wormbase->species;

if (!$checkonly) {
  $log->write_to("  Running dump.pl . . .\n");
  
  my $sort_file_pre = join(".", $species, "junk");
  my $sort_file_out = "$species.srt";

  my $dump_one_cmd = "${Bin}/dump_one_new.pl";
  my $dump_all_cmd = "perl ${Bin}/dump.pl -db worm_ensembl_${species} -dump_script $dump_one_cmd -dumploc $SORT_DUMP_DIR -prefix $sort_file_pre";

  $wormbase->run_command($dump_all_cmd, $log)
      and $log->log_and_die("Failed to successfully run command - stopping ($dump_all_cmd)\n");

  $log->write_to("  Sorting $species.srt . . .\n");

  my @sort_input = glob("$SORT_DUMP_DIR/${sort_file_pre}.*.srt");
  my $sort_cmd = "sort -m -S 2G -T $SORT_DUMP_DIR " . $ENV{SORT_OPTS} . " -o $SORT_DUMP_DIR/${sort_file_out} @sort_input";
  $wormbase->run_command($sort_cmd, $log)
      and $log->log_and_die("Failed to successfully run command - stopping ($sort_cmd)\n");;
  unlink @sort_input;

  $log->write_to("  Running dump_blastp_from_file.pl -file $SORT_DUMP_DIR/${sort_file_out} -matches -database worm_${species} -dumpdir $dumpdir\n");
  $wormbase->run_script("${Blast_scripts}/dump_blastp_from_file.pl -file $SORT_DUMP_DIR/${sort_file_out} -matches -database worm_${species} -dumpdir $dumpdir", $log);

  $log->write_to("  Running dump_motif.pl -database worm_ensembl_${species} -dumpdir $dumpdir\n");
  $wormbase->run_script("${Blast_scripts}/dump_motif.pl -database worm_ensembl_${species} -dumpdir $dumpdir", $log);

  $log->write_to("  Running dump_interpro_motif.pl -database worm_ensembl_${species}\n");
  $wormbase->run_script("${Blast_scripts}/dump_interpro_motif.pl -database worm_ensembl_${species}", $log);

  $log->write_to("  Running dump_repeats.pl  -database worm_ensembl_${species} -dumpdir $acedir\n");
  $wormbase->run_script("${Blast_scripts}/dump_repeats.pl  -database worm_ensembl_${species} -dumpdir $acedir", $log);

}


##################
# Check the files
################## 
$wormbase->check_file($dumpdir."/${species}_blastp.ace", 
                      $log,
                      minsize => 1000000,
                      );

$wormbase->check_file($dumpdir."/worm_ensembl_${species}_motif_info.ace", 
                      $log,
                      minsize => 1000000,
                      );

$wormbase->check_file($dumpdir."/worm_ensembl_${species}_interpro_motif_info.ace", 
                      $log,
                      minsize => 1000000,
                      );

$wormbase->check_file($acedir."/repeat_homologies.ace", 
                      $log,
                      minsize => 1000000,
                      );


$log->mail();
print "\nFinished.\n";
exit(0);


