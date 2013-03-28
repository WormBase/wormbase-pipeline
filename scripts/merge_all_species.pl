#/software/bin/perl -w
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-03-28 13:46:53 $

#################################################################################
# Variables                                                                     #
#################################################################################

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;

##############################
# command-line options       #
##############################

my $debug;      # Debug mode, verbose output to runner only
my $test;        # If set, script will use TEST_BUILD directory under ~wormpub
my $basedir;
my $store;
my $no_dump;

GetOptions (	"debug=s"    => \$debug,
		"test"       => \$test,
		"store:s"    => \$store,
                "nodump"     => \$no_dump,
	   	);
#this script is always run as elegans so no species option req.

my $wormbase;
if( $store ) {
    $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
    $wormbase = Wormbase->new( -debug   => $debug,
			       -test    => $test,
			       );
}

my $log = Log_files->make_build_log($wormbase);

my %accessors = $wormbase->species_accessors;
my $WS_name         = $wormbase->get_wormbase_version_name();

if (not $no_dump) {

  # move all the old MERGE files out of the way
  $log->write_to("\nRemove the old MERGE ace files . . .\n");
  foreach my $spDB (values %accessors) {
    my $dir = $spDB->acefiles."/MERGE/".$spDB->species."/";
    if (-e $dir) {
      my @files = glob("$dir/*.*");
      foreach my $file (@files) {
        if ($file !~ /\.old$/) {
          $wormbase->run_command("mv $file $file.old", $log);
        }
      }
    }
  }
  
  # dump out the files in parallel.
  my $lsf =  LSF::JobManager->new();
  $log->write_to("\nDumping acefile from . . .\n");
  foreach my $spDB (values %accessors) {
    $log->write_to("\t".$spDB->full_name('-short' => 1)."\n");
    my @bsub_options = (
      -o => '/dev/null',
      -M => 2000000,
      -R => sprintf("select[mem>%d] rusage[mem=%d]", 2000, 2000),
      -J => $spDB->species . "_make_acefiles_merge");
    my $cmd = $spDB->build_cmd("make_acefiles.pl -merge");
    
    $lsf->submit(@bsub_options, $cmd);
  }
  
  $lsf->wait_all_children( history => 1 );
  foreach my $job ($lsf->jobs) {
    if ($job->history->exit_status != 0) {
      $log->error(sprintf("make_acefiles job failed: %s\n", $job->history->command));
    }
  }
  if ($log->report_errors) {
    $log->log_and_die("At least one make_acefiles job died - exiting\n");
  }

  $log->write_to("\nFinished writing acefiles\n");
}


$log->write_to("\nDelete the homol_data files . . .\n");
foreach my $spDB (values %accessors) {
  my $file = sprintf("%s/MERGE/%s/%s_Homol_data.ace", $spDB->acefiles, $spDB->species, $spDB->species);
  if (-e $file) {
    $log->write_to("Removing file prior to loading: $file\n");
    unlink $file;
  }
}

$log->write_to("\nAbout to load . . .\n");
#and then load then one after another.
foreach my $spDB (values %accessors) {
  my $species = $spDB->species;
  $log->write_to("  Load from $species . . .\n");
  my @loaded;
  my $dir = $spDB->acefiles."/MERGE/".$spDB->species."/";
  next unless -e $dir;
  push(@loaded,$spDB->species);
  foreach my $file ( &read_dir($dir) ) {
    $log->write_to("    loading $file\n");
    $wormbase->load_to_database($wormbase->orgdb, $file, "merge_all_species", $log);
  }
  $log->write_to("  Loaded ".join(', ',@loaded)." in to ".$wormbase->orgdb."\n");
}

# look to see which BLAT files were loaded into the BUILD/species database and load them into autoace
$log->write_to("\nNow loading BLAT data . . .\n");
foreach my $spDB (values %accessors) {
  my $species = $spDB->species;
  $log->write_to("  Load BLAT from $species . . .\n");

  my $logfile = $spDB->orgdb."/database/log.wrm";
  open (LOG, "< $logfile") || die "Can't open $logfile : $!";
  my @acelist;
  while (my $line = <LOG>) {
    chomp $line;
    if ($line =~ /Parsing\s+file\s+(.+)/) {
      my $acefile = $1;
      # we only want the BLAT files and not duplicate files
      if ($acefile =~ /\/BLAT\// && !grep /$acefile/, @acelist) {
	if (-e $acefile) {
	  push @acelist, $acefile;
	} else {
	  $log->error("ERROR: The file $acefile does not exist!\n");
	}
      }
    }
  }
  close (LOG);

  foreach my $file ( @acelist ) {
    $log->write_to("    loading $file\n");
    $wormbase->load_to_database($wormbase->orgdb, $file, "merge_all_species", $log);
  }
}




$log->write_to("\nNow loading BLAST, protein and repeat data . . .\n");
foreach my $spDB (values %accessors) {
  my $species = $spDB->species;

  my @blastfiles = qw( SPECIES_blastp.ace worm_ensembl_SPECIES_interpro_motif_info.ace worm_ensembl_SPECIES_motif_info.ace repeat_homologies.ace inverted_repeats.ace pepace.ace);
  foreach my $f (@blastfiles){
    my $file = $f;		# don't use $f as it is a reference to the array element
    $file =~ s/SPECIES/$species/;
    if (-e $spDB->acefiles."/$file") {
      $log->write_to("    loading $file\n");
      $wormbase->load_to_database($wormbase->orgdb, $spDB->acefiles."/$file", "merge_all_species", $log);
    } 
    else {
      $log->error("ERROR: Can't find $file\n");
    }
  }
}


$log->write_to("\nNow loading briggsae TEC-RED homol data . . .\n");
my $briggsaeDB = Wormbase->new( 
			       -debug   => $wormbase->debug,
			       -test     => $wormbase->test,
			       -organism => 'briggsae'
			       );	
$wormbase->load_to_database($wormbase->orgdb, $briggsaeDB->acefiles."/misc_briggsae_TEC_RED_homol.ace", "merge_all_species", $log);
$wormbase->load_to_database($wormbase->orgdb, $briggsaeDB->acefiles."/misc_briggsae_TEC_RED_homol_data.ace", "merge_all_species", $log);


$log->write_to("\nNow loading briggsae BAC end data . . .\n");
#
# briggsae BAC end data
#
my $brigdb = $wormbase->database('briggsae');
my $bac_end_dir = "BAC_ENDS";
if (-d  "$brigdb/$bac_end_dir") {
  $log->write_to("\nLoading briggsae BAC ends from $brigdb/$bac_end_dir\n===========================\n");
  
  foreach my $be_file ("briggsae_BAC_ends.fasta",
                    "briggsae_BAC_end_homols.ace") {
    $log->write_to("\tLoading $be_file\n");
    $wormbase->load_to_database($wormbase->autoace,"$brigdb/$bac_end_dir/$be_file","BAC_ends", $log);
  }
} else {
  $log->write_to("\tWARNING : $brigdb/$bac_end_dir dir not found; will NOT load BAC_END data\n");
}

$log->mail;

exit;
	
###################################################################################################

sub read_dir {
  my $dir = shift;
  opendir (DIR,$dir) or $log->log_and_die("cant open directory $dir\n");
  $log->write_to("    Reading $dir\n");
  my @files = readdir DIR;
  my @to_load;
  foreach my $file ( @files ) {
    next if( $file eq '.' or $file eq '..');
    if( (-T $dir."/".$file) and substr($file,-3,3 ) eq "ace" ) {
      push (@to_load, "$dir"."$file");
    }
  }
  close DIR;
  return @to_load;
}
	
	
	
