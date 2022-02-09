#/software/bin/perl -w
#
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2013-08-13 12:45:19 $

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
my ($do_not_load, %do_not_load, $dry_run);

GetOptions ( "debug=s"     => \$debug,
	     "test"        => \$test,
	     "store:s"     => \$store,
             "nodump"      => \$no_dump,
             "donotload=s" => \$do_not_load,
             "dryrun"      => \$dry_run,
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

if (defined $do_not_load) {
  open(my $fh, $do_not_load) or $log->log_and_die("Could not open file $do_not_load\n");
  while(<$fh>) {
    /^(\S+)/ and $do_not_load{$1} = 1;
  }
}

# look to see which BLAT files were loaded into the BUILD/species database and load them into autoace
$log->write_to("\nLoading missed BLAT data . . .\n");
foreach my $spDB (values %accessors) {
    my $species = $spDB->species;
    $log->write_to("  Load BLAT from $species . . .\n");
    
    my $logfile = $spDB->orgdb."/database/log.wrm";
    open (my $log_fh, $logfile) || die "Can't open $logfile : $!";
    my @acelist;
    while (my $line = <$log_fh>) {
	chomp $line;
	if ($line =~ /Parsing\s+file\s+(.+)/) {
	    my $acefile = $1;
	    # Logs for species that haven't yet been buid on codon still contain old file paths
	    # We only want to load these files in this script for WS284 as all other files have been loaded
	    if ($acefile =~ /^\/nfs\/panda\/ensemblgenomes\/wormbase(.+)$/) {
		$acefile = $ENV{'BUILD_HOME'} . $1;
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
    }
    close ($log_fh);

    foreach my $file ( @acelist ) {
	&load( $file );
    }
}

$log->mail;

exit;
	
###################################################################################################

	
sub load {
  my ($file, $tsuser) = @_;

  $tsuser = "merge_all_species" if not defined $tsuser;
  
  if (not exists $do_not_load{$file}) {
    $log->write_to("   loading $file\n");
    if (not $dry_run) {
      $wormbase->load_to_database($wormbase->orgdb, $file, "merge_all_species", $log);
    }
    return $file;
  } else {
    $log->write_to("   NOT loading $file\n");
    return ();
  }
}
