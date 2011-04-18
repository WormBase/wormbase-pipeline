#!/software/binn/perl -w 
#
# initiate_build.pl
#
# Last edited by: $Author: klh $
# Last edited on: $Date: 2011-04-18 08:36:02 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use File::Copy;
use File::Spec;
use Storable;

my ($test,$debug,$database, $version, $species);
my ($store, $wormbase, $update);
GetOptions (
	    'test'       => \$test,
	    'debug:s'    => \$debug,
	    'database:s' => \$database,
	    'version:s'  => \$version,
	    'store:s'    => \$store,
	    'species:s'  => \$species,
	   );

if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    	=> $test,
			     -organism => $species,
			   );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$log->log_and_die( "version to build not specified\n") unless $version;

if ($version =~ /^WS(\d+)/) {
  $version = $1;
}
if ($version !~ /^\d\d\d$/) {
  $log->log_and_die("The version should be given as three digits\n");
}
$wormbase->version($version);

if (-e "${\$wormbase->autoace}/database") {
  $log->log_and_die( "There appears to still be data left over from the previous Build in autoace\nPlease check that finish_build.pl has completed.\n" );
}

# update the main scripts, autoace_contig and wgf for elegans only
if ($wormbase->species eq 'elegans') {
  $wormbase->run_command("cd $ENV{CVS_DIR}; cvs update .", $log)
      and $log->log_and_die("Failed to cvs update scripts dir; stopping\n");
  
  # update autoace_config
  $wormbase->run_command("cd ". $wormbase->basedir . "/autoace_config ; cvs update .", $log)
      and $log->log_and_die("Failed to cvs update autoace_config dir; stopping\n");
  
  #copy the genefinder files 
  $wormbase->run_command('cp -r '.$wormbase->build_data.'/wgf '.$wormbase->autoace.'/wgf', 'no_log');
}

#################################################################################
# initiate autoace build                                                        #
#################################################################################

## update CVS wspec, wquery and autoace_config from CVS
$wormbase->run_command("cd ".$wormbase->autoace."; cvs -d :pserver:cvsuser\@cvs.sanger.ac.uk:/cvsroot/ensembl checkout -d wspec -r WS${version} wormbase/wspec", $log)
    and $log->log_and_die("Failed to checkout wspec dir; not proceeding\n");
$wormbase->run_command("cd ".$wormbase->autoace."; cvs -d :pserver:cvsuser\@cvs.sanger.ac.uk:/cvsroot/ensembl checkout -d wquery wormbase/wquery", $log)
    and $log->log_and_die("Failed to checkout wquery directory; not proceeding\n");

#$wormbase->run_command("cd ".$wormbase->autoace."; cvs -d :ext:cvs.sanger.ac.uk:/cvsroot/ensembl checkout -d wspec -r WS${version} wormbase/wspec", $log);
#$wormbase->run_command("cd ".$wormbase->autoace."; cvs -d :ext:cvs.sanger.ac.uk:/cvsroot/ensembl checkout -d wquery wormbase/wquery", $log);

## update database.wrm using cvs
eval {
  my $dbwrm = $wormbase->autoace .  "/wspec/database.wrm";
  $wormbase->run_command("chmod u+w $dbwrm") 
      and die "Could not make database.wrm writable in order to update the version number\n";

  my $tmp_dbwrm = "/tmp/database.wrm";

  open my $new_fh, ">$tmp_dbwrm" or die "Could not open $tmp_dbwrm for writing\n";
  open my $old_fh, $dbwrm or die "Could not open $dbwrm for reading\n";

  $species = $wormbase->species;
  my $short_name = $wormbase->full_name('-short'=>1);
  
  while(<$old_fh>) {
    s/WS0/WS${version}/;
    s/species/${short_name}/;
    print $new_fh $_;
  };

  close($old_fh);
  close($new_fh);

  move($tmp_dbwrm, $dbwrm) or die "Could not move $tmp_dbwrm to $dbwrm\n";
};
$@ and $log->log_and_die("Failed to inject version number into database.wrm; stopping\n");

# Dump the sequence data from the species primary database being build.
$log->write_to("Dumping sequence data to file for ".$wormbase->species."\n");
$wormbase->run_script("dump_primary_seq_data.pl -organism $species", $log)
    and $log->log_and_die("Failed to successfully dump primary sequence data; stopping\n");

# Mask the sequences ready for BLATting
$log->write_to("Masking sequence data for ".$wormbase->species."\n");
$wormbase->run_script("BLAT_controller.pl -mask -qspecies $species", $log)
    and $log->log_and_die("Failed to successfully mask sequence data; stopping\n");

# add lines to the logfile
my $msg = "Updated ".$wormbase->species." version number to WS".$wormbase->version."\n";
$msg .= "You are ready to build another WormBase release\n";
$msg .= "Please tell camace and geneace curators to update their database to use the new models!!!\n\n";

$log->write_to($msg);
$log->mail;
exit(0);

__END__
