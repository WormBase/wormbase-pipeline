#!/software/bin/perl -w
#
# prepare_primary_databases.pl
#
# Last edited by: $Author: gw3 $
# Last edited on: $Date: 2015-05-08 10:36:56 $

use strict;
my $scriptdir = $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use File::Path;


#################################################################################
# prepare primary databases                                                     #
#################################################################################
my %search_places = (
 elegans => [['citace', 'citace', 'citace*'],],
);

my ($test,$debug,$database, $store, $wormbase, $species);

GetOptions ( 
	    "test"       => \$test,
	    'debug:s'    => \$debug,
	    'database:s' => \$database,
	    'store:s'    => \$store,
	    'species:s'  => \$species
	   );

if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new(-debug   => $debug,
                            -test    => $test,
                            -organism => $species
                            );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);
$species = $wormbase->species;

my (%databases);

if ($species eq 'elegans') {
  &FTP_versions();
  &last_versions();
}


my (@errors, @options, @delete_from);

foreach my $primary (keys %databases){
  print "PRIMARY: $primary\n";
  next if (defined $database and ($database ne $primary));
  
  if( not $databases{$primary}->{new_date}) {
    push @errors, "Could not find latest data for $primary in FTP uploads\n";
  } elsif ($databases{$primary}->{new_date} le $databases{$primary}->{last_date}) {
    push @errors, "New version of $primary does not have newer date stamp than previous\n";
  } else {
    push @options, sprintf("-%s %s", $primary, $databases{$primary}->{new_date});
    push @delete_from, $primary;
    
    $log->write_to("For $primary, version on FTP site is newer than current - will delete and unpack\n");
  }
}

  
if (@errors) {
  foreach my $error (@errors) {
    $log->write_to("$error\n");
    print "ERROR: $error\n";
  }
  #$log->log_and_die("Did not find expected set of new data files; bailing\n");
}

# confirm unpack_db details and execute
foreach my $prim (@delete_from) {
  $log->write_to("Deleting old files from $prim\n");
  $wormbase->delete_files_from($wormbase->primary($prim),'*','+');
}
  
$log->write_to(" running unpack_db.pl @options\n\n");
$wormbase->run_script("unpack_db.pl @options", $log);


# New non FTP code for staging the species canonical databases
my (@refs,$ref);
if($species eq 'elegans') {
  push @refs, 'camace';
  push @refs, 'geneace';
}
else {
  push @refs, $species;
}

foreach $ref (@refs) {
my $primary_path = $wormbase->primary($ref);
$log->write_to("Transfering $primary_path to PRIMARIES\n");
print "Transfering $primary_path to PRIMARIES $ref\n";

# Check if the primary database is updated
my $test_file = "$primary_path/database/block1.wrm";
if (-e $test_file) {
  ($databases{$ref}->{last_date}) = $wormbase->find_file_last_modified($test_file);
}

# Delete the old files
$wormbase->delete_files_from($wormbase->primary($ref),'*','+');

# Transfer the new files to the build location
$wormbase->run_script("TransferDB.pl -start ".$wormbase->database($ref). " -end $primary_path -database -wspec", $log);

($databases{$ref}->{new_date}) = $wormbase->find_file_last_modified($test_file);
$wormbase->run_command("ln -sf ".$wormbase->autoace."/wspec/models.wrm ".$wormbase->primary($ref)."/wspec/models.wrm", $log);
my $test_file2 = $wormbase->primary($ref) . "/database/block1.wrm";
($databases{$ref}->{new_date}) = $wormbase->find_file_last_modified($test_file2);
}

my $log_msg = "\nDatabase versions used in the build (Previous => New):\n";

foreach my $db (keys %databases) {
  $log_msg .= sprintf("\t%s\t%s => %s\n", 
                      $db, 
                      exists($databases{$db}->{last_date}) ? $databases{$db}->{last_date} : "-",
                      exists($databases{$db}->{new_date}) ? $databases{$db}->{new_date} : "-");
}
$log->write_to("$log_msg\n\n");
my $record = join("/", $wormbase->reports, "Databases_used_in_build");
open my $rfh, ">$record";
foreach my $db (keys %databases) {
  if (exists $databases{$db}->{new_date}) {
    printf $rfh "%s\t%s\n", $db, $databases{$db}->{new_date};
  }
}
close($rfh);


$log->mail;
exit(0);


##################
# FTP_versions   #
##################

sub FTP_versions {
  $log->write_to("Getting FTP versions . . \n");

  my $ftp = $wormbase->ftp_upload;
  foreach my $location (@{$search_places{$species}}) {
    my ($db, $dir, $prefix_pat) = @$location;

    my @files = glob("$ftp/$dir/${prefix_pat}*.tar.gz");

    my @dates;
    foreach my $file (@files) {
      if ($file =~ /\/[^\/]+_(\d{4}\-\d{2}\-\d{2})\./) {
        push @dates, $1;
      }
    }
    if (@dates) {
      @dates = sort { $b cmp $a } @dates;
      $databases{$db}->{new_date} = $dates[0];
    } else {
      $databases{$db}->{new_date} = 0;
      $log->write_to("Could not find any version of $db in ftp_uploads\n");
    }
  }
}

#####################################################################################################################

sub last_versions {
  $log->write_to("Getting previous versions . . \n");

  foreach my $db (keys %databases) {
    my $primary = $wormbase->primary($db);
    my $unpack_dir = "$primary/temp_unpack_dir";
    if (-d $unpack_dir) {
      my @current_dates;
      foreach my $f (glob("$unpack_dir/${db}*.tar.gz")) {
        my ($date) = $f =~ /_(\d{4}\-\d{2}\-\d{2})/;
        push @current_dates, $date;
      }
      @current_dates = sort @current_dates;
      if (@current_dates) {
        $databases{$db}->{last_date} = $current_dates[-1];
      }
    }

    unless (defined $databases{$db}->{last_date} ) {
      $databases{$db}->{last_date} = '0000-00-00';
      $log->write_to("Could not find last update date for $db - assuming 0000-00-00\n");
    }
  }
}


__END__
