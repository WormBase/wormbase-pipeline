#/software/bin/perl -w
#
# make_interpro.pl
# 
# by Gary Williams
#
# This detects whether there is a new version of the Interpro databases available
# It copies over the latest file, nupacks it into place and runs the InterPro
# indexing program on it.
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2008-04-17 08:28:06 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use DBI;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($species, $nodownload, $dbpass);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species, # for testing puposes, comma-delimted list of species to delete old results from, all species are done if none are specified
	    "nodownload" => \$nodownload,  # for testing purposes, inhibit the downloading and indexing of the databases even if they should be downloaded
	    "dbpass:s"   => \$dbpass, # the worm_ensembl_* mysql databases password
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

# if no species are specified, then assume that we need to clean out
# the old results from the database of all species
# We may wish to specify the name of a species not in this list
# for testing purposes, e.g. 'brugia'
my @species_list = qw(
		      elegans
		      briggsae
		      remanei
		      );
my @species;			# species databases to remove old results from
if (! $species) {
  @species = @species_list;
} else {
  @species = split /[\,\s+]/, $species;
}


##########################
# MAIN BODY OF SCRIPT
##########################


# get the current database release values
my ($current_main, $current_panther) = &get_current_versions;

# get the latest version numbers from the web site
my ($latest_main, $latest_panther) = &get_versions();

# if the main set of databases has changed, get them
my ($get_main, $get_panther);
$log->write_to("Latest versions:  Main='$latest_main', Panther='$latest_panther'\n");
$log->write_to("Current versions: Main='$current_main', Panther='$current_panther'\n");

if ($latest_main eq $current_main) {
  $log->write_to("No new release required for Main InterPro database\n");
} else {
  $log->write_to("Getting new release for Main InterPro database\n");
  $get_main = 1;
}

if ($latest_panther eq $current_panther) {
  $log->write_to("No new release required for Panther InterPro database\n");
} else {
  $log->write_to("Getting new release for Panther InterPro database\n");
  $get_panther=1;
}

# if the main databases changed, get them
if (!$nodownload) {
  if ($get_main) {
    $log->write_to("Downloading the main Interpro databases\n");
    &get_file('main', $latest_main);
  }

# if the panther database changed, get it
  if ($get_panther) {
    $log->write_to("Downloading the Panther database\n");
    &get_file('panther', $latest_panther);
  }
} else {
  $log->write_to("Not downloading any database files because -nodownload is specified\n");
}

# index the databases if they have changed
if ($get_main || $get_panther) {
  if (!$nodownload) {
    $log->write_to("Indexing the databases for interproscan\n");
    $wormbase->run_command("perl /software/worm/iprscan/bin/index_data.pl -p /data/blastdb/Worms/interpro_scan/iprscan/data/ -inx -bin -v", $log);
  }

# for each species, delete the results of previous runs of the
# analysis pipeline for interpro
  foreach my $species_to_delete (@species) {
    if ($test) {
      $log->write_to("In test mode - the Interpro results for $species_to_delete will not be changed\n");
    } else {
      $log->write_to("Deleting old Interpro analysis results\n");
      &delete_results($get_main, $get_panther, $species_to_delete);
    }
  }

# save the version number of the new database
  $log->write_to("Saving record of latest database versions\n");
  &save_versions($latest_main, $latest_panther);

}

$log->mail();
exit(0);






##############################################################
#
# Subroutines
#
##############################################################

##########################################
# get the current versions
#my ($current_main, $current_panther) = &get_current_versions;

sub get_current_versions {
  my ($current_main, $current_panther);

  my $file = $wormbase->compare . "/interpro_versions.dat";

  open (FILE, "< $file") || $log->log_and_die("Can't open $file\n");
  while (my $line = <FILE>) {
    if ($line =~ /Main:\s+(\d+\.\d+)/) {
      $current_main = $1;
    }
    if ($line =~ /Panther:\s+(\d+\.\d+)/) {
      $current_panther = $1;
    }
  }
  close(FILE);

  return ($current_main, $current_panther);
}

##########################################
# get the current release versions from the web page

sub get_versions {

  my $page = "/tmp/interpro.page";
  my $main;
  my $panther;
  `wget --quiet -O $page ftp://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/DATA/`;

  open (PAGE, "<$page") || $log->log_and_die("Can't open the InterPro web page file $page\n");
  while (my $line = <PAGE>) {
    chomp $line;
    if ($line =~ /latest_nopthr/) {
      ($main) = $line =~ /iprscan_DATA_(\d+\.\d+)/;
    }
    if ($line =~ /latest_pthr/) {
      ($panther) = $line =~ /iprscan_PTHR_DATA_(\d+\.\d+)/;
    }
  }
  close(PAGE);

  if ($main !~ /\d+\.\d+/ || $panther !~ /\d+\.\d+/) {
    $log->log_and_die("Can't read the database versions from the FTP site ftp://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/DATA/\nMain: $main\nPanther: $panther\n");
  }

  return ($main, $panther);
}

##########################################
# get the specified database file and unpack it

sub get_file {

  my ($type, $version) = @_;
  my $filename;
  my $dir = "/data/blastdb/Worms/interpro_scan"; 
  if ($type eq 'panther') {
    $filename = "iprscan_PTHR_DATA_${version}.tar.gz";

    $log->write_to("Delete the old Panther files (this takes a while)\n");
    $wormbase->run_command("rm -rf $dir/iprscan/data/Panther", $log);

  } else {
    $filename = "iprscan_DATA_${version}.tar.gz";

    $log->write_to("Delete the old index files\n");
    $wormbase->run_command("rm -rf $dir/iprscan/data/*.inx", $log);
    $wormbase->run_command("rm -rf $dir/iprscan/data/*.bin", $log);
    $wormbase->run_command("rm -rf $dir/iprscan/data/*.psq", $log);
    $wormbase->run_command("rm -rf $dir/iprscan/data/*.phr", $log);
    $wormbase->run_command("rm -rf $dir/iprscan/data/*.pin", $log);

  }

  $log->write_to("Downloading file $filename\n");

#  `wget --quiet -O $dir/$filename ftp://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/DATA/$filename`;
  `wget -O $dir/$filename ftp://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/DATA/$filename`;

  $log->write_to("Gunzipping file $filename\n");
  `gunzip -c "$dir/$filename" | tar -xvf -`;
  unlink "$dir/$filename";
}

##########################################
# delete the interpro results from the specied species' database
# &delete_results($species_to_delete);

sub delete_results {

  my ($get_main, $get_panther, $species) = @_;


  # mysql database parameters
  my $dbhost = "ia64d";
  my $dbuser = "wormadmin";          
  my $dbname = "worm_ensembl_" . $species; # construct the database name for this species
#  my $dbpass = "XXXXXX";

  # connect to the mysql database
  $log->write_to("connect to the mysql database $dbname on $dbhost as $dbuser\n\n");
  my $dbh = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
      || $log->log_and_die("cannot connect to db, $DBI::errstr");

  
  # logical names in the Ensemble pipeline for the analyses used in interpro
  my @methods = qw(scanprosite Prints pfscan blastprodom hmmpanther Smart Tigrfam Pfam PIRSF Superfamily gene3d);

  # get the mapping of method 2 analysis id
  my %method2analysis;
  $log->write_to("get mapping of method to analysis id \n");
  my $sth = $dbh->prepare ( q{ SELECT analysis_id
				   FROM analysis
				   WHERE logic_name = ?
				 } );
  
  foreach my $method (@methods) {
    $sth->execute ($method);
    (my $anal) = $sth->fetchrow_array;
    $method2analysis{$method} = $anal; 
#    $log->write_to("$method  $anal\n");
#    print "$method  $anal\n" if ($verbose);
  }

  # if we only have a new panther database then only delete the panther results
  if ($get_panther && !$get_main) {
    @methods = qw(hmmpanther);
  }
  # if we don't have a new panther database then don't delete the panther results
  if ($get_main && !$get_panther) {
    @methods = grep !/hmmpanther/, @methods;
  }

  $sth = $dbh->prepare ( q{ DELETE FROM protein_feature
				WHERE analysis_id = ?
			      } );
  my $sth2 = $dbh->prepare ( q{ DELETE FROM input_id_analysis 
				    WHERE analysis_id = ?
				  } );

# delete the results of previous runs of the analysis pipeline for interpro
  foreach my $method (@methods) {
    foreach my $method (@methods) { # these are the analysis.analysis_id values for the Main InterPro analyses
      my $analysis_id = $method2analysis{method};
      if ($analysis_id) {
	#print "Deleting method $method, analysis_id $analysis_id\n";
	$sth->execute ($analysis_id);
	$sth2->execute ($analysis_id);
      }
    }
  }

  $sth->finish;
  $sth2->finish;
  $dbh->disconnect;
}


##########################################
# save the version number of the database we now have
#&save_versions($latest_main, $latest_panther);

sub save_versions {

  my ($latest_main, $latest_panther) = @_;

  my $file = $wormbase->compare . "/interpro_versions.dat";

  open (FILE, "> $file") || $log->log_and_die("Can't open $file\n");
  print FILE "This holds the version numbers of the current InterPro databases\n";
  print FILE "See: BLAST_scripts/make_interpro.pl\n";
  print FILE "\n";
  print FILE "Main:\t$latest_main\n";
  print FILE "Panther:\t$latest_panther\n";
  close(FILE);

}

##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

=back

This script does...blah blah blah

script_template.pl MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Keith Bradnam (krb@sanger.ac.uk)

=back

=cut
