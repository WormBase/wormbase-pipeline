#!/usr/local/bin/perl5.8.0 -w
#
# make_interpro.pl
# 
# by Gary Williams
#
# This detects whether there is a new version of the Interpro databases available
# It copies over the latest file, nupacks it into place and runs the InterPro
# indexing program on it.
#
# Last updated by: $Author: mh6 $     
# Last updated on: $Date: 2007-12-17 17:18:27 $      

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

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);





##########################
# MAIN BODY OF SCRIPT
##########################

# mysql database parameters
my $dbhost = "ia64c";
my $dbuser = "wormro";          # worm read-only access
my $dbname = "worm_pep";
my $dbpass = "";

my $hmmpfam = 11;		# the pfam analysis_id number
my $panther = 25;		# the panther analysis_id number

# connect to the mysql database
$log->write_to("connect to the mysql database $dbname on $dbhost as $dbuser\n\n");
print "connect to the mysql database $dbname on $dbhost as $dbuser\n\n" if ($verbose);
my $dbh = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
    || $log->log_and_die("cannot connect to db, $DBI::errstr");


# read the current database release values
my $sth = $dbh->prepare ( q{ SELECT db_version
                               FROM analysis
                              WHERE analysis_id = ?
                           } );

$sth->execute ($hmmpfam);	# the Main InterPro database release value is stored here
(my $current_main) = $sth->fetchrow_array;
$sth->execute ($panther);	# the Panther InterPro database release value is stored here
(my $current_panther) = $sth->fetchrow_array;


# get the latest version numbers from the web site
my ($latest_main, $latest_panther) = &get_versions();

my ($get_main, $get_panther);
print "latest versions:  Main='$latest_main', Panther='$latest_panther'\n";
print "current versions: Main='$current_main', Panther='$current_panther'\n";
if ($latest_main eq $current_main) {
  print "No new release for Main InterPro database\n";
} else {
  print "Getting new release for Main InterPro database\n";
  $get_main = 1;
}
if ($latest_panther eq $current_panther) {
  print "No new release for Panther InterPro database\n";
} else {
  print "Getting new release for Panther InterPro database\n";
  $get_panther=1;
}

# get the database files, unpack them and update the record of which evrsion we have
$sth = $dbh->prepare ( q{ UPDATE analysis
			      SET db_version = ?
			      WHERE analysis_id = ?
			    } );
if ($get_main) {
  &get_file('main', $latest_main);
  $sth->execute ($latest_main, $hmmpfam); # update the db_version in the analysis table
}
if ($get_panther) {
  &get_file('panther', $latest_panther);
  $sth->execute ($latest_panther, $panther); # update the db_version in the analysis table
}

# index the databases if they have changed
if ($get_main || $get_panther) {
  $wormbase->run_script("perl /software/worm/iprscan/bin/index_data.pl -p /data/blastdb/Worms/interpro_scan/iprscan/data/ -inx -bin -v", $log);
}

# delete the results of previous runs of the analysis pipeline for interpro
$sth = $dbh->prepare ( q{ DELETE FROM protein_feature
			      WHERE analysis_id = ?
			    } );
my $sth2 = $dbh->prepare ( q{ DELETE FROM input_id_analysis 
				  WHERE analysis_id = ?
				} );

if ($get_main) {
  foreach my $analysis_id (11, 16..24) { # these are the analysis.analysis_id values for the Main InterPro analyses
    $sth->execute ($analysis_id);
    $sth2->execute ($analysis_id);
  }
}
if ($get_panther) {		# we only need to delete the previous results of the panther anaysis
  foreach my $analysis_id (25) { # this is the analysis.analysis_id value for the Panther analysis
    $sth->execute ($analysis_id);
    $sth2->execute ($analysis_id);
  }
}

# now do the same for worm_brigpep
$dbname = "worm_brigpep";

$log->write_to("connect to the mysql database $dbname on $dbhost as $dbuser\n\n");
print "connect to the mysql database $dbname on $dbhost as $dbuser\n\n" if ($verbose);
$dbh = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
    || $log->log_and_die("cannot connect to db, $DBI::errstr");

$sth = $dbh->prepare ( q{ DELETE FROM protein_feature
			      WHERE analysis_id = ?
			    } );
$sth2 = $dbh->prepare ( q{ DELETE FROM input_id_analysis 
				  WHERE analysis_id = ?
				} );

if ($get_main) {
  foreach my $analysis_id (11, 16..24) { # these are the analysis.analysis_id values for the Main InterPro analyses
    $sth->execute ($analysis_id);
    $sth2->execute ($analysis_id);
  }
}
if ($get_panther) {		# we only need to delete the previous results of the panther anaysis
  foreach my $analysis_id (25) { # this is the analysis.analysis_id value for the Panther analysis
    $sth->execute ($analysis_id);
    $sth2->execute ($analysis_id);
  }
}



$sth->finish;
$sth2->finish;
$dbh->disconnect;

$log->mail();
exit(0);






##############################################################
#
# Subroutines
#
##############################################################


##########################################
# get the current release versions from the web page

sub get_versions {

  my $page = "/tmp/interpro.page";
  my $main;
  my $panther;
  `wget --quiet -O $page ftp://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/DATA/`;

  open (PAGE, "<$page") || die "Can't open the InterPro web page file $page\n";
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
  return ($main, $panther);
}

##########################################
# get the specified database file and unpack it

sub get_file {

  my ($type, $version) = @_;
  my $filename;
  my $dir = "/data/blastdb/Worms/interpro_scan/"; 
  if ($type eq 'panther') {
    $filename = "iprscan_PTHR_DATA_${version}.tar.gz";
  } else {
    $filename = "iprscan_DATA_${version}.tar.gz";
  }

  `wget --quiet -O $dir/$filename ftp://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/DATA/$filename`;
  `gunzip -c "$dir/$filename" | tar -xvf -`;
  unlink "$dir/$filename";
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
