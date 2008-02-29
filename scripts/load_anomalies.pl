#!/software/bin/perl -w
#
# find_anomalies.pl                           
# 
# by Gary Williams                        
#
# This read curation anomalies from a data file prepared by Sanger and
# stores the results in the mysql database 'worm_anomaly'
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2008-02-29 11:35:37 $      

use strict;                                      
use Getopt::Long;
use Carp;
use DBI;

######################################
# variables and command-line options # 
######################################

my ($help, $test, $verbose, $store, $wormbase);
my ($datafile, $species);

GetOptions ("help"       => \$help,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "datafile:s" => \$datafile,
	    "species:s"  => \$species,
	    );

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);
}

##########################
# MAIN BODY OF SCRIPT
##########################

# mysql database parameters - change this to your local naming scheme.
# This is also the 'base' name of any non-elegans anomaly curation
# databases you may create in the future, e.g. "worm_anomaly_briggsae"

my $sqldb = "worm_anomaly";

# if the species is anything other than 'elegans' (the default) then
# we use the database that has a name formed from the 'base' database
# name and that species' name

if ($species && $species ne 'elegans') {$sqldb = $sqldb . "_" . lc $species;}

##########################
# CHANGE YOUR DETAILS HERE
##########################
my $dbsn = "DBI:mysql:database=$sqldb;host=ia64d";
my $dbuser = "wormadmin";
my $dbpass = "worms";

my $mysql = DBI -> connect($dbsn, $dbuser, $dbpass, {RaiseError => 1})
      || die "ERROR cannot connect to database $sqldb\n, $DBI::errstr";

# get the last used anomaly_id key value
# This can be done easily in mysql 5.x, but older versions of mysql
# may require this explicit method

my $array_ref = $mysql->selectcol_arrayref("select max(anomaly_id) from anomaly;");
my $db_key_id = $array_ref->[0];
my $go_faster_by_ignoring_db_checks = 0;
if (defined $db_key_id) {
  #print "db_key_id=$db_key_id\n";
} else {
  # reset the database key value
  $db_key_id = 0; 
  #print "db_key_id has been reset to 0\n";
  $go_faster_by_ignoring_db_checks = 1;	# don't need to check records in the database because there are none
}


# now read the data file
&read_data($datafile, $species);


# disconnect from the mysql database
$mysql->disconnect || die "ERROR disconnecting from database\n", $DBI::errstr;

print "Finished.\n" if ($verbose);
exit(0);

##########################################
# 
# Subroutines
# 
##########################################

# read in the data and parse the action to be performed in each line

sub read_data {
  my ($datafile, $species) = @_;

  open(DAT, "< $datafile") || die "Can't open $datafile: $!\n";
  while (my $line = <DAT>) {
    chomp $line;
    next if ($line =~ /^#/);         # comment line
    next if ($line =~ /^\s*$/);      # blank line

    my @line = split /\t/, $line;    # the data is TAB delimited
    my $key = shift @line;           # get the action key for this line

    if ($key eq 'SPECIES') {

      # check that the species that this data was prepared from
      # matches the species of the database that we are writing into

      if (!defined $species || $species eq '') {$species = 'elegans';} # set the default name explicitly for the error message to show it
      next if ($species eq $line[0]);
      die "ERROR: This data file was prepared from the '$line[0]' genome\nbut the '$species' mysql database has been specified\n";

    } elsif ($key eq 'DELETE') {

      # now delete things that have not been updated in this run that you
      # would expect to have been updated like protein-homology-based
      # anomalies that might have gone away.  This also deletes anomalies
      # that we are no longer putting into the database and which can be
      # removed.

      &delete_anomalies($line[0]);
    

    } elsif ($key eq 'INSERT') {

      &output_to_database(@line);

    } else {

      die "ERROR unknown type of data line in file $datafile:\n$line\n";

    }

  }
  close(DAT);
}

##########################################
# checks if there is sense = + or - and if not, outputs two anomaly records, one for each sense.
# else it just outputs one record for the specified sense.

sub output_to_database {
  my ($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, $explanation, $clone, $clone_start, $clone_end, $lab) = @_;

  # calculate the window value as blocks of 10 kb
  my $window =  int($chrom_start/10000);

  # if there is no strand information, put one anomaly in for each strand
  if ($chrom_strand ne '+' && $chrom_strand ne '-') {
    &tidy_up_senseless_records($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, '.', $window, $anomaly_score, $explanation, $clone, $clone_start, $clone_end, $lab);
    &put_anomaly_record_in_database($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, '+', $window, $anomaly_score, $explanation, $clone, $clone_start, $clone_end, $lab);
    &put_anomaly_record_in_database($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, '-', $window, $anomaly_score, $explanation, $clone, $clone_start, $clone_end, $lab);

  } else {
    &put_anomaly_record_in_database($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, $chrom_strand, $window, $anomaly_score, $explanation, $clone, $clone_start, $clone_end, $lab);

  }
}

##########################################
# output the record to the database.  Checks if there is a record
# there already and doesn't change the 'ignore' flag if there is
# already a record in existence.

sub put_anomaly_record_in_database {

  my ($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, $chrom_strand, $window, $anomaly_score, $explanation, $clone, $clone_start, $clone_end, $lab) = @_;


  # ignore this if the score is less than 0.01
  if ($anomaly_score < 0.01) {return;}

  ########################################
  # write the anomaly data to the database
  ########################################

  # need to do a check for a very similar previous record that may
  # need to be overwritten, preserving the status.

  my $nearest_db_key_id = -1;
  my $nearest = 21;	# this default size of distance will cause a new record to be inserted if it is not changed in the test below to be <= 20

  if (! $go_faster_by_ignoring_db_checks) {
    # see if there is something there already
    my $db_query = $mysql->prepare ( qq{ SELECT anomaly_id, chromosome_start, chromosome_end FROM anomaly WHERE type = "$anomaly_type" AND chromosome = "$chromosome" AND sense = "$chrom_strand" AND thing_id = "$anomaly_id" AND window = $window});

    $db_query->execute();
    my $ref_results = $db_query->fetchall_arrayref;

    # find the nearest one to the current data
    #print "\tstart search for $anomaly_id\n";
    foreach my $result_row (@$ref_results) {
      my $dist_start = abs($chrom_start - $result_row->[1]);
      my $dist_end = abs($chrom_end - $result_row->[2]);
      if ($dist_start + $dist_end < $nearest) {
	$nearest = $dist_start + $dist_end;
	$nearest_db_key_id = $result_row->[0];
      }
    }
  }

  # is the distance in $nearest less than 20 bases, rather than the default size of 21?
  # if it is not zero this is probably a move of the anomaly 
  # as a result of genome sequence changes or
  # changes in the blast database size.
  # so we should update the existing record
  if ($test) {
    print "In test mode - so not updating the mysql database\n";
  } else {
    if ($nearest <= 20) {
      $mysql->do(qq{ UPDATE anomaly SET clone="$clone", clone_start=$clone_start, clone_end=$clone_end, centre="$lab", chromosome_start=$chrom_start, chromosome_end=$chrom_end, thing_score=$anomaly_score, explanation="$explanation"   WHERE anomaly_id = $nearest_db_key_id; });
      # NB we do not write the status record for this anomaly_id
      print "*** updating existing record for $anomaly_id in $clone\n" if ($verbose);

    } else {

      # we want a new record inserted
      # write the data to the database
      $db_key_id++;
      $mysql->do(qq{ insert into anomaly values ($db_key_id, "$anomaly_type", "$clone", $clone_start, $clone_end, "$lab", "$chromosome", $chrom_start, $chrom_end, "$chrom_strand", "$anomaly_id", $anomaly_score, "$explanation", $window, 1, NULL); });
      print "*** inserting new record for $anomaly_id in $clone\n" if ($verbose);
    }
  }

}
##########################################
# there have been a lot of records put in the database that have a sense
# as '.'  we want to change these to be two records with a sense '+
# and a sense '-' (this routine will not be required after the first
# time it has been run and cleaned up the database)

sub tidy_up_senseless_records {

  my ($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, $chrom_strand, $window, $anomaly_score, $explanation, $clone, $clone_start, $clone_end, $lab) = @_;

  my $old_db_key_id;
  my $active;

  my $db_query = $mysql->prepare ( qq{ SELECT anomaly_id, active FROM anomaly WHERE type = "$anomaly_type" AND chromosome = "$chromosome" AND sense = "." AND chromosome_start = $chrom_start AND chromosome_end = $chrom_end AND thing_id = "$anomaly_id" });

  $db_query->execute();
  my $ref_results = $db_query->fetchall_arrayref;

  # find the nearest one to the current data
  foreach my $result_row (@$ref_results) {
    $old_db_key_id = $result_row->[0];
    $active = $result_row->[1];
    
    # we want two new records inserted
    $db_key_id++;
    if ($test) {
      print "In test mode - so not updating the mysql database\n";
    } else {

      $mysql->do(qq{ INSERT INTO anomaly VALUES ($db_key_id, "$anomaly_type", "$clone", $clone_start, $clone_end, "$lab", "$chromosome", $chrom_start, $chrom_end, "+", "$anomaly_id", $anomaly_score, "$explanation", $window, $active, NULL); });
      $db_key_id++;
      $mysql->do(qq{ INSERT INTO anomaly VALUES ($db_key_id, "$anomaly_type", "$clone", $clone_start, $clone_end, "$lab", "$chromosome", $chrom_start, $chrom_end, "-", "$anomaly_id", $anomaly_score, "$explanation", $window, $active, NULL); });
      
    
      # and we want to delete the old record with no sense
      $mysql->do(qq{ DELETE FROM anomaly WHERE anomaly_id = $old_db_key_id; });
    }
  }
}

##########################################
##########################################
# now delete things that have not been updated in this run that you
# would expect to have been updated like protein-homology-based
# anomalies that might have gone away.  This also deletes anomalies
# that we are no longer putting into the database and which can now be
# removed.
#
# This means that an anomaly based on something like a similarity to a
# protein where the protein is no longer existing will be removed from
# the database

sub delete_anomalies{

  my ($type) = @_;

  # Delete anything that is still active (not being ignored)
  # and that is of the required type

  if ($test) {
    print "In test mode - so not updating the mysql database\n";
  } else {
    $mysql->do(qq{ DELETE FROM anomaly WHERE type = "$type" AND active = 1 });
    print "*** deleting old $type anomalies\n" if ($verbose);
  }
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
# This should expand on your brief description above and add details of any options
# that can be used with the program.  Such documentation can be viewed using the perldoc
# command.


__END__

=pod

=head2 NAME - load_anomalies.pl

=head1 USAGE

=over 4

=item load_anomalies.pl  [-options]

=back

This script populates the worm_anomaly mysql database with data describing some types of anomalies.

It is intended to be run by St. Louis using database prepared by Sanger during the Build.

It reads the database file 'anomalies.dat' produced by the script find_anomalies.pl

These anomalies can be inspected using the script:

history_maker.pl -anomalies -chromosome X

script_template.pl MANDATORY arguments:

=over 4

=item -datafile name of datafile to read in

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -species

This specifies a database other than the default elegans database in the event that we start making curation anomaly data for other species.

The default is 'elegans'.

The default name of the elegans mysql database is 'worm_anomaly'.

The name of the mysql database used by other species is formed by adding an underscore and their name to this, e.g. 'worm_anomaly_briggsae'.

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Verbose/Debug mode
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything in the mysql database.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back
                                                                                             

=head1 REQUIREMENTS

=over 4

=item There must be a mysql database server running.

=back

=head1 AUTHOR

=over 4

=item Gary Williams (gw3@sanger.ac.uk)

=back

=cut
