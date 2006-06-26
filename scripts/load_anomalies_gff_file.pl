#!/usr/local/bin/perl5.8.0 -w
#
# load_anomalies_gff_file.pl                           
# 
# by Gary Williams                        
#
# This looks for anomalous thnigs such as protein homologies not
# matching a CDS and stores the ersults in the mysql database
# 'worm_anomaly'
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2006-06-26 12:32:20 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
#use Sequence_extract;
use Coords_converter;
use DBI;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($database, $input, $gff_source, $gff_type, $id_after, $anomaly_type);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "database:s" => \$database,	    # use the specified database instead of autoace
	    "input:s"   => \$input, # the name of the input GFF file
	    "gff_source:s"   => \$gff_source, # the 'source' column of the GFF file
	    "gff_type:s"     => \$gff_type, # the 'type' column of the GFF file
	    "id_after:s"       => \$id_after, # the regexp of the data in the last column that precedes the ID name of this object, e.g. "CDS_name\s+",
	    "anomaly_type:s" => \$anomaly_type, # type of anomaly to store in database (can be any informative name)
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

#################################
# test for missing command-line arguments
#################################

if (!defined $input || $input eq "") {
  die "no input GFF file specified with -input\n";
}

if (!defined $gff_source || $gff_source eq "") {
  die "no GFF source (second column of GFF file) specified with -gff_source\n";
}

if (!defined $gff_type || $gff_type eq "") {
  die "no GFF type (third column of GFF file)  specified with -gff_source\n";
}

if (!defined $id_after || $id_after eq "") {
  print "WARNING - no tag for finding ID name in the last GFF column specified with -id_after - will use all of the last GFF column, or will construct IDs from $anomaly_type\n";
}

if (!defined $anomaly_type || $anomaly_type eq "") {
  die "no type for identifying this sort of anomaly specified with -anomaly_type\n";
}


#################################
# Set up some useful paths      #
#################################

# Set up top level base directories (these are different if in test mode)
my $basedir         = $wormbase->basedir;     # BASE DIR
my $ace_dir         = $wormbase->autoace;     # AUTOACE DATABASE DIR
my $wormpep_dir     = $wormbase->wormpep;     # CURRENT WORMPEP
my $wormrna_dir     = $wormbase->wormrna;     # CURRENT WORMRNA
my $common_data_dir = $wormbase->common_data; # AUTOACE COMMON_DATA
my $chromosomes_dir = $wormbase->chromosomes; # AUTOACE CHROMSOMES
my $reports_dir     = $wormbase->reports;     # AUTOACE REPORTS
my $gff_dir         = $wormbase->gff;         # AUTOACE GFF
my $gff_splits_dir  = $wormbase->gff_splits;  # AUTOACE GFF SPLIT
my $logs_dir        = $wormbase->logs;        # AUTOACE LOGS

# some database paths
my $geneace   = $wormbase->database('geneace');
my $camace    = $wormbase->database('camace');
my $currentdb = $wormbase->database('current');
my $stlace    = $wormbase->database('stlace');
my $citace    = $wormbase->database('citace');
my $cshace    = $wormbase->database('cshace');
my $brigace   = $wormbase->database('brigace');

# other paths
my $ftp_upload_dir  = $wormbase->ftp_upload;  # "/nfs/ftp_uploads/wormbase"
my $tace            = $wormbase->tace;        # TACE PATH
my $giface          = $wormbase->giface;      # GIFACE PATH



##########################
# MAIN BODY OF SCRIPT
##########################

# mysql database parameters
my $dbsn = "DBI:mysql:database=worm_anomaly;host=ecs1f";
my $dbuser = "wormadmin";
my $dbpass = "worms";

my $mysql = DBI -> connect($dbsn, $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to database, $DBI::errstr";

# get the last used anomaly_id key value
my $array_ref = $mysql->selectcol_arrayref("select max(anomaly_id) from anomaly;");
my $db_key_id = $array_ref->[0];
my $go_faster_by_ignoring_db_checks = 0;
if (defined $db_key_id) {
  print "db_key_id=$db_key_id\n";
} else {
  # reset the daatbase key value
  $db_key_id = 0; 
  print "db_key_id has been reset to 0\n";
  $go_faster_by_ignoring_db_checks = 1;	# don't need to check records in the database because there are none
}

$database = $wormbase->autoace if (!defined $database || $database eq "");

my $coords = Coords_converter->invoke($database, 0, $wormbase);
my %clonelab = $wormbase->FetchData("clone2centre", undef, "$database/COMMON_DATA");


my %GFF_data =   (
  file		=> $input,
  gff_source	=> $gff_source,
  gff_type	=> $gff_type,
  anomaly_type	=> $anomaly_type,
  ID_after	=> $id_after,
);

&read_GFF_file(\%GFF_data);


# disconnect from the mysql database
$mysql->disconnect || die "error disconnecting from database", $DBI::errstr;

$log->mail();

print "Finished.\n" if ($verbose);
exit(0);


##############################################################
#
# Subroutines
#
##############################################################

##########################################
# use the data in the hash-ref to read in GFF files
# various actions can be performed on the lines to filter and process them

sub read_GFF_file {
  my ($GFF_data) = @_;

  my @files = glob($GFF_data->{'file'});

  my $score = 1.0;		# default score
  my $id;
  my $id_count = 1;		# count to construct an ID if these are not held in the GFF file

  foreach my $file (@files) {
    open (GFF, "< $file") || die "Can't open $file\n";
    while (my $line = <GFF>) {
      chomp $line;
      if ($line =~ /^\s*$/) {next;}
      if ($line =~ /^#/) {next;}
	  my @f = split /\t/, $line;
	  my ($chromosome, $source, $type, $start, $end, $sense) = ($f[0], $f[1], $f[2], $f[3], $f[4], $f[6]);
	  if ($GFF_data->{'gff_source'} ne "" && $GFF_data->{'gff_source'} ne $source) {next;}
	  if ($GFF_data->{'gff_type'} ne "" && $GFF_data->{'gff_type'} ne $type) {next;}
	  if ($GFF_data->{'ID_after'} ne "") {
	    ($id) = ($f[8] =~ /$GFF_data->{'ID_after'}(\S+)/);
	  } else {
	    # there is no tag for the ID - see if there is anything in the last column we could use
	    if ($f[8] ne "") {
	      $id = $f[8];
	    } else {
	      # OK, simply construct a fake ID name
	      $id = $GFF_data->{'anomaly_type'} . "_" . $id_count++;
	    }
	  }
	  if (! defined $id) {next;}
	  $id =~ s/\"//g;	# remove quotes
	  if ($chromosome =~ /CHROMOSOME_(\S+)/) {$chromosome = $1;} # abbreviate chromosome
	  
	  #print "read GFF: $chromosome, $id, $start, $end, $sense\n";

	  # save this line directly to the database
	  &output_to_database($GFF_data->{'anomaly_type'}, $chromosome, $id, $start, $end, $sense, $score, '');

    }
    close (GFF);
  }

}

##########################################

sub output_to_database {
  my ($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, $explanation ) = @_;

  # get the clone and lab for this location
  my ($clone, $clone_start, $clone_end) = $coords->LocateSpan($chromosome, $chrom_start, $chrom_end);
  my $lab = $clonelab{$clone};          # get the lab that sequenced this clone
  if ($clone =~ /CHROMOSOME/) {
    $lab = "HX/RW";
  } elsif ($clone =~ /SUPERLINK_(\S\S)/) {
    if ($1 eq "CB") {
      $lab = "HX";
    } else {
      $lab = "RW";
    }
  }
 
  # calculate the window value as blocks of 10 kb
  my $window =  int($chrom_start/10000);

  # need to do a check for a very similar previous record that may
  # need to be overwritten, preserving the status.

  my $nearest_db_key_id = -1;
  my $nearest = 100;	# this size of distance will cause a new record to be inserted if it is not changed to <= 20

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
      #print "db_id=$result_row->[0] chrom_start=$result_row->[1] chrom_end=$result_row->[2]\n";
      #print "searching distance of start pos = $dist_start end pos = $dist_end, nearest = $nearest\n";
      if ($dist_start + $dist_end < $nearest) {
	$nearest = $dist_start + $dist_end;
	$nearest_db_key_id = $result_row->[0];
	#print "got a new best distance $nearest\n";
      }
    }
  }



  # is the distance less than 20 bases?
  # if it is not zero this is probably a move of the anomaly 
  # as a result of genome sequence changes or
  # changes in the blast database size.
  # so we should update the existing record
  if ($nearest <= 20) {
    $mysql->do(qq{ UPDATE anomaly SET   clone="$clone", clone_start=$clone_start, clone_end=$clone_end, centre="$lab", chromosome_start=$chrom_start, chromosome_end=$chrom_end, thing_score=$anomaly_score, explanation="$explanation"   WHERE anomaly_id = $nearest_db_key_id; });
    # NB we do not write the status record for this anomaly_id

  } else {

    # we want a new record inserted
    # write the data to the database
    $db_key_id++;
    $mysql->do(qq{ insert into anomaly values ($db_key_id, "$anomaly_type", "$clone", $clone_start, $clone_end, "$lab", "$chromosome", $chrom_start, $chrom_end, "$chrom_strand", "$anomaly_id", $anomaly_score, "$explanation", $window, 1, NULL); });
    #print "*** inserting new record\n";

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

=head2 NAME - load_anomalies_gff_file.pl

=head1 USAGE

=over 4

=item  load_anomalies_gff_file [-options]

=back

This script loads the contents of a GFF file into the database of anomalies to be curated.
It loads only those lines in the GFF file that match the specied GFF source and GFF type names.

Each entry in the anomalies database has an anomaly type that explains
what the anomaly is, e.g. "MANGLED_RNAI_FEATURE" Each GFF file should
have a different anomaly type. It can be up to 32 characters long - it
is truncated if it is longer.

The last line of the GFF typically holds an ID for the thing that is
being used as evidence for this anomaly, usually an accession for a
protein or EST, for Feature ID.  This is often preceded by a tag name,
e.g. 'Transcript "xyz"'. The -id_after argument is used to define the
regular expression for the tag and the spaces (if any),
e.g. 'Transcript\s+'. The quote marks are always stripped out - there
is no need to put these in the regular expression. If there is no data
in the last column then an ID will be constructed out of the anomaly
type and a count.

Each entry loaded gets the maximum anomaly score - there may be provision in the future to use the GFF score column.

Subsequent running of this script with the same file will 

script_template.pl MANDATORY arguments:

=over 4

=item -database

=back

use the specified database instead of autoace

=over 4

=item -input

=back

the name of the input GFF file

=over 4

=item -gff_source

=back

the 'source' column of the GFF file

=over 4

=item -gff_type

=back

the 'type' column of the GFF file

=over 4

=item -id_after

=back

the regexp of the data in the last column that precedes the ID name of this object, e.g. "CDS_name\s+",

=over 4

=item -anomaly_type

=back

type of anomaly to store in database (can be any informative name)




load_anomalies_gff.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Verbose/Debug mode
 
=back

=over 4

=item -test, Test mode, generate the acefile but do not upload themrun the script, but don't change anything

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back
                                                                                             

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can write to the mysql database used to hold the 'worm_anomaly' database

=back

=head1 AUTHOR

=over 4

=item Gary Williams (gw3@sanger.ac.uk)

=back

=cut
