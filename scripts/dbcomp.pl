#!/usr/local/bin/perl5.8.0 -w
#
# dbcomp.pl
#
# Counts the number of objects in an ACEDB database for each Class stated in the config file
# Compares this number to those from a second database.
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2003-06-10 16:45:10 $      


use strict;
use lib '/wormsrv2/scripts';
use Wormbase;
use IO::Handle;
use Getopt::Long;
$|=1;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $database, $database2);
our ($log,$errfile,$outfile); 
our ($db_1, $db_2, $dbname_1, $dbname_2);
our $WS_current  = &get_wormbase_version;
our $WS_previous = $WS_current - 1;
my $maintainers = "All";
my $exec        = &tace;

GetOptions ("help"          => \$help,
            "debug=s"       => \$debug,
	    "database=s"    => \$database,
	    "database2=s"   => \$database2);


############################################
# Check command-line options
############################################

# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}



#################################################################
# Compare autoace to previous build unless -database specified
#################################################################


$dbname_1    = "WS${WS_previous}";
$db_1        = "/wormsrv2/$dbname_1"; 

$dbname_2    = "WS${WS_current}";
$db_2        = "/wormsrv2/autoace";


# First alternative datatbase specified?
if($database){
  $dbname_1  = "$database";
  $db_1      = "$database"; 
}

# Second alternative datatbase specified?
if($database2){
  $dbname_2  = "$database2";
  $db_2      = "$database2"; 
}




&create_log_files;


#########################################################################
# Read list of classes
# The complete set of classes to dump is between the DATA and END tokens
#########################################################################

my @TotalClasses;

READARRAY:   while (<DATA>) {
  chomp $_;
  last READARRAY if $_ =~ /END/;
  push (@TotalClasses,$_);
}




#######################
# Main loop
#######################

my $counter = 1; # for indexing each class to be counted

foreach my $query (@TotalClasses) {
  next if ($query eq "");
  LOG->autoflush();

  print LOG  " Counting '$query'\n";
  printf OUT " | %18s |", $query;

  ##################################################
  # Get class counts from both databases           #
  ################################################## 

  my ($class_count_1,$class_count_2) = &count_class($query,$counter);

  print LOG " Counting $dbname_1 : $class_count_1\n";  
  print LOG " Counting $dbname_2 : $class_count_2\n\n";
  
  printf OUT "%8s",$class_count_1;
  print OUT " | ";  
  printf OUT "%8s",$class_count_2;


  ##################################################
  # Calculate difference between databases         #
  ##################################################
  
  my ($added, $removed) = &diff($counter);
  my $diff = $added - $removed;

#  printf OUT "| %6s |\n",$diff;

  printf OUT "| %7s ",$added;
  printf OUT "| %7s ",$removed;
  printf OUT "| %7s |\n",$diff;

  $counter++;
}


print OUT  " +--------------------+---------+---------+---------+---------+---------+\n";





# create symbolic links to current.out and current.dbcomp
print LOG "\nCreating 'current.out' and 'current.dbcomp'\n";
system("rm -f /wormsrv2/autoace/COMPARE/current.out") && die "Couldn't remove 'current.out' symlink\n";
system("rm -f /wormsrv2/autoace/COMPARE/current.dbcomp") && die "Couldn't remove 'current.dbcomp' symlink\n";
system("ln -s $errfile /wormsrv2/autoace/COMPARE/current.out") && die "Couldn't create new symlink\n";
system("ln -s $outfile /wormsrv2/autoace/COMPARE/current.dbcomp") && die "Couldn't create new symlink\n";


close (OUT);
close (ERR);
close (LOG);

# Email log file
&mail_maintainer("WormBase Report: dbcomp.pl",$maintainers,$log);

#write to the release letter - subroutine in Wormbase.pm
&release_databases;

exit (0);




##############################################################
#
# Subroutines
#
##############################################################

sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate     = `date +%y%m%d`; chomp $rundate;
  $log        = "/wormsrv2/logs/$script_name.$rundate.$$";

  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name\n";
  print LOG "started at ",`date`,"\n";
  print LOG "=============================================\n";
  print LOG "\n";

  print LOG "\n";
  print LOG "Previous db      : $dbname_1 '$db_1'\n";
  print LOG "Current db       : $dbname_2 '$db_2'\n";
  print LOG "\n\n";

  # open two main output files to store results
  $errfile = "/wormsrv2/autoace/COMPARE/WS${WS_previous}-WS${WS_current}.out";
  $outfile = "/wormsrv2/autoace/COMPARE/WS${WS_previous}-WS${WS_current}.dbcomp";

  open (OUT, ">$outfile") || die "Couldn't write to out file\n";
  open (ERR, ">$errfile") || die "Couldn't write to err file\n";
  
  print OUT  " +----------------------------------------+\n";
  print OUT  " | Class              |   ACEDB database  |\n";
  print OUT  " |                    +---------+---------+---------+---------+---------+\n";
  printf OUT " |                    | %7s | %7s |    +    |    -    |   Net   |\n", $dbname_1,$dbname_2;
  print OUT  " +--------------------+---------+---------+---------+---------+---------+\n";




  LOG->autoflush();
}

##########################################
sub count_class{
  my $query  = shift;
  my $counter = shift;
  my $out;
  my $class_count1;
  my $class_count2;

  # Formulate query
  my $command=<<EOF;
query find $query 
list -a 
quit
EOF

  ####################################
  # Count objects in first database
  ####################################

  # open temp output file
  $out = "/tmp/dbcomp_A_${counter}";
  open (COUNT, ">$out") || die "Couldn't write to tmp file: $out\n";

  # open tace connection and count how many objects in that class
  open (TACE, "echo '$command' | $exec $db_1 | ");
  while (<TACE>) {
    ($class_count1 = $1) if (/^\/\/ (\d+) Active Objects/);
    (print COUNT "$_") if (/\:/); # Add list of object to temp file
  }
  close (TACE);
  close (COUNT);


  ####################################
  # Count objects in second database
  ####################################

  $out = "/tmp/dbcomp_B_${counter}";
  open (COUNT, ">$out") || die "Couldn't write to tmp file: $out\n";

  # open tace connection and count how many objects in that class
  open (TACE, "echo '$command' | $exec $db_2 | ");
  while (<TACE>) {
    ($class_count2 = $1) if (/^\/\/ (\d+) Active Objects/);
    (print COUNT "$_") if (/\:/); # Add list of object to temp file
  }
  close (TACE);
  close (COUNT);


  return ($class_count1, $class_count2);


}


##################################################

sub diff {

  # look at the differences between the two sets of tmp files (one pair of files
  # for each class being queried)
  my $counter = shift;

  my $added   = 0; 
  my $removed = 0;

  system ("cat /tmp/dbcomp_A_${counter} | sort > /tmp/look-1");
  system ("cat /tmp/dbcomp_B_${counter} | sort > /tmp/look-2");
  open (COMM, "comm -3 /tmp/look-1 /tmp/look-2 |");
  while (<COMM>) {
    if (/^(\S+.+)/){
      print ERR " <- $dbname_1 $1\n";
      $removed++;
    }
    if (/^\s+(\S+.+)/){
      print ERR " -> $dbname_2 $1\n";
      $added++;
    }
  }   
  close (COMM);
  system ("rm -f /tmp/look-1");
  system ("rm -f /tmp/look-2");

  system ("rm -f /tmp/dbcomp_A_${counter}");
  system ("rm -f /tmp/dbcomp_B_${counter}");
 
  # Add break symbol to output file to separate classes
  print ERR "\/\/ -----------------------------\n";
  return($added, $removed);
}

###############################################



sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}





__DATA__
2_point_data
Accession_number
Allele
Author
Cell
Cell_group
Class
Clone
Comment
Contig
Database
Display
DNA
Expr_pattern
Expr_profile
Feature
Feature_data
Gene_Class
Gene_name
GO_term
Homol_data
Journal
Keyword
Laboratory
Life_stage
Lineage
Locus
Map
Method
Microarray_aff
Microarray_result
Motif
Movie
Multi_pt_data
Oligo
Operon
Paper
PCR_product
Peptide
Person
Person_name
Phenotype
Picture
Pos_neg_data
Protein
Rearrangement
Reference
Repeat_Info
RNAi
Sequence
Session
SK_map
Species
Strain
Table
Tag
Transgene
Transcript
Url
__END__


=pod

=head2   NAME - dbcomp.pl

=head1 USAGE

=over 4

=item dbcomp.pl [-options]

=back

This script will (by default) perform a class by class comparison of autoace and 
the previous WormBase release database.  It will calculate what database this
will be but you can force the choice of the second database by using the -database
option.  Classes to be compared are specified within the script and should be amended
as necessary when classes are added/removed to the database.


=over 4

=item MANDATORY arguments: none

=back

=over 4

=item OPTIONAL arguments: -debug, -help, -database


-debug and -help are standard Wormbase script options.

-database allows you to specify the path of a database which will
then be compared to /wormsrv2/autoace.

-database2 allows you to specify a second database to compare to that specified
by -database (rather than compare with previous build)


=back

=head1 AUTHOR - Dan Lawson (but completely rewritten by Keith Bradnam)


Email krb@sanger.ac.uk



=cut

