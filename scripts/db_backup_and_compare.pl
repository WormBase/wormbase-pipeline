#!/usr/local/bin/perl5.6.1 -w
#
# db_backup_and_compare.pl
#
# backup database and compare to last backed up database to look for lost data
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2003-04-10 11:00:32 $      

use strict;
use lib '/wormsrv2/scripts';
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Carp;

######################################
# variables and command-line options # 
######################################

our ($help,$debug,$log,@backups);
my $db;
our $backup_dir = "/nfs/disk100/wormpub/DATABASES/BACKUPS";
our $maintainers = "All";
our $date = `date +%y%m%d`; chomp $date;
my $exec        = &tace;

GetOptions ("help"    => \$help,
	    "db=s"    => \$db,
            "debug=s" => \$debug);



# Check for mandatory/correct command line options
&check_options($db);


# Create log files
&create_log_files;


# Look to see what backups are already there, make new backup if appropriate
&find_and_make_backups($db);


# Compare new backup to previous backup, eeport any data that has disappeared
&compare_backups($db);


# tidy up, email log
close (LOG);
&mail_maintainer("$db backup and comparison",$maintainers,$log);
exit (0);

# C'est la fin




##############################################################
#
#              T H E   S U B R O U T I N E S
#
##############################################################


############################################
# Check command-line options
############################################

sub check_options{
  my $db = shift;

  # Display help if required
  &usage("Help") if ($help);
  &usage("db") if (!$db);
  &usage("db") unless (($db eq "camace") || ($db eq "geneace"));

  # Use debug mode?
  if($debug){
    print "DEBUG = \"$debug\"\n\n";
    ($maintainers = $debug . '\@sanger.ac.uk');
  }
}


#####################################################
# Find what backups have been made for the database
#####################################################

sub find_and_make_backups{
  my $db = shift;
  my $backup_dbs = "/${backup_dir}/${db}_backup.*";
  open (BACKUP_DB, "/bin/ls -d -1 -t $backup_dbs |")  || croak "cannot open $backup_dbs\n";
  while (<BACKUP_DB>) { 
    chomp;
    (/$db\_backup\.(\d+)$/); 
    # All database dates get added to @backups array, first element is last backup
    push(@backups,$1);
  }
  close(BACKUP_DB);


  # quit if backup has already been made for today
  if($date eq $backups[0]){
    print LOG "Last backup is dated today which means no backup is needed\n";
    print LOG "This script is ending early.  Goodbye\n\n";
    &mail_maintainer("$db backup and comparison",$maintainers,$log);    
    exit(0);
  }
  # otherwise make a new backup
  else{
    # keep TransferDB logs in backup directory
    chdir("$backup_dir") || print LOG "Couldn't cd to $backup_dir\n";
    print LOG "Making new backup - ${db}_backup\.${date}\n";
    my $return = system("TransferDB.pl -start /wormsrv1/$db -end ${backup_dir}/${db}_backup\.${date} -database -wspec -name ${db}\.${date}"); 
    if($return != 0){
      print LOG "ERROR: Couldn't run TransferDB.pl correctly.  Check log\n";
      close(LOG);
      &mail_maintainer("$db backup and comparison",$maintainers,$log);    
      exit;
    }
    # Now need to remove the oldest database (assuming that there are now five backups).
    if (scalar(@backups) == "4"){
      print LOG "Removing oldest backup - ${db}_backup\.${backups[3]}\n\n";
      system("rm -rf ${backup_dir}/${db}_backup\.${backups[3]}");
    }
  }
}


################################################################
# Compare last pair of backups
###############################################################

sub compare_backups{
  my $db = shift;
  
  # $db1 = most recent database, $db2 = next recent
  my $db1 = "${backup_dir}/${db}_backup\.${date}";
  my $db2 = "${backup_dir}/${db}_backup\.${backups[0]}";
  my $db_name1 = "${db}\.${date}";
  my $db_name2 = "${db}\.${backups[1]}";


  print LOG "First database:  $db_name1\n";
  print LOG "Second database: $db_name2\n\n";
  print LOG "Objects lost from the first database (in comparison to second database):\n";
  print LOG "------------------------------------------------------------------------\n\n";

  my $counter = 1; # for indexing each class to be counted

  # Read list of all classes to compare (listed at bottom of script)
 READARRAY:   while (<DATA>) {
   chomp $_;
   last READARRAY if $_ =~ /END/;
   my $query = $_; 
   next if ($query eq "");
   
   # Get class counts from both databases        
   &count_class($query,$counter,$db1,$db2); 
   $counter++;
 }
  
}

#############################################

sub create_log_files{
  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate     = `date +%y%m%d`; chomp $rundate;
  $log        = "/wormsrv2/logs/$script_name.$rundate.$$";

  open (LOG, ">$log") or croak "cant open $log";
  print LOG "=============================================\n";
  print LOG "$script_name\n";
  print LOG "started at ",`date`;
  print LOG "=============================================\n\n";

  # open two main output files to store results
  LOG->autoflush();
}

##########################################
sub count_class{
  my $query  = shift;
  my $counter = shift;
  my $db1 = shift;
  my $db2 = shift;
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
  open (COUNT, ">$out") || croak "Couldn't write to tmp file: $out\n";

  # open tace connection and count how many objects in that class
  open (TACE, "echo '$command' | $exec $db1 | ");
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
  open (COUNT, ">$out") || croak "Couldn't write to tmp file: $out\n";

  # open tace connection and count how many objects in that class
  open (TACE, "echo '$command' | $exec $db2 | ");
  while (<TACE>) {
    ($class_count2 = $1) if (/^\/\/ (\d+) Active Objects/);
    (print COUNT "$_") if (/\:/); # Add list of object to temp file
  }
  close (TACE);
  close (COUNT);

  #########################################################
  # Calculate difference between databases    
  #########################################################
  system ("cat /tmp/dbcomp_A_${counter} | sort > /tmp/look-1"); 
  system ("cat /tmp/dbcomp_B_${counter} | sort > /tmp/look-2");
  open (COMM, "comm -3 /tmp/look-1 /tmp/look-2 |");
  while (<COMM>) {
    next if (/\/\//);
    # Only write to log file where db1 has lost data in respect to db2 (older database)
    print LOG "$1\n" if (/^\s+(\S+.+)/);

  }
  
  close (COMM);
  system ("rm -f /tmp/look-1")  && carp "Couldn't remove /tmp/look-1 file\n";
  system ("rm -f /tmp/look-2")  && carp "Couldn't remove /tmp/look-2 file\n";
  system ("rm -f /tmp/dbcomp_*") && carp "Couldn't remove /tmp/dbcomp* files\n";

}



###############################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }

  elsif ($error eq "db") {
    print "\nYou must specify -db camace OR -db geneace\n";
    exit (0);
  }
}

########################################################################



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

=head2   NAME - db_backup_and_compare.pl

=head1 USAGE

=over 4

=item db_backup_and_compare.pl [-options]

=back

This script will attempt to backup a database (either camace or geneace) and the compare
that backup with the previous backup (whatever is the most recent before that) and mail
a list of data that has been lost from the newest database (in comparison to the next oldest).

This script will keep cycling through 4 backup databases, such that after a new backup is 
made, the oldest backup will be removed, leaving four backup databases.


=over 4

=item MANDATORY arguments: none

=back

=over 4

=item OPTIONAL arguments: -debug, -help, -database


-debug and -help are standard Wormbase script options.

-database allows you to specify the path of a database which will
then be compared to /wormsrv2/autoace.


=back

=head1 AUTHOR - Dan Lawson (but completely rewritten by Keith Bradnam)


Email krb@sanger.ac.uk



=cut

