#!/usr/local/bin/perl5.8.0 -w
#
# prepare_primary_databases.pl
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2004-11-24 15:06:33 $

use strict;
my $scriptdir = glob("~ar2/wormbase/rebuild");#$ENV{'CVS_DIR'};
use lib "/nfs/team71/worm/ar2/wormbase/rebuild";#$ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;


#################################################################################
# prepare primary databases                                                     #
#################################################################################
#
# Requirements : [01] - Presence of the Primary_databases_used_in_build file
#
# Checks       : [01] - Fail if the build_in_progess flag is absent.
#              : [02] - Fail if the Primary_databases_used_in_build file is absent
#
# Does         : [01] - checks the Primary_database_used_in_build data

my ($test,$debug,$database);

GetOptions ( 
	    "test"       => \$test,
	    'debug:s'    => \$debug,
	    'database:s' => \$database
	   );

$test = 1 if ( defined $ENV{'TEST_BUILD'} and ($ENV{'TEST_BUILD'} == 1));

# establish log file.
my $log = Log_files->make_build_log($debug);

# set paths to take account of whether -test is being used
my $basedir = glob("~wormpub");
$basedir     .= "/TEST_BUILD" if ($test);

# exit if the Primary_databases_used_in_build is absent
#  &usage(13) unless (-e "$logdir/Primary_databases_used_in_build");
 
local (*LAST_VER);
my ($stlace_date,$brigdb_date,$citace_date,$cshace_date) = &FTP_versions;
my ($stlace_last,$brigdb_last,$citace_last,$cshace_last) = &last_versions;
my $options = "";

# use test mode if autoace_minder -test was specified
$options .= " -test" if ($test);

# stlace
print "\nstlace : $stlace_date last_build $stlace_last";
unless ($stlace_last eq $stlace_date) {
  $options .= " -stlace $stlace_date";
  print "  => Update stlace";
}
  
# brigdb
print "\nbrigdb : $brigdb_date last_build $brigdb_last";
unless ($brigdb_last eq $brigdb_date) {
  $options .= " -brigace $brigdb_date";
  print "  => Update brigdb";
}
  
# citace
print "\ncitace : $citace_date last_build $citace_last";
unless ($citace_last eq $citace_date) {
  $options .= " -citace $citace_date";
  print "  => Update citace";
}
  
# cshace
print "\ncshace : $cshace_date last_build $cshace_last";
unless ($cshace_last eq $cshace_date) {
  $options .= " -cshace $cshace_date";
  print "  => Update cshace";
}

print "\n\nrunning unpack_db.pl $options\n";

# confirm unpack_db details and execute
unless ($options eq "") {
  print "Do you want to unpack these databases ?\n";
  my $answer=<STDIN>;
  &usage(2) if ($answer ne "y\n");

  $log->write_to(" running unpack_db.pl $options\n\n");
  &run_command("$scriptdir/unpack_db.pl $options", $log);
}

# make a unpack_db.pl log file in /logs

if ($test) {
  $log->write_to("WARNING: Can't transfer geneace from /wormsrv1.  You will have to do that by hand!\n");  
} else {
  my $camace_orig = glob("~wormpub/camace_orig");
  # transfer /wormsrv1/camace to $basedir/camace 
  $log->write_to("Transfering geneace and camace\n");
  &run_command("$scriptdir/TransferDB.pl -start $camace_orig -end $basedir/camace -database");
  # transfer /wormsrv1/geneace to $basedir/geneace 
  &run_command("$scriptdir/TransferDB.pl -start /wormsrv1/geneace -end $basedir/geneace -database");
}

  
#################################################
# Check that the database have unpack correctly #
#################################################

$log->write_to("writing Primary_databases_used_in_build\n");
# rewrite Primary_databases_used_in_build
open (LAST_VER, ">$database/Primary_databases_used_in_build");
print LAST_VER "stlace : $stlace_date\n"; 
print LAST_VER "brigdb : $brigdb_date\n"; 
print LAST_VER "citace : $citace_date\n"; 
print LAST_VER "cshace : $cshace_date\n"; 
close LAST_VER;
  
# make a unpack_db.pl log file in /logs
#  system("touch $logdir/$flag{'A4'}");



$log->mail;
exit(0);


##################
# FTP_versions   #
##################

sub FTP_versions {

  $log->write_to("\tgetting FTP versions . . \n");
  local (*STLACE_FTP,*BRIGDB_FTP,*CITACE_FTP,*CSHACE_FTP);

  my $stlace_FTP = "/nfs/privateftp/ftp-wormbase/pub/incoming/stl/stlace_*";
  my $brigdb_FTP = "/nfs/privateftp/ftp-wormbase/pub/incoming/stl/brigdb_*";
  my $citace_FTP = "/nfs/privateftp/ftp-wormbase/pub/incoming/caltech/citace_*";
  my $cshace_FTP = "/nfs/privateftp/ftp-wormbase/pub/incoming/csh/cshl_*";
  my ($stlace_date,$brigdb_date,$citace_date,$cshace_date);

  # stlace
  open (STLACE_FTP, "/bin/ls -t $stlace_FTP |")  || die "cannot open $stlace_FTP\n";
  while (<STLACE_FTP>) {
    chomp; (/\_(\d+)\-(\d+)\-(\d+)\./); $stlace_date = substr($1,-2).$2.$3; last;
  }
  close STLACE_FTP; 

  # brigdb
  open (BRIGDB_FTP, "/bin/ls -t $brigdb_FTP |") || die "cannot open $brigdb_FTP\n";
  while (<BRIGDB_FTP>) {
    chomp; (/\_(\d+)\-(\d+)\-(\d+)\./); $brigdb_date = substr($1,-2).$2.$3; last;
  }
  close BRIGDB_FTP; 
  
  # citace
  open (CITACE_FTP, "/bin/ls -t $citace_FTP |") || die "cannot open $citace_FTP\n";
  while (<CITACE_FTP>) {
    chomp; (/\_(\d+)\-(\d+)\-(\d+)\./); $citace_date = substr($1,-2).$2.$3; last;
  }
  close CITACE_FTP; 
  
  # cshace
  open (CSHACE_FTP, "/bin/ls -t $cshace_FTP |") || die "cannot open $cshace_FTP\n";
  while (<CSHACE_FTP>) {
    chomp; (/\_(\d+)\-(\d+)\-(\d+)\./); $cshace_date = substr($1,-2).$2.$3; last;
  }
  close CSHACE_FTP; 
  
  # return current dates as 6-figure string
  return($stlace_date,$brigdb_date,$citace_date,$cshace_date);
}

#####################################################################################################################

sub last_versions {
    
  $log->write_to("\tgetting last versions . . \n");
  local (*LAST_VER);
  my ($stlace_last,$brigdb_last,$citace_last,$cshace_last);

  open (LAST_VER, "<$database/Primary_databases_used_in_build") || usage("Primary_databases_file_error");
  while (<LAST_VER>) {
    $stlace_last = $1 if /^stlace \: (\d+)$/;
    $brigdb_last = $1 if /^brigdb \: (\d+)$/;
    $citace_last = $1 if /^citace \: (\d+)$/;
    $cshace_last = $1 if /^cshace \: (\d+)$/;
  }
  close LAST_VER;

  foreach my $database ( ($stlace_last,$brigdb_last,$citace_last,$cshace_last) ) {
    unless (defined $database ){
      $database = 000000;
      $log->write_to("Database not defined - using FTP version\n");
    }
  }

  # return last version dates as 6-figure string
  return($stlace_last,$brigdb_last,$citace_last,$cshace_last);

}
sub usage {
  print "usage\n";
}
__END__

_BUILD_INFO_

<p class=indent><a name="unpack">
<a href="http://intweb.sanger.ac.uk/Projects/C_elegans/PERL/autoace_minder.txt.shtml">
<font color=red>autoace_minder.pl</a></font><font color=red> -unpack</font> (-test)<br>
+ Checks the date stamps of database dumps on the FTP site and compares with those used in the last release<BR> 
+ Runs <a href="#unpack_db.pl"><font color=red>unpack_db.pl</font></a> on new databases
to copy & unpack database files to relevant place (e.g. <font color = green>/wormsrv2/stlace</font>)<BR>
+ Re-initialise the databases, and then loads new data<BR>
+ Copies camace and geneace from  <font color=green>/wormsrv1/</font> to <font color=green>/wormsrv2/</font>.<br>
+ Writes <font color=green>Primary_databases_used_in_build</font> lock file (contains dates of databases
used in build).<BR> 
+ Writes <font color=green>/wormsrv2/autoace/logs/A3:Unpack_FTP_databases</font> lock file<BR>
+ Writes <font color=green>/wormsrv2/autoace/logs/A4:Primary_databases_on_wormsrv2</font> lock file<BR>
<b>~1 hour</b>


<BR><BR><font color=blue><B>CHECK:</B></font><BR>
Copying camace and geneace is done by TransferDB, which will generate two log files <font color=green>/wormsrv2/logs/TransferDB.yymmdd.pid</font>. Make sure that the last line reads "SUCCESS ...".
</p>
