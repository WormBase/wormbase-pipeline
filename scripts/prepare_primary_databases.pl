#!/usr/local/bin/perl5.8.0 -w
#
# prepare_primary_databases.pl
#
# Last edited by: $Author: gw3 $
# Last edited on: $Date: 2006-01-05 09:50:38 $

use strict;
my $scriptdir = $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;


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

my ($test,$debug,$database, $store, $wormbase);

GetOptions ( 
	    "test"       => \$test,
	    'debug:s'    => \$debug,
	    'database:s' => \$database,
	    'store:s'    => \$store
	    
	   );

if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

my %databases;
$databases{'stlace'}->{'search'} = 'stl/stlace*';
$databases{'brigdb'}->{'search'} = 'stl/brigdb*';
$databases{'brigdb'}->{'option'} = 'brigace'; # for this database the option passed to unpack is not db name <sigh>
$databases{'citace'}->{'search'} = 'caltech/citace*';
$databases{'cshace'}->{'search'} = 'csh/cshace*';

&FTP_versions;
&last_versions;
my $options = "";

# use test mode if autoace_minder -test was specified
$options .= " -test" if ($test);

foreach my $primary (keys %databases){
  if( $databases{$primary}->{'ftp_date'} ) {
    unless ( $databases{$primary}->{'last_date'} == $databases{$primary}->{'ftp_date'} ) {
      $options .= " -".($databases{"$primary"}->{'option'} or $primary )." ".$databases{$primary}->{'ftp_date'};
      print " => Update $primary";
    }
  }else {
    $log->write_to("no version of $primary on FTP site");
    
  }
}

print "\n\nrunning unpack_db.pl $options\n";

# confirm unpack_db details and execute
unless ($options eq "") {
  print "Do you want to unpack these databases ?\n";
  my $answer=<STDIN>;
  &usage(2) if ($answer ne "y\n");
  
  $log->write_to(" running unpack_db.pl $options\n\n");
  $wormbase->run_script("unpack_db.pl $options", $log);
}

# make a unpack_db.pl log file in /logs


# transfer camace and geneace to correct PRIMARIES dir
$log->write_to("Transfering geneace and camace\n");
$wormbase->run_script("TransferDB.pl -start ".$wormbase->database('camace'). " -end ".$wormbase->primary("camace") ." -database -wspec");
$wormbase->run_script("TransferDB.pl -start ".$wormbase->database('geneace')." -end ".$wormbase->primary("geneace")." -database -wspec");
#system("cp -R misc_static $autoace/acefiles/primary  #check whats happened here - looks like partial edit

#################################################
# Check that the database have unpack correctly #
#################################################

$log->write_to("writing Primary_databases_used_in_build\n");
# rewrite Primary_databases_used_in_build
my $new_primary = $wormbase->basedir."Primary_databases_used_in_build";
open (LAST_VER, ">$new_primary");
foreach my $primary ( keys %databases){
  print LAST_VER "$primary : ".$databases{$primary}->{'ftp_date'}."\n";
}
close LAST_VER;

$log->mail;
exit(0);


##################
# FTP_versions   #
##################

sub FTP_versions {

  $log->write_to("\tgetting FTP versions . . \n");

  my $ftp = $wormbase->ftp_upload;
  foreach my $db (keys %databases){
    open(my $dir, "/bin/ls -t $ftp/$databases{$db}->{'search'} |") or $log->log_and_die("Cant open $ftp/$databases{$db}->{'search'}"."\n");
    while (<$dir>) {
      unless (/\_\d+\-\d+\-\d+\./) {next;}
      chomp; 
      (/\_(\d+)\-(\d+)\-(\d+)\./); 
      $databases{$db}->{'ftp_date'} = substr($1,-2).$2.$3; 
      last;
    }
    close $dir;
  }
}

#####################################################################################################################

sub last_versions {
    
  $log->write_to("\tgetting last versions . . \n");

  my $file = $wormbase->basedir."/Primary_databases_used_in_build";
  open (LAST_VER, "<$file") or $log->log_and_die("cant find $file\n");
  while (<LAST_VER>) {
    if (/^(\w+) \: (\d+)$/) {
      $databases{"$1"}->{'last_date'} = $2;
    }
  }
  close LAST_VER;

  foreach my $db ( keys %databases ) {
    unless (defined $databases{$db}->{'last_date'} ) {
      $databases{$db}->{'last_date'} = '000000';
      $log->write_to("Database not defined - using FTP version\n");
    }
  }
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
+ Copies camace and geneace from  <font color=green>/nfs/disk100/wormpub/DATABASES/</font> to <font color=green>/wormsrv2/</font>.<br>
+ Writes <font color=green>Primary_databases_used_in_build</font> lock file (contains dates of databases
used in build).<BR> 
+ Writes <font color=green>/wormsrv2/autoace/logs/A3:Unpack_FTP_databases</font> lock file<BR>
+ Writes <font color=green>/wormsrv2/autoace/logs/A4:Primary_databases_on_wormsrv2</font> lock file<BR>
<b>~1 hour</b>


<BR><BR><font color=blue><B>CHECK:</B></font><BR>
Copying camace and geneace is done by TransferDB, which will generate two log files <font color=green>~wormpub/logs/TransferDB.yymmdd.pid</font>. Make sure that the last line reads "SUCCESS ...".
</p>
