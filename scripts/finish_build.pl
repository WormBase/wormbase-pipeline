#!/usr/local/bin/perl5.8.0 -w
#
# finish_build.pl
# 
# by Keith Bradnam aged 12 and a half
#
# Usage : finish_build.pl [-options]
#
# 1) checks to see if there are three existing (and unpacked) WS releases in /wormsrv2
#    If there are, then it archives the oldest release away into /nfs/wormarchive
# 2) Does a similar thing with Wormpep releases in /wormsrv2/WORMPEP
# 3) Archives old GFF_SPLITS directory
# 4) Makes current_DB (copy of latest release) in ~wormpub/DATABASES
#
# Last updated by: $Author: pad $
# Last updated on: $Date: 2005-01-20 16:19:55 $


use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Cwd;
use File::Glob ':glob';
use File::Path;


##############################
# command-line options       #
##############################

my $help;          # Help/Usage page
my $test;          # use test environment in ~wormpub/TEST_BUILD/
my $verbose;       # Verbose mode

GetOptions ("help"         => \$help,
            "test"         => \$test);


&usage if ($help);



#################################################################################
# variables                                                                     #
#################################################################################

my $maintainers = "All";
my $log;

my $basedir     = "/wormsrv2";
$basedir        = glob("~wormpub")."/TEST_BUILD" if ($test); 

# need to change this if in test mode
my $WS_name;
my $WS_current;

if ($test) {
    $WS_current = "666";
    $WS_name    = "WS666";
}
else {
    $WS_current = &get_wormbase_version;
    $WS_name    = &get_wormbase_version_name;
}

my $WS_new      = $WS_current + 1;
my $WS_new_name = "WS".$WS_new;
my $WS_oldest   = $WS_current - 3; # the version that *should* be the oldest in /wormsrv2
my $WS_old_name = "WS".$WS_oldest;
my $WS_old_path = "$basedir"."/$WS_old_name";
my $old_wormpep = "$basedir/WORMPEP/wormpep".($WS_current-3);

#####################################################################################

&create_log_files;

&archive_old_releases;

# update all Common_data files - see Commom_data.pm
system("update_Common_data.pl -build -all") && die "Couldn't run update_Common_data.pl -update -in_build -all\n";


# Transfer autoace to WSxx
print LOG "Transferring autoace into /wormsrv2/$WS_name\n";
system("TransferDB.pl -start /wormsrv2/autoace -end /wormsrv2/$WS_name -database -release -wspec -chromosomes -acefiles -name $WS_name") 
  && die "couldn't run TransferDB for autoace\n";

# Transfer autoace to ~wormpub/DATABASES/current_DB - first remove existing files
print LOG "Removing ~wormpub/DATABASES/current_DB/database/\n";
&delete_files_from("/nfs/disk100/wormpub/DATABASES/current_DB/database","*","+");

print LOG "Removing ~wormpub/DATABASES/current_DB/database/CHROMOSOMES/\n";
&delete_files_from("/nfs/disk100/wormpub/DATABASES/current_DB/CHROMOSOMES","*","+");



print LOG "Running TransferDB.pl to copy autoace to ~wormpub/DATABASES/current_DB\n";
system("TransferDB.pl -start $basedir/autoace -end /nfs/disk100/wormpub/DATABASES/current_DB -database -chromosomes -wspec -name $WS_name")  && die "couldn't run TransferDB for wormpub\n";
print LOG "Unzipping any gzipped chromosome files\n";
system("/bin/gunzip /nfs/disk100/wormpub/DATABASES/current_DB/CHROMOSOMES/*.gz") && die "Couldn't gunzip CHROMOSOMES/*.gz\n";


# Remove redundant files and directories in /wormsrv2/autoace/
print LOG "Removing old files in /wormsrv2/autoace/release/\n";
&delete_files_from("/wormsrv2/autoace/release","*","-");

print LOG "Removing old files in /wormsrv2/autoace/acefiles/\n";
&delete_files_from("/wormsrv2/autoace/acefiles","*","-");

print LOG "Removing old CHROMOSOME files in /wormsrv2/autoace/CHROMOSOMES/\n";
&delete_files_from("/wormsrv2/autoace/CHROMOSOMES","*","-");

print LOG "Removing *.wrm files in /wormsrv2/autoace/database/\n";
&delete_files_from("/wormsrv2/autoace/database",".\.wrm","-") or print LOG "ERROR: Problems removing files from autoace/database\n";

print LOG "Removing files in /wormsrv2/autoace/database/new/\n";
&delete_files_from("/wormsrv2/autoace/database/new","*","+");
print LOG "Removing files in /wormsrv2/autoace/database/touched/\n";
&delete_files_from("/wormsrv2/autoace/database/touched","*","+");

print LOG "Removing old files in /wormsrv2/autoace/logs\n";
&delete_files_from("$basedir/autoace/logs",":","-");
# Exception needed because we need to keep one file (Primary_databases) and this log file uses different name, so
# *:* won't remove it.
unlink("$basedir/autoace/logs/UTR_gff_dump");             



# archive old GFF splits directory'
print LOG "Archiving GFFsplits directory using GFFsplitter.pl -a\n\n";
system("GFFsplitter.pl -archive") && die "Couldn't run GFFsplitter.pl -a\n";


##################
# End
##################

print LOG "hasta luego\n";
close(LOG);
&mail_maintainer("WormBase Report: finish_build.pl",$maintainers,$log);


exit(0);




#################################################################################
# set up log file                                                               #
#################################################################################
sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch $basedir/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate = &rundate;
  $log        = "$basedir/logs/$script_name.$WS_name.$rundate.$$";
  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name\n";
  print LOG "started at ",&rundate,"\n";
  print LOG "=============================================\n";
  print LOG "\n";

}


#################################################################################
# remove non-essential files from old database directory                        #
#################################################################################



sub archive_old_releases{
  print LOG "\n\n";

  my @list;
  my $file;

  # remove some unnecessary files before archiving
  if (-d "$WS_old_path") {

      if (-d "$WS_old_path/database") {
	  print LOG "Removing $WS_old_path/database\n";
	  
	  # remove contents of the database folder
	  &delete_files_from("$WS_old_path/database/new","*","+");
	  &delete_files_from("$WS_old_path/database/oldlogs","*","+");
	  &delete_files_from("$WS_old_path/database/readlocks","*","+");
	  &delete_files_from("$WS_old_path/database/touched","*","+");
	  &delete_files_from("$WS_old_path/database","*","+");
      }
      
      # remove wspec, wgf etc.
      print LOG "Removing $WS_old_path/wspec\n";
      &delete_files_from("$WS_old_path/wspec","*","+");
            
      if (-d "$WS_old_path/pictures"){
	  print LOG "Removing $WS_old_path/pictures\n";
	  &delete_files_from("$WS_old_path/pictures","*","+");
      }    
  }
  

  # turn the old release into a tarball, move into /nfs/wormarchive and remove old directory
  print LOG "\nCreating $WS_old_path.tar.gz\n";
  system ("tar -P /$basedir/ -cvf $WS_old_path.tar $WS_old_path/") && die "Couldn't create tar file\n";
  system ("gzip $WS_old_path.tar") && die "Couldn't create gzip file\n";
  print LOG "Moving archive to /nfs/wormarchive and removing $WS_old_path\n";
  system ("mv $WS_old_path.tar.gz /nfs/wormarchive/") && die "Couldn't move to /nfs/wormarchive\n";
  &delete_files_from("$WS_old_path","*","+");
 
  # archive old wormpep version
  if (-d $old_wormpep) {
      print LOG "\nCreating $old_wormpep archive\n";
      system ("tar -P /$basedir/ -cvf $old_wormpep.tar $old_wormpep") && die "Couldn't create wormpep tar file\n";
      system ("gzip $old_wormpep.tar") && die "Couldn't create gzip wormpep tar file\n";
      &delete_files_from("$old_wormpep","*","+");
  }
  
}

#################################################################################
# Prints help and disappears                                                    #
#################################################################################

sub usage {
    exec ('perldoc',$0);
}


__END__

=pod

=head2 NAME - finish_build.pl

=head1 USAGE

=over 4

=item finish_build.pl  [-options]

=back

This script:

 1) checks to see if there are three existing (and unpacked) WS releases 
 in /wormsrv2. If there are, then it archives the oldest release away into 
 /wormsrv2/wormbase_archive
 2) Does a similar thing with Wormpep releases in /wormsrv2/WORMPEP
 but 
 3) Runs GFFsplitter.pl -a to archive away the last GFF_SPLITS directory
 4) Copies autoace into a separate WSxx directory
 5) updates the /wormsrv2/current_DB symlink to point to the directory created
    in step 4.

finish_build.pl MANDATORY arguments:

=over 4

=item none

=back

finish_build.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Keith Bradnam (krb@sanger.ac.uk)

=back

=cut
