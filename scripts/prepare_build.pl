#!/usr/local/bin/perl5.6.0 -w
#
# prepare_build.pl
# 
# by Keith Bradnam aged 12 and a half
#
# Usage : prepare_build.pl [-options]
#
# This script replaces archive_dbs.pl, the script that would be run at the start of build.
# It does what that script used to do, i.e.
# 1) checks to see if there are three existing (and unpacked) WS releases in /wormsrv2
#    If there are, then it archives the oldest release away into /wormsrv2/wormbase_archive
# 2) Does a similar thing with Wormpep releases in /wormsrv2/WORMPEP
# but it also does a few more things that have to be done before the build proper can start.

$| = 1;
use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use IO::Handle;
use Getopt::Std;
use vars qw($opt_h);

#################################################################################
# variables                                                                     #
#################################################################################

my $maintainers = "dl1\@sanger.ac.uk krb\@sanger.ac.uk kj2\@sanger.ac.uk";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
our $log        = "/wormsrv2/logs/prepare_build.$rundate";
my $cvs_version = &get_cvs_version("$0");

my $db_path     = "/wormsrv2";
my $WS_name     = &get_wormbase_version_name;
my $WS_current  = &get_wormbase_version;
my $WS_new      = $WS_current + 1;
my $WS_new_name = "WS".$WS_new;
my $WS_oldest   = $WS_current - 3; # the version that *should* be the oldest in /wormsrv2
my $WS_old_name = "WS".$WS_oldest;
my $WS_old_path = "$db_path"."/$WS_old_name";

my $old_wormpep = "$db_path/WORMPEP/wormpep".($WS_current-3);

 ##############################
 # command-line options       #
 ##############################

$opt_h = "";   # Help/Usage page
getopts ('h');
&usage if ($opt_h);

#####################################################################################

&create_log_file;
&archive_old_releases;

# Transfer autoace to WSxx
print LOG "Transferring autoace into /wormsrv2/$WS_name\n";
system("TransferDB -start /wormsrv2/autoace -end /wormsrv2/$WS_name -all -name $WS_name") 
  && die "couldn't run TransferDB for autoace\n";

# Remove redundant files from /wormsrv2/autoace/release and /wormsrv2/autoace/CHROMOSOMES
print LOG "Removing old files in /wormsrv2/autoace/release/\n";
system("rm -f $db_path/autoace/release/*") && die "Couldn't remove old release files\n";
print LOG "Removing old files in /wormsrv2/autoace/CHROMOSOMES/\n";
system("rm -f $db_path/autoace/CHROMOSOMES/*") && die "Couldn't remove old CHROMOSOME files\n";

# Transfer /wormsrv1/camace to /wormsrv2/camace
print LOG "Transferring /wormsrv1/camace into /wormsrv2/camace\n";
system("TransferDB -start /wormsrv1/camace -end /wormsrv2/camace -database -wspec -name camace")
  && die "Couldn't run TransferDB for camace\n";

# update symbolic link for 'current_DB'
print LOG "Updating symbolic link to point to current_DB\n\n";
system("rm -f $db_path/current_DB") && die "Couldn't remove 'current_DB' symlink\n";
system("ln -s $db_path/$WS_name/ $db_path/current_DB") && die "Couldn't create new 'Current_DB' symlink\n";

# update database.wrm using cvs
my $cvs_file = "$db_path/autoace/wspec/database.wrm";
print LOG "Updating $cvs_file to include new WS number - using CVS\n\n";
system("cvs -d '/nfs/ensembl/cvsroot/' edit $cvs_file") && die "Couldn't 'cvs edit' $cvs_file\n";
system("sed 's/$WS_name/$WS_new_name/' < $cvs_file > ${cvs_file}.new") && die "Couldn't edit $cvs_file\n";
system("mv /wormsrv2/autoace/wspec/database.wrm.new $cvs_file") && die "Couldn't update $cvs_file\n";
system("cvs -d '/nfs/ensembl/cvsroot/' commit -m 'updating $cvs_file to $WS_new_name' $cvs_file") && die "Couldn't 'cvs commit' $cvs_file\n";


##################
# End
##################

print LOG "C'est finis.\n";
close(LOG);
&mail_maintainer("WormBase Report: prepare_build.pl",$maintainers,$log);


exit(0);




#################################################################################
# set up log file                                                               #
#################################################################################

sub create_log_file{

  open (LOG,">$log") || die "Cannot open logfile $!\n";
  LOG->autoflush();
  
  print LOG "# prepare_build.pl\n\n";     
  print LOG "# run details    : $rundate $runtime\n";
  print LOG "# WormBase version : $WS_name\n";
  print LOG "# cvs version      : $cvs_version\n";
  print LOG "\n\n";

}
#################################################################################
# remove non-essential files from old database directory                        #
#################################################################################

sub archive_old_releases{
  print LOG "\n\n";

  # remove some unnecessary files before archiving
  if (-d "$WS_old_path"){
    if (-d "$WS_old_path/database"){
      print LOG "Removing $WS_old_path/database\n";
      system("rm -rf $WS_old_path/database") && die "Couldn't remove $WS_old_path/database\n";
    }
    # remove wspec, wgf etc.
    my @files = glob("$WS_old_path/w*");
    foreach my $file (@files){
      print LOG "Removing $WS_old_path/$file\n";
      system("rm -rf $file") && die "Couldn't remove $file\n";
    }
    if (-d "$WS_old_path/pictures"){
      print LOG "Removing $WS_old_path/pictures\n";
      system("rm -rf $WS_old_path/pictures") && die "Couldn't remove $WS_old_path/pictures\n";
    }    
  }

  # turn the old release into a tarball, move into wormbase_archive and remove old directory
  print LOG "\nCreating $WS_old_path.tar.gz\n";
  system ("tar -cvf $WS_old_path.tar $WS_old_path/") && die "Couldn't create tar file\n";
  system ("gzip $WS_old_path.tar") && die "Couldn't create gzip file\n";
  print LOG "Moving archive to /wormsrv1/WORMBASE_ARCHIVE and removing $WS_old_path\n";
  system ("mv $WS_old_path.tar.gz /wormsrv1/WORMBASE_ARCHIVE/") && die "Couldn't move to /wormsrv1/WORMBASE_ARCHIVE\n";
  system ("rm -rf $WS_old_path") && die "Couldn't remove old directory\n";
  
  # archive old wormpep version
  if(-d $old_wormpep){
    print LOG "\nCreating $old_wormpep archive\n";
    system ("tar -cvf $old_wormpep.tar $old_wormpep") && die "Couldn't create wormpep tar file\n";
    system ("gzip $old_wormpep.tar") && die "Couldn't create gzip wormpep tar file\n";
    system ("rm -rf $old_wormpep") && die "Couldn't remove old directory\n";
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

=head2 NAME - prepare_build.pl

=head1 USAGE

=over 4

=item prepare_build.pl  [-options]

=back

This script replaces archive_dbs.pl, the script that would be run at the 
start of the build.

 It does what that script used to do, i.e.

 1) checks to see if there are three existing (and unpacked) WS releases 
 in /wormsrv2. If there are, then it archives the oldest release away into 
 /wormsrv2/wormbase_archive
 2) Does a similar thing with Wormpep releases in /wormsrv2/WORMPEP
 but it also does a few more things that have to be done before the build 
 proper can start.

prepare_build.pl MANDATORY arguments:

=over 4

=item none

=back

prepare_build.pl  OPTIONAL arguments:

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
