#!/usr/local/bin/perl5.6.1 -w
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
# Last updated by: $Author: krb $
# Last updated on: $Date: 2003-06-03 16:26:46 $



$| = 1;
use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use IO::Handle;
use Getopt::Std;
use vars qw($opt_h);
use Common_data;

#################################################################################
# variables                                                                     #
#################################################################################

my $maintainers = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
our $log        = "/wormsrv2/logs/finish_build.$rundate";

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
system("TransferDB.pl -start /wormsrv2/autoace -end /wormsrv2/$WS_name -database -wspec -chromosomes -name $WS_name") 
  && die "couldn't run TransferDB for autoace\n";

# Remove redundant files from /wormsrv2/autoace/release and /wormsrv2/autoace/CHROMOSOMES
print LOG "Removing old files in /wormsrv2/autoace/release/\n";
system("rm -f $db_path/autoace/release/*") && die "Couldn't remove old release files\n";
print LOG "Removing old files in /wormsrv2/autoace/CHROMOSOMES/\n";
system("rm -f $db_path/autoace/CHROMOSOMES/*") && die "Couldn't remove old CHROMOSOME files\n";

# Remove redundant files from /wormsrv2/autoace/logs
 print LOG "Removing old files in /wormsrv2/autoace/logs\n";
system("rm -f $db_path/autoace/logs/*:*") && die "Couldn't remove old log files\n";

# Transfer autoace to ~wormpub/DATABASES/current_DB
print LOG "Transferring autoace to ~wormpub/DATABASES/current_DB\n";
system("TransferDB.pl -start /wormsrv2/autoace -end /nfs/disk100/wormpub/DATABASES/current_DB -all -name $WS_name")  && die "couldn't run TransferDB for wormpub\n";
system("/bin/gunzip /nfs/disk100/wormpub/DATABASES/current_DB/CHROMOSOMES/*.gz") && die "Couldn't gunzip CHROMOSOMES/*.gz\n";


# archive old GFF splits directory'
print LOG "Archiving GFFsplits directory using GFFsplitter.pl -a\n\n";
system("GFFsplitter.pl -a") && die "Couldn't run GFFsplitter.pl -a\n";

# run locus2seq.pl
print LOG "Running locus2seq.pl\n\n";
system("locus2seq.pl -camace ") && die "Couldn't run locus2seq.pl -a\n";


# update all Common_data files - see Commom_data.pm
system("update_Common_data.pl -update -in_build -all") && die "Couldn't run update_Common_data.pl -update -in_build -all\n";

# update "Confirmed Introns" webpage (introns to be addressed)
system("/nfs/intweb/cgi-bin/wormpub/confirmed_introns/parse_gff.pl") && warn "Couldn't run parse_gff.pl\n";

##################
# End
##################

print LOG "C'est finis.\n";
close(LOG);
&mail_maintainer("WormBase Report: finish_build.pl",$maintainers,$log);


exit(0);




#################################################################################
# set up log file                                                               #
#################################################################################

sub create_log_file{

  open (LOG,">$log") || die "Cannot open logfile $!\n";
  LOG->autoflush();
  
  print LOG "# finish_build.pl\n\n";     
  print LOG "# run details    : $rundate $runtime\n";
  print LOG "# WormBase version : $WS_name\n";
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

  # turn the old release into a tarball, move into /nfs/wormarchive and remove old directory
  print LOG "\nCreating $WS_old_path.tar.gz\n";
  system ("tar -cvf $WS_old_path.tar $WS_old_path/") && die "Couldn't create tar file\n";
  system ("gzip $WS_old_path.tar") && die "Couldn't create gzip file\n";
  print LOG "Moving archive to /nfs/wormarchive and removing $WS_old_path\n";
  system ("mv $WS_old_path.tar.gz /nfs/wormarchive/") && die "Couldn't move to /nfs/wormarchive\n";
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
