#!/usr/local/bin/perl5.6.1 -w
#
# backup_geneace.pl
# v 0.1
# 
# by dl
#
# Usage : backup_geneace
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2003-04-13 16:30:46 $

# touch logfile for run details
$0 =~ m/\/*([^\/]+)$/; system("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

#################################################################################
# variables                                                                     #
#################################################################################

$|=1;
use Cwd;
use IO::Handle;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Getopt::Long;
use strict;

#############################
# Command line options      #
#############################

my ($debug,$help);
GetOptions (
           "debug=s"   => \$debug,
           "help"      => \$help
           );

# help page
&usage("Help") if ($help);

# no debug name
&usage("Debug") if ((defined $debug) && ($debug eq ""));


 ##############################
 # Script variables (run)     #
 ##############################

my $maintainers = "All";
my $rundate    = `date +%y%m%d`;   chomp $rundate;
my $runtime    = `date +%H:%M:%S`; chomp $runtime;
my $logfile = "/wormsrv2/logs/$1.$rundate.$$";

 ##############################
 # paths for I/O files        #
 ##############################

my $dbpath     = "/wormsrv1/geneace";
my $backuppath = "/nfs/disk100/wormpub/acedb/ace4/geneace";
my $tace       = &tace;
my $giface     = &giface;

 ##############################
 # open logfile               #
 ##############################

open (LOGFILE,">$logfile");
LOGFILE->autoflush;

print LOGFILE "# backup_geneace\n#\n";
print LOGFILE "# run details         : $rundate $runtime\n";
print LOGFILE "# Source directory    : $dbpath\n";
print LOGFILE "# Target directory    : $backuppath\n#\n";
print LOGFILE "\n";

 ##########################################
 # copy the tar.gz file from the ftp site #
 ##########################################

# move to database directory
chdir $dbpath;
my $dir = cwd();
print "Move to directory: '$dir'\n";

# tace command to dump the database
my $command = "dump -T\nquit\n";
&DbWrite($command,$tace,$dbpath,"Dump");
print "Dump files to be compressed:\n";

my $file;
opendir (LIST, $dbpath) or die ("Could not opendir $dbpath");
while (defined($file=readdir(LIST))) {
    next unless ($file =~ /dump\S+.+\.ace/);
    print "compressing $file ... \n";
    print LOGFILE "compressing $file ... \n";
    system ("/bin/gzip $file");
    system ("/usr/bin/mv $file.gz $backuppath/");
    print " and moving to /nfs/disk100/wormpub ..\n";
    print LOGFILE " and moving to /nfs/disk100/wormpub ..\n";
}
close LIST;

print LOGFILE "\narchive procedure complete\n";
print LOGFILE "\n\n";
close LOGFILE;

###############################
# Mail log to curator(s)      #
###############################

&mail_maintainer("Backup_geneace",$maintainers,$logfile);

###############################
# hasta luego                 #
###############################

exit(0);

###################################################
# Subroutine for writing to a given database      #   
###################################################

sub DbWrite {
  my ($command,$exec,$dir,$name)=@_;
  open (WRITEDB,"| $exec $dir >> $logfile") or do {print LOGFILE "$name DbWrite failed\n";close LOGFILE; die();};
  print WRITEDB $command;
  close WRITEDB;
}

############################################################
# Subroutine for displaying messages of debug options      #   
############################################################

sub usage {
     my $error = shift;

     if ($error eq "Help") {
         # Normal help menu
         system ('perldoc',$0);
         exit (0);
     }
     elsif ($error eq "Debug") {
         # No debug bod named
         print "You haven't supplied your name\nI won't run in debug mode
         until i know who you are\n";
        exit (0);
    }
}

__END__

=pod

=head2   NAME - backup_geneace

=head1 USAGE

=over 4

=item backup_geneace

=back

backup_geneace makes an acedump of the geneace database from /wormsrv1 to
/nfs/disk100/wormpub.

=back

=head1 AUTHOR

=over 4

=item Dan Lawson (dl1@sanger.ac.uk)

=back

=cut

