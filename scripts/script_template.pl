#!/usr/local/bin/perl5.6.1 -w                  
#
# script_template.pl                           
# 
# by Keith Bradnam                         
#
# This is a example of a good script template   
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2003-01-24 15:47:53 $      

use strict;                                      
use lib "/wormsrv2/scripts/";                    
use Wormbase;
use Getopt::Long;


######################################
# variables and command-line options # 
######################################

my ($help, $debug);
my $maintainers = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
our $log;

GetOptions ("help"      => \$help,
            "debug=s"   => \$debug);


# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

&create_log_files;



##########################
# MAIN BODY OF SCRIPT
##########################






# Close log files and exit

&mail_maintainer("script template",$maintainers,$log);
close(LOG);
exit(0);






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





# Add perl documentation in POD format
# This should expand on your brief description above and add details of any options
# that can be used with the program.  Such documentation can be viewed using the perldoc
# command.


__END__

=pod

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

=back

This script:

 1) checks to see if there are three existing (and unpacked) WS releases 
 in /wormsrv2. If there are, then it archives the oldest release away into 
 /wormsrv2/wormbase_archive
 2) Does a similar thing with Wormpep releases in /wormsrv2/WORMPEP
 but 
 3) Runs GFFsplitter -a to archive away the last GFF_SPLITS directory
 4) Copies autoace into a separate WSxx directory
 5) updates the /wormsrv2/current_DB symlink to point to the directory created
    in step 4.

script_template.pl MANDATORY arguments:

=over 4

=item none

=back

script_template.pl  OPTIONAL arguments:

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
