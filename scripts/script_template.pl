#!/usr/local/bin/perl5.6.0 -w                    # perl5.6.0 and -w flag
#
# script_template.pl                             # insert name of script
# 
# by Keith Bradnam                               # original author
#
# This is a example of a good script template    # Brief description of script
#
# Last updated by: $Author: krb $                      # These lines will get filled in by cvs and helps us
# Last updated on: $Date: 2002-06-17 12:21:08 $                        # quickly see when script was last changed and by whom


use strict;                                      # always use
use lib "/wormsrv2/scripts/";                    # Try to use Wormbase.pm where it useful to do so
use Wormbase;



# Try to keep different parts of code cleanly separated using comments...

##############
# variables  #                                                                   #
##############

# Most checking scripts should produce a log file that is a) emailed to us all 
# and b) copied to /wormsrv2/logs

my $maintainers = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
our $log        = "/wormsrv2/logs/finish_build.$rundate";


# Always good to cleanly exit from your script
exit(0);


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
