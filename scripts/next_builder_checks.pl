#!/usr/local/bin/perl5.8.0 -w
#
# next_builder_checks.pl                           
# 
# by Keith Bradnam                         
#
# A simple script to send a check list to the person who will be performing the next
# build to check the current build
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2003-11-28 10:28:44 $      

use strict;
use lib "/wormsrv2/scripts/";   
use Wormbase;
use Getopt::Long;

my $WS_current    = &get_wormbase_version;
my $log         = "/tmp/next_builder_check.log";

my ($help, $debug, $user,$maintainers);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "user=s"     => \$user);


# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

# Use -user setting to send email
if($user){
  ($maintainers = $user . '\@sanger.ac.uk');
}


#######################
#
# Main part of script
#
#######################




open (LOG, ">$log");
print LOG "\nGreetings $user.  You are meant to be building the next release of WormBase, so\n";
print LOG "you should take some time to check the current release.\n";
print LOG "====================================================================================\n\n";

print LOG "THINGS TO CHECK:\n\n";

print LOG "1) The following 12 clones are representative of the whole genome in that\n";
print LOG "they include one Sanger and one St. Louis clone for each chromosome.  Check\n";
print LOG "each clone to ensure that it contains BLAT data (EST and mRNA), BLAST data,\n";
print LOG "waba data, gene models, UTRs etc.\n\n";
print LOG "i)    C25A1\n";
print LOG "ii)   F56A3\n";
print LOG "iii)  C04H5\n";
print LOG "iv)   B0432\n";
print LOG "v)    C07A9\n";
print LOG "vi)   F30H5\n";
print LOG "vii)  C10C6\n";
print LOG "viii) B0545\n";
print LOG "ix)   C12D8\n";
print LOG "x)    K04F1\n";
print LOG "xi)   C02C6\n";
print LOG "xii)  AH9\n\n";

print LOG "2) If the genome sequence has changed, you should inspect the clones containing\n";
print LOG "those changes to see if there are any strange errors (e.g. duplicate sets of data\n";
print LOG "which are slightly out of sync.)\n\n";

print LOG "3) Check /wormsrv2/autoace/CHROMOSOMES/composition.all - are there any non-ATCGN\n";
print LOG "characters\n\n";


print LOG "4a) Check that the latest WormPep proteins have proper protein and motif homologies\n";
print LOG "This has been a problem in some builds where all new WormPep proteins don't get any\n";
print LOG "BLAST analysis.  Pick a few random Wormpep proteins and especially check that all of\n";
print LOG "the various blastp homologies are there (human, fly, worm, yeast etc.) and try to\n";
print LOG "check at least one protein from the /wormsrv2/WORMPEP/wormpepXXX/new_entries.WSXXX file\n\n";

print LOG "4b) Now that we have a curated set of brigpep, should do this periodically for\n";
print LOG "C. briggase protein objects too...these now have their own set of blastp hits\n\n";


print LOG "5) Open a keyset of all genome sequences and then query them for:\n";
print LOG " 'Transcript_child AND NEXT AND NOT NEXT'\n";
print LOG "I.e. tRNAs not attached properly to parent sequence.  This has happened before and you\n";
print LOG "should notify the responsible group (HX or RW) to fix them for next build\n\n";

print LOG "6) Check all Protein objects have a Species tag set\n\n";

print LOG "\nThat's all...for now!  If you are satisfied the build is ok, please inform the person\n";
print LOG "building the database. Please continue to add to this list as appropriate\n";

print LOG "\n========================================================================================\n\n";
close(LOG);

&mail_maintainer("Please check the ongoing build of WS${WS_current}",$maintainers,$log);


exit(0);


##############################################################
#
# Subroutines
#
##############################################################


sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##################################################################




__END__

=pod

=head1 NAME - next_builder_checks.pl

=head1 USAGE

=over 4

=item next_builder_checks.pl --user <user>

=back

This script simply sends a list of check items to the next person doing the build.
They should 'sign off' on each build and hand back to the main person when they
are happy all is ok.

=item MANDATORY arguments: -user <valid unix username>

Needed to send email

=back

=over 4

=item OPTIONAL arguments: -help, -debug <user>


 
=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk


=cut
