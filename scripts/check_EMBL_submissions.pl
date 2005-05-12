#!/usr/local/bin/perl5.8.0 -w
#
# check_EMBL_submissions.pl
# dl
#
# Checks the return mail from the EBI against the list of submitted
# clones. Entries which have failed to load or return are highlighted
# and changes in sequence version are notified.

# Last updated on: $Date: 2005-05-12 11:09:46 $
# Last updated by: $Author: dl1 $

# touch logfile for run details
$0 =~ m/\/*([^\/]+)$/; system("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

use strict;
use Getopt::Std;
use IO::Handle;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use vars qw/ $opt_d $opt_h $opt_f/;

 ########################################
 # Script variables (run)               #
 ########################################

my $maintainer = "All";
my $rundate = &rundate;
my $runtime = &runtime;
my $log="/wormsrv2/logs/check_EMBL_submissions.$rundate.$$";

 ########################################
 # Paths/variables etc                  #
 ########################################

my $dbdir = glob "/wormsrv1/camace";
my $tace = &tace;
my $submitted_file = "/nfs/disk100/wormpub/analysis/TO_SUBMIT/submitted_to_EMBL";

my (%clone2id,%id2sv,%embl_acc,%embl_status,%embl_sv);

 ########################################
 # command-line options                 #
 ########################################

getopts ('dhf:');

my $debug = $opt_d;
my $embl_file = $opt_f;

 ########################################
 # simple consistency checks            #
 ########################################

&usage(0) if ($opt_h);
&usage(1) unless (-e $embl_file);
&usage(2) unless (-e $submitted_file);

 ########################################
 # Open logfile                         #
 ########################################

system ("/bin/touch $log");
open (LOG,">>$log");
LOG->autoflush;

print LOG "# check_EMBL_submissions\n\n";     
print LOG "# run details    : $rundate $runtime\n";
print LOG "\n";

 ########################################
 # query autoace for EMBL ID/AC lines   #
 ########################################

$ENV{'ACEDB'} = $dbdir;

my $command = "Table-maker -p $dbdir/wquery/SCRIPT:check_EMBL_submissions.def\nquit\n";

my ($sv,$id,$name);

open (TACE, "echo '$command' | $tace | ");
while (<TACE>) {
    s/acedb\> //g;
    print "$_\n" if ($debug);
    chomp;
    next if ($_ eq "");
    next if (/\/\//);

    s/\"//g;
    if (/(\S+)\s+(\S+)\s+(\S+)/) {
	($clone2id{$1} = $3) if ($2 eq "NDB_ID"); # ACEDB_clone => NDB_ID  
	if ($2 eq "NDB_SV") { # NDB_ID => NDB_SV
	    $id   = $1;
	    $sv   = $3;
	    $name = $sv =~ (/\S+\.(\d+)/);
	    $id2sv{$id} = $name;
	  } 

    }
}
close TACE;

print LOG "Populated hashes with data for " . scalar (keys %clone2id) . " clones\n\n";

 ########################################
 # returned entries from EMBL           #
 ########################################

open (EMBL, "<$embl_file") || die "Cannot open EMBL returns\n";
while (<EMBL>) {
    next unless (/^CE/);
    if (/NOT LOADED/) {
	(/^(\S+)/);
	$embl_status{$1} = "NOT_LOADED";
	next;
    }
    (/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)$/); 
    $embl_acc{$1}    = $2;
    $embl_status{$1} = $4;
    $embl_sv{$1}     = $5;
    print "processed $1\n" if ($debug);
}
close EMBL;

print LOG "Populated hashes with data for " . scalar (keys %embl_status) . " entries\n\n\n";

 #############################################################
 # submitted clones from ~wormpub/analysis/TO_SUBMIT         #
 #############################################################

 
open (SUBMITTED, "<$submitted_file") || die "Cannot open submit_TO_SUBMIT log file\n";;
while (<SUBMITTED>) {
    ($name) = (/^(\S+)/);
    if (!defined $clone2id{$name}) {print LOG "eek .. no ID for clone $name\n";}
    $id = $clone2id{$name};

    print LOG "# $name   \tSubmitted_to_EMBL\t";

    if (!defined $embl_status{$id}) {
	print LOG "not returned\n";
	next;
    }
    elsif ($embl_status{$id} eq "Finished") {
	print LOG "loaded  \t";
    }
    elsif ($embl_status{$id} eq "NOT_LOADED") {
	print LOG "failed to load\n";
	next;
    }
    
    if ($embl_sv{$id} ne $id2sv{$name}) {
	print LOG "Update sequence version\n";
    }
    else {
	print LOG "Sequence version unchanged\n";
    }
}
close SUBMITTED;
close LOG;

###############################
# Mail log to curator         #
###############################

&mail_maintainer("check_EMBL_submissions",$maintainer,$log);

###################################################
# hasta luego                                     #
###################################################

exit(0);

sub usage {
    my $error = shift;
    if ($error == 0) {
        exec ('perldoc',$0);
        exit (0);
    }
    # Error  1 - EMBL return filename is not valid
    elsif ($error == 1) {
	print "The EMBL return summary file does not exist\n\n";
        exit(0);
    }
    # Error  2 - submitted_to_EMBL file is not present
    elsif ($error == 2) {
	print "The submitted_to_EMBL summary file does not exist\n\n";
        exit(0);
    }
}

__END__

=pod

=head2 NAME - check_EMBL_submissions.pl

=head1 USAGE:

=over 4

=item check_EMBL_submissions.pl

=back

check_EMBL_submissions.pl correlates the returned summary e-mail from the EBI
with the list of clone submitted by submit_TO_SUBMIT. The output will highlight
which entries are unaccounted for and whether the sequence version field has
changed.

check_EMBL_submissions.pl mandatory arguments:

=over 4

=item -f <filename>, EMBL summary e-mail (exported from Pine)

=back

check_EMBL_submissions.pl optional arguments:

=over 4

=item -h, Help page

=item -d, Verbose/Debug mode


=back

=head1 RUN REQUIREMENTS:

=back

The submitted_to_EMBL file produced by the submit_TO_SUBMIT script is needed
as reference for the entries submitted to the EBI.

=head1 RUN OUTPUT:

=back

The log file at /wormsrv2/logs/check_EMBL_submissions.$$ contains the output
from the script. An example of this output is given below.

C01G10        Submitted_to_EMBL       loaded          Sequence version unchanged
C01G12        Submitted_to_EMBL       loaded          Update sequence version
C01G6         Submitted_to_EMBL       not returned
C01H6         Submitted_to_EMBL       loaded          Sequence version unchanged
C02B4         Submitted_to_EMBL       not returned

=head1 EXAMPLES:

=over 4

check_EMBL_submissions.pl -f /nfs/griffin2/dl1/EMBL_020123

=head1 AUTHOR: 

Daniel Lawson

Email dl1@sanger.ac.uk

=cut


