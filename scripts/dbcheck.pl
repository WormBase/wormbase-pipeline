#!/usr/local/bin/perl5.6.1 -w
#
# dbcheck.pl
# v 0.2
#
# Cronjob integrity check controls for generic ACEDB database.
#
# Usage: dbcheck [-options]
#
# Last updated on: $Date: 2002-12-18 15:46:23 $
# Last updated by: $Author: ck1 $

# touch logfile for run details
$0 =~ m/\/*([^\/]+)$/; system("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");


#################################################################################
# variables                                                                     #
#################################################################################

$|=1;
BEGIN {
  unshift (@INC,"/nfs/disk92/PerlSource/Bioperl/Releases/bioperl-0.05");
}

use strict;
use Bio::Seq;
use IO::Handle;
use Getopt::Std;
use Ace;
require "/nfs/disk100/wormpub/analysis/scripts/babel.pl";

 ##############################
 # Script variables (run)     #
 ##############################

my $maintainer = "dl1\@sanger.ac.uk";
my $rundate = `date +%y%m%d`; chomp $rundate;
my $runtime = `date +%H:%M:%S`; chomp $runtime;
my $version = get_script_version("dbcheck.pl");

 ##############################
 # command-line options       #
 ##############################

my $opt_d="";   # Verbose debug mode
my $opt_h="";   # Help/Usage page


my $opt_s="";   # check Status tags




getopts ('hd');
&usage(0) if ($opt_h);
my $debug = $opt_d;

 ##############################
 # Paths etc                  #
 ##############################

my $tace = glob("~wormpub/ACEDB/bin.ALPHA_4/tace");    # tace executable path
my $dbpath = "/wormsrv2/autoace";                                  # Database path

 ########################################
 # Open logfile                         #
 ########################################

my $log="/wormsrv2/logs/dbcheck.$rundate";

open (LOG,">$log");
LOG->autoflush();

print LOG "# dbcheck\n";     
print LOG "# version        : $version\n";
print LOG "# run details    : $rundate $runtime\n";
print LOG "\n";

 ########################################
 # Connect with acedb server            #
 ########################################

my $db = Ace->connect(-path=>$dbpath);

print LOG "DBCheck run STARTED at $runtime\n\n";

 #####################################################################
 # Do the following for every sequence in current.cosmid             #
 #####################################################################

my($i, $obj);

$i = $db->fetch_many(-query=> 'find Genome_Sequence "*"');  
while ($obj = $i->next) {
    my $clone = $obj;
    
    print "[$clone]    \t" if ($debug);

 ######################################################################
 # Retrieve the sequence and date FROM ACEDB                          #
 ######################################################################

    $obj = $db->fetch(Sequence=>$clone);
    if (!defined ($obj)) {
	print LOG "Could not fetch sequence $clone\n";
	next;
    }
    
  #####################################################################
  # Push the sequence as string in $seq                               #
  #####################################################################
    
    my $sequence=$obj->asDNA();
    if (!$sequence) {
	print LOG "NO DNA in ACEDB $clone\n" ;
	next;
    }
    $sequence =~ s/\>\w+//mg;
    $sequence =~ tr/a-z/A-Z/;
    $sequence =~ s/\W+//mg;
    
	
 ######################################################################
 # Iterative checks for each clone                                    #
 ######################################################################

    print " [";

 ######################################################################
 # Check for N's in FINISHED sequences                                #
 ######################################################################
	
    print "N's" if ($debug);
    &checkchars;

 ######################################################################
 # Check correctness of gene structure                                #
 ######################################################################

    print " | CDS_coords" if ($debug);
    &checkgenes;

 ######################################################################
 # last check complete, tidy up                                       #
 ######################################################################

    print "]\n" if ($debug);

 ######################################################################
 # Get rid of this sequence object                                    #
 ######################################################################

    $obj->DESTROY();
}


##################
# LINK objects 
###################

$i = $db->fetch_many(-query=> 'find Sequence "*LINK*"');  
while ($obj = $i->next) {
    my $link = $obj;
    
   
    print "[$link]    \t" if ($debug);
    print " | CDS_coords" if ($debug);
    &checkgenes;
    print "]\n" if ($debug);
}

$runtime = `date +%H:%M:%S`; chomp $runtime;
print LOG "\nDBCheck run ENDED at $runtime\n\n";

close LOG;

 ##############################
 # mail $maintainer report    #
 ##############################

open (OUTLOG,"|/usr/bin/mailx -s dbcheck_report $maintainer ");
open (READLOG, "<$log");
while (<READLOG>) {
    print OUTLOG "$_";
}
close READLOG;
close OUTLOG;

 ##############################
 # Write log to wormpub intweb
 ##############################

&writehtml;


 ##############################
 # hasta luego                #
 ##############################
$db->close;
exit(0);

########################################################################################
####################################   Subroutines   ###################################
########################################################################################
   
#######################################################################
# Odd chars and N's in  finished sequences                            #
#######################################################################

sub checkchars {

    my ($seq_ace, $clone);

    if ($seq_ace =~ /[^ACGTUMRWSYKVHDBXN]/img) {
	print LOG "ACEDBSEQ for $clone contains bad characters\n";
	$seq_ace =~s/[^ACGTUMRWSYKVHDBXN]//img;
    }
    if ($seq_ace =~ /N/g) { 
	print LOG "ACEDBSEQ FINISHED SEQUENCE for $clone contains N \n";
    }
}


#######################################################################
# Gene length as declared is subsequence and in exons list            #
#######################################################################

sub checkgenes {

    foreach my $child ($obj->Subsequence) {
	undef my @num;
	undef my ($method);
	undef my ($source);
	
	my ($seq, $start, $end, $link) = $child->row();
	
	unless ($seq =~ /\./) {next;}

	if ($link =~ /^SUPERLINK/) {print LOG "SUPERLINK $link contains CDS gene models\n";}
	
	my $diff = $end - $start;
	if ($diff < 0) {
	    $diff = $start - $end;
	}
	
	my $subseq = $db->fetch(Sequence => "$child");
	if (!defined ($subseq)) {
	    print LOG "Cannot fetch subsequence $child\n";
	    next;
	}

	# Method
	$method = $subseq->Method(1);
	if ((!defined ($method))) {
	    print LOG "The subsequence $child has no method\n";
	}
	
	#Species
	my $species = $subseq->Species(1);
	if ((!defined ($species))) {
	    print LOG "The subsequence $child [$method] has no Species tag\n";
	}

	@num = $subseq->Source_Exons(2);
	if (!@num) {
	    print LOG "The subsequence $child [$method]  has no Source_Exons\n";
	}
	my $clone;
	my $index = $#num;
	my $length = ($num[$index])-1;
	if ($diff != $length) {
	    print LOG "The subsequence $child [$method]  belonging to $clone has diff=$diff and length=$length\n";
	}

	
	if ($child =~ /\S+\.t\d+/) {
	    my $tag = $subseq->Transcript;
	    if ((!defined ($tag))) {
		print LOG "The subsequence $child  [$method] has no Transcript tag\n";
	    }
	}
	else {
	    my $tag = $subseq->CDS;
	    if ( (!defined ($tag)) && ($method ne "Pseudogene") ) {
		print LOG "The subsequence $child  [$method] has no CDS tag\n";
	    }
	}
	
	$source = $subseq->Source(1);
	if ((!defined ($source))) {
	    print LOG "The subsequence $child  [$method] has no Source\n";
	    next;
	}
	
	$subseq->DESTROY();
	$child->DESTROY();
	$diff="";
	$length="";
    }
}

#######################################################################
# Write HTML page with maintenance job results
#######################################################################

sub writehtml {

    my (@finished);
    my (@annotated);
    my (@date);
    my (@sequence);
    my (@subsequence);
    my (@link);
    
    my $logdir="/nfs/disk100/wormpub/LocalWWW";

my $HTML_START=<<START;
<HTML>
<HEAD>
<TITLE>Automated db Maintenance log</TITLE>
</HEAD>
<BODY BGCOLOR="WHITE">
START

my $HTML_END=<<END;
</BODY>
</HTML>
END

    open (OUTHTML,">$logdir/dbchecklog.html");
    print OUTHTML $HTML_START;
    print OUTHTML "<TABLE BORDER=1 WIDTH=100%>\n";

    open (READLOG, "<$log");
    while (<READLOG>) {
	push (@finished,$_)    if (/^NOT_FINISHED/);
	push (@annotated,$_)   if (/^FINISHED_BUT_NOT_ANNOTATED/);
	push (@date,$_)        if (/^DATE/);
	push (@sequence,$_)    if (/^SEQUENCE/);
	push (@subsequence,$_) if (/^The subsequence/);
    }
    close READLOG;

    # Not Finished
    print OUTHTML "<TABLE BORDER=\"0\" WIDTH=\"100%\">\n";
    print OUTHTML "<TH BGCOLOR=\"darkblue\"><TD><FONT COLOR=\"white\">Not Finished</FONT></TD></TH>\n";
    foreach (@finished) {
	print OUTHTML "<TR><TD><FONT SIZE=\"-1\">$_</FONT></TD></TR>\n";
    }
    print OUTHTML "</TABLE></P>\n";


    # Not Annotated
    print OUTHTML "<TABLE BORDER=\"0\" WIDTH=\"100%\">\n";
    print OUTHTML "<TH BGCOLOR=\"darkblue\"><TD><FONT COLOR=\"white\">Not Annotated</FONT></TD></TH>\n";
    foreach (@annotated) {
	print OUTHTML "<TR><TD><FONT SIZE=\"-1\">$_</FONT></TD></TR>\n";
    }
    print OUTHTML "</TABLE></P>\n";

    # Date Mismatch (not uploaded into camace)
    print OUTHTML "<TABLE BORDER=\"0\" WIDTH=\"100%\">\n";
    print OUTHTML "<TH BGCOLOR=\"darkblue\"><TD><FONT COLOR=\"white\">Date mismatches</FONT></TD></TH>\n";
    foreach (@date) {
	print OUTHTML "<TR><TD><FONT SIZE=\"-1\">$_</FONT></TD></TR>\n";
    }
    print OUTHTML "</TABLE></P>\n";

    # Sequence Mismatch (problem)
    print OUTHTML "<TABLE BORDER=\"0\" WIDTH=\"100%\">\n";
    print OUTHTML "<TH BGCOLOR=\"darkblue\"><TD><FONT COLOR=\"white\">Genomic sequences</FONT></TD></TH>\n";
    foreach (@sequence) {
	print OUTHTML "<TR><TD><FONT SIZE=\"-1\">$_</FONT></TD></TR>\n";
    }
    print OUTHTML "</TABLE></P>\n";

    # Subsequence problem
    print OUTHTML "<TABLE BORDER=\"0\" WIDTH=\"100%\">\n";
    print OUTHTML "<TH BGCOLOR=\"darkblue\"><TD><FONT COLOR=\"white\">Subsequences</FONT></TD></TH>\n";
    foreach (@subsequence) {
	print OUTHTML "<TR><TD><FONT SIZE=\"-1\">$_</FONT></TD></TR>\n";
    }
    print OUTHTML "</TABLE></P>\n";

    print OUTHTML $HTML_END;
    close OUTHTML;

    undef (@finished);
    undef (@annotated);
    undef (@date);
    undef (@sequence);
    undef (@subsequence);
    undef (@link);
    
}

#######################################################################
# Help and error trap outputs                                         #
#######################################################################

sub run_details {
    print "# dbcheck\n";     
    print "# version        : $version\n";
    print "# run details    : $rundate $runtime\n";
    print "\n";
}


sub usage {
    my $error = shift;

    my $WormBase_release_file;

    if ($error == 1) {
        my $Wormbase_release_file;
        # No WormBase release number file
        print "The WormBase release number cannot be parsed\n";
        print "Check File: '$Wormbase_release_file'\n\n";
        &run_details;
        exit(0);
    }
    elsif ($error == 0) {
        # Normal help menu
        exec ('perldoc',$0);
        exit (0);
    }
}


__END__

=pod

=head1 NAME - dbcheck

=back

=head1 USAGE

=over 4

=item dbcheck [-options]

dbcheck performs a number of integrity/consistency checks against
the dbace database. The script is based on an iterative loop across
all Genome_sequences and LINK* objects.

=back

=head2 dbcheck MANDATORY arguments:

=over 4

=item none

=back

=head2 dbcheck OPTIONAL arguments:

=over 4

=item -h, Help

=item -d, Debug/Verbose mode

=back

=head1 DOCUMENTATION

=over 4

=back

The following checks have been incorporated into dbcheck:

=head2 Status tags

=head3 Genome sequences which are not Finished.

Genome sequences which do not have a Finished tag.

=head3 Genome sequences which are Finished but not Annotated.

Genome sequences which are finished but not annotated.

=head2 File storage on /analysis/cosmids

=head3 Date mismatch between the file system and dbace.

Inconsistent Date_directory tag in ACEDB with respect to the file
system (/analyis/cosmids/current.versions).

For details of how the date dirtectory structure works:
 http://intweb.sanger.ac.uk/Projects/C_elegans/MANUAL

=head3 Sequence mismatch between the file system and dbace. 

This is based on a GCG checksum calculation for the .seq file in 
the date directory and the sequence extracted from ACEDB.

For details of how the date dirtectory structure works:
 http://intweb.sanger.ac.uk/Projects/C_elegans/MANUAL

=head2 Gene Models

=head3 Absence of Source in ?Sequence

All Subsequence Gene Models MUST have a parent ?Sequence.

    i.e. Sequence "ZK637.5"
         Source "ZK637"

=head3 Absence of Source_Exons in ?Sequence

All Subsequence Gene Models MUST have Source_Exons.

    i.e. Sequence "ZK637.5"
         Source_Exons      1   434
                         483   741
                         950  1159
                        1288  1413

=head3 Absence of Method in ?Sequence

All Subsequence Gene Models MUST have a Method tag. Method tags 
are used in two ways: Firstly, as a means of tagging Gene Models
into classes (e.g. 'curated' for active CDS prediction, 'RNA' for
RNA gene) and secondly as an internal ACEDB description of how
to display the object in the F-map.

    i.e. Sequence "ZK637.5"
         Method "curated"

For details of permissible Method tags see:
 http://intweb.sanger.ac.uk/Projects/C_elegans/MANUAL

=head3 Absence of Species tag in ?Sequence

All Subsequence Gene Models MUST have Species tag.

    i.e. Sequence "ZK637.5"
         Species "Caenorhabditis elegans"

=head3 Absence of CDS tag in ?Sequence

All Subsequence Gene Models MUST have a Coding CDS tag.

    i.e. Sequence "ZK637.5"
         Coding CDS

=head3 Mismatch between Parent coordinates and CDS span length

The coordinate span of the Gene Model (1 => n) based on the 
Source_Exon tags MUST be in agreement with the coordinate span
of the Gene Model as a Subsequence of the parent ?Sequence.

    i.e. Sequence "ZK637"
         Subsequence ZK637.5 11124 12536

         Sequence "ZK637.5"
         Source_Exons      1   434
                         483   741
                         950  1159
                        1288  1413

         Span in CDS (ZK637.5) = ( 1413 -     1) + 1 = 1413
         Parent (ZK637)        = (12536 - 11124) + 1 = 1413
                                                       ----


=head2 LINK objects

=head3 CDS Gene Models on SUPERLINK objects

CDS Gene Models are not allowed as children of SUPERLINK objects.

=head1 AUTHOR - Daniel Lawson

Email dl1@sanger.ac.uk

=cut







