#!/usr/local/bin/perl
#
# tRNApatch.pl v0.1
# dl
# 1999-11-24
#
# Aceperl script to assign Source_exons based on parent Subsequence coordinates
#
###############################################################################
#
# 991124 dl : PP version
# 000315 rd : added $path variable, getopt, and -p option for named $path


########################################
# iniatialise                          #
########################################

use Ace;
use Getopt::Std;

########################################
# command-line parsing                 #
########################################

getopts ('csdl:');

if ($opt_c) { $path = glob("~/acedb/ace4/cam") ; }
elsif ($opt_s) { $path = glob("~/acedb/ace4/stl") ; }
elsif ($opt_l) { $path = glob($opt_l) ; }
else { warn "You must use option -c or -s or -l\n" ; &usage ; }

$debug = $opt_d ;

my $file = shift;
if ($file eq "") {&usage;}

########################################
# usage subroutine                     #
########################################

sub usage {
    warn "Usage: tRNApatch.pl [-options] <filename>\n";
    warn "  Options:\n";
    warn "    -c  Use Cambridge database  \n";
    warn "    -s  Use St Louis database   \n";
    warn "    -l  <path> Use local database path  \n";
    warn "    -d  Debug/Verbose mode\n\n";
    exit;
}

########################################
# Connect with acedb database          #
########################################

$|=1;
warn  "Opening database $path ....\n";

$db = Ace->connect(-path=>$path) || 
    do { warn "Connection failure: ",Ace->error; &maillog; die(); } ;

warn "Connection OK.\n\n";

########################################
# Main Loop                            #
########################################

open (FILE, "$file") || die "Can't open input file\n";
while (<FILE>) {

    if (/Sequence : \"(\S+)\"/) {
	$gene = $1;
	($project,$n) = split (/\./, $gene);
	if ($debug == 1) {warn "Parsing tRNA $gene\n";}
	&get_subsequence_coordinates;
    }
}
close FILE;
exit;

########################################
# Subroutines                          #
########################################

sub get_subsequence_coordinates {

    $obj = $db->fetch(Sequence=>$project);
    if (!defined ($obj)) {
	warn "Could not fetch sequence $project\n";
	next;
    }

    foreach $child ($obj->Subsequence) {
	($seq, $start, $end) = $child->row();
	if ($debug) {warn "$seq\t$start -> $end\n";}
	$diff = $end - $start;
	if ($diff < 0) {
	    $diff = $start - $end + 1;
	}

	if ($seq eq $gene) {
	    print "Sequence : \"$gene\"\nSource_Exons 1 $diff\n";
	    print "Method tRNAscan-SE-1.11\n\n";
	}
    
    $child->DESTROY();
    $obj->DESTROY();
    }
}

__END__

=pod

=head2   NAME - tRNApatch.pl

=head1 USAGE

tRNApatch will reconstruct the Source_exons for tRNA genes based on the 
Subsequence data in the parent object.


tRNApatch.pl Mandatory arguments:

=over 3

=item  -c  Use Cambridge database

=item  -s  Use St Louis database

=item  -l  <path> Use local database path

=back


tRNApatch.pl OPTIONAL arguments:

=over 4

=item  -d  Debug/Verbose mode

=back



