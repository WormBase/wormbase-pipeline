#!/usr/local/bin/perl5.6.1 -w

##########################################################
#
# makeChromLinks is hard-coded to make the worm
# CHROMOSOME* objects from databases which are 
# specified from the command line
#
# RD 990719
#
# use global $pos to say where you are
# ag3 21012000 Added path choice and doc
#
# 010307 : dl  : Added Method tag to chromosome objects
#
##########################################################
#
# Last updated by: $Author: dl1 $                     
# Last updated on: $Date: 2002-12-09 14:56:20 $       

$|=1;
use strict;
use lib "/wormsrv2/scripts/";   
use Wormbase;
use Ace;
use Getopt::Long;
use Cwd;

 ##############################
 # Script variables (run)     #
 ##############################

my $maintainers = "All";
my $rundate     = `date +%y%m%d`;   chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;

# touch logfile for run details
$0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");
my $logfile = "/wormsrv2/logs/$1.`date +%y%m%d`.$$";


 ##############################
 # command-line options       #
 ##############################

my $help;       # Help perdoc
my $database;   # Database name for single db option
my $debug;      # Debug mode, verbose output to runner only

GetOptions (
	    "database:s"  => \$database,
	    "debug=s"     => \$debug,
	    "help"        => \$help
	    );

# help page
&usage("Help") if ($help);

# no debug name
&usage("Debug") if ((defined $debug) && ($debug eq ""));

# assign $maintainers if $debug set
($maintainers = $debug . '\@sanger.ac.uk') if ($debug);

# where am i
my $CWD = cwd;
$ENV{PATH}="/nfs/disk100/wormpub/ACEDB/bin.ALPHA_4:$ENV{PATH}";

if (!defined $database) {
    $database = "/wormsrv2/autoace";
}

# AcePerl connection to $database
my $db = Ace->connect(-path=>'$database') or die ("Could not connect with $database\n");

my ($pos,$i);

print "\nSequence CHROMOSOME_I\nMethod Link\n" ; $pos = 1 ;
&add ("SUPERLINK_RW1");     &overlap ("C30F12") ;
&add ("SUPERLINK_CB_I");    &overlap ("H10E24") ;
&add ("SUPERLINK_RW1R");    &overlap ("F49D11") ;
&add ("SUPERLINK_CB_IR");

print "\nSequence CHROMOSOME_II\nMethod Link\n" ; $pos = 1 ;
&add ("SUPERLINK_RW2");     &overlap ("C06A8") ;
&add ("SUPERLINK_CB_II");   &overlap ("Y53F4B");
&add ("SUPERLINK_RW2R");

print "\nSequence CHROMOSOME_III\nMethod Link\n" ; $pos = 1 ;
&add ("SUPERLINK_RW3A");    &overlap ("Y53G8AR") ;   
&add ("SUPERLINK_CB_IIIL"); &overlap ("C38D4") ;
&add ("SUPERLINK_RW3B");    &overlap ("PAR3") ;
&add ("SUPERLINK_CB_IIIR");

print "\nSequence CHROMOSOME_IV\nMethod Link\n" ; $pos = 1 ;
&add ("SUPERLINK_RW4");     &overlap ("H23L24") ;
&add ("SUPERLINK_CB_IV");

print "\nSequence CHROMOSOME_V\nMethod Link\n" ; $pos = 1 ;
&add ("SUPERLINK_RW5");     &overlap ("H24G06") ;
&add ("SUPERLINK_CB_V");

print "\nSequence CHROMOSOME_X\nMethod Link\n" ; $pos = 1 ;
&add ("SUPERLINK_RWXL");    &overlap ("C23F12") ;
&add ("SUPERLINK_CB_X");    &overlap ("C11G6") ;
&add ("SUPERLINK_RWXR"); 

$db->close;

 ###############
 # hasta luego #
 ###############
exit 0 ;

############################################################

sub add {
    my ($seq) = @_ ;
    my $obj = $db->fetch(Sequence=>$seq) ;
    my $length = 0 ;
    $obj || die "can't find $seq\n" ;
    if ($seq =~ /LINK/) {
	foreach $i ($obj->Subsequence(2)) {
	    if ($i > $length) { $length = $i ; }
	}
	foreach $i ($obj->Subsequence(3)) {
	    if ($i > $length) { $length = $i ; }
	}
    } else {
	$length = $obj->DNA(2) ;
    }
    $length || die "no length for $seq\n" ;
    my $end = $pos + $length - 1 ;
    if ($obj->Flipped(0)) {
	print "Subsequence $seq $end $pos // Flipped\n" ;
    } else {
	print "Subsequence $seq $pos $end\n" ;
    }
    $pos = $end + 101 ;		# NB modify global
}

sub overlap {
    my ($seq) = @_ ;
    my $obj = $db->fetch(Sequence=>$seq) ;
    my $olap = $obj->Overlap_right(2) ;
    if ($olap) {
	my $length = $obj->DNA(2) ;
	$pos += ($olap - $length - 101) ;
    }
}

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

=head1 NAME - makeChromLinks.pl

=head2 DESCRIPTION

makeChromLinks.pl will build and display an ace file for CHROMOSOME* 
objects from autoace or cgcace, depending from the chosen switch.

=head2 MANDATORY arguments (one of the following):

=over 4

=item -database [path] reads from supplied path, default is autoace

=item -debug [name], verbose report 

=item -help, this help page

=back

The database directory can be written in one of the following manners:

=over 2

=item ~username/physical/path,

=item /physical/path 

=item directoryname, when it is a subdirectory of the current user directory

=cut






