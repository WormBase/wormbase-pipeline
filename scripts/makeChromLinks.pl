#!/usr/local/bin/perl5.8.0 -w

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
# Last updated by: $Author: ar2 $                     
# Last updated on: $Date: 2004-06-28 13:12:17 $       

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Getopt::Long;
use Cwd;



##############################
# command-line options       #
##############################

my $help;       # Help perdoc
my $database;   # Database name for single db option
my $debug;      # Debug mode, verbose output to runner only
my $test;       # test mode, uses ~wormpub/TEST_BUILD
my $out;
GetOptions (
	    "database:s"  => \$database,
	    "debug=s"     => \$debug,
	    "help"        => \$help,
	    "test"        => \$test,
	    "out:s"       => \$out
	    );

# help page
&usage("Help") if ($help);


my $maintainers = "All";

# debug
if($debug){
  print "// makeChromLinks run with debug as $debug\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

# database/file paths and locations
my $basedir     = "/wormsrv2";
$basedir        = glob("~wormpub")."/TEST_BUILD" if ($test); 



 ##############################
 # Script variables (run)     #
 ##############################


my $rundate     = &rundate;
my $runtime     = &runtime;

# touch logfile for run details
$0 =~ m/\/*([^\/]+)$/; system ("touch $basedir/logs/history/$1.`date +%y%m%d`");
my $logfile = "$basedir/logs/$1.`date +%y%m%d`.$$";



# where am i
my $CWD = cwd;

if (!defined $database) {
    $database = "$basedir/autoace";
}

print "// Using $database as source of data for chromosomes\n" if ($debug);
print STDERR "\t\tconnecting to $database\n";
# AcePerl connection to $database
my $db = Ace->connect(-path=>$database,
                      -program =>&tace) or die ("Could not connect with $database\n");
print STDERR "\t\tConnected\n";
print "Connected to database\n" if ($debug);

my ($pos,$i);

our $fh;
open ($fh, ">$out") or die "$out asfgasdg";


print $fh "\nSequence CHROMOSOME_I\nMethod Link\n" ; $pos = 1 ;
&add ("SUPERLINK_RW1");     &overlap ("C30F12") ;
&add ("SUPERLINK_CB_I");    &overlap ("H10E24") ;
&add ("SUPERLINK_RW1R");    &overlap ("F49D11") ;
&add ("SUPERLINK_CB_IR");

print $fh "\nSequence CHROMOSOME_II\nMethod Link\n" ; $pos = 1 ;
&add ("SUPERLINK_RW2");     &overlap ("C06A8") ;
&add ("SUPERLINK_CB_II");   &overlap ("Y53F4B");
&add ("SUPERLINK_RW2R");

print  $fh "\nSequence CHROMOSOME_III\nMethod Link\n" ; $pos = 1 ;
&add ("SUPERLINK_RW3A");    &overlap ("Y53G8AR") ;   
&add ("SUPERLINK_CB_IIIL"); &overlap ("C38D4") ;
&add ("SUPERLINK_RW3B");    &overlap ("PAR3") ;
&add ("SUPERLINK_CB_IIIR");
 
print  $fh "\nSequence CHROMOSOME_IV\nMethod Link\n" ; $pos = 1 ;
&add ("SUPERLINK_RW4");     &overlap ("H23L24") ;
&add ("SUPERLINK_CB_IV");

print  $fh "\nSequence CHROMOSOME_V\nMethod Link\n" ; $pos = 1 ;
&add ("SUPERLINK_RW5");     &overlap ("H24G06") ;
&add ("SUPERLINK_CB_V");

print $fh "\nSequence CHROMOSOME_X\nMethod Link\n" ; $pos = 1 ;
&add ("SUPERLINK_RWXL");    &overlap ("C23F12") ;
&add ("SUPERLINK_CB_X");    &overlap ("C11G6") ;
&add ("SUPERLINK_RWXR"); 

$db->close;
print STDERR "\tclosing database\n"; 
 ###############
 # hasta luego #
 ###############
exit 0 ;

############################################################

sub add {
    my ($seq) = @_ ;
    my $obj = $db->fetch(Sequence=>$seq) ;
    my $length = 0 ;
    print STDERR "\t\t\tadding $seq\n";

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
	print $fh "Subsequence $seq $end $pos // Flipped\n" ;
    } else {
	print $fh "Subsequence $seq $pos $end\n" ;
    }
    $pos = $end + 101 ;		# NB modify global
    print STDERR "\t\t\tadded $seq\n";
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

=item -test, uses test environment in ~wormpub/TEST_BUILD

=back

The database directory can be written in one of the following manners:

=over 2

=item ~username/physical/path,

=item /physical/path 

=item directoryname, when it is a subdirectory of the current user directory

=cut






