#!/usr/local/bin/perl

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

use Ace ;
use Getopt::Std;
use Cwd;

getopts ('cal:');

my $CWD = cwd;
$ENV{PATH}="/nfs/disk100/wormpub/ACEDB/bin.ALPHA_4:$ENV{PATH}";

if ($opt_c) {
  $db = Ace->connect(-path=>'/nfs/disk100/wormpub/acedb/ace4/cgc') or die ("Could not connect with cgcace\n");
} elsif ($opt_a) {
  $db = Ace->connect(-path=>'/wormsrv2/autoace') or die ("Could not connect with autoace\n");
} elsif ($opt_l) {
  if ($opt_l =~ /^(\~\w+)\//){
    $TARGETDIR=glob("$1");
    $TARGETDIR =~ s/\/tmp_mnt//;
    $FILENAME=$';
    $TARGETDIR="$TARGETDIR"."/"."$FILENAME";
 } elsif ($opt_l =~ /^(\w+)/) {
    $TARGETDIR="$CWD"."/"."$opt_l";
  } elsif ($opt_l =~ /\/\w+/) {
    $TARGETDIR=$opt_l;
  } else {
    print "Incorrect database path\n";
  }
  if (!-d $TARGETDIR) {
    die ("Directory $TARGETDIR not existent\n");
  }
  $db = Ace->connect(-path=>"$TARGETDIR");
} else {
  &PrintHelp;
}


print "\nSequence CHROMOSOME_I\nMethod Link\n" ; $pos = 1 ;
&add ("cTel33B") ; &overlap ("cTel33B") ;
&add ("SUPERLINK_RW1") ; &overlap ("C30F12") ;
&add ("SUPERLINK_CB_I") ; &overlap ("H10E24") ;
&add ("SUPERLINK_RW1R") ; &overlap ("F49D11") ;
&add ("SUPERLINK_CB_IR") ;

print "\nSequence CHROMOSOME_II\nMethod Link\n" ; $pos = 1 ;
&add ("LINK_cTel52S") ; &overlap ("2L52") ;
# [010221 dl] Modified to reflect new cTel clone on the left end
&add ("SUPERLINK_RW2") ; &overlap ("C06A8") ;
&add ("SUPERLINK_CB_II") ;

print "\nSequence CHROMOSOME_III\nMethod Link\n" ; $pos = 1 ;
&add ("cTel54X") ; &overlap ("cTel54X") ;
&add ("SUPERLINK_RW3A") ; &overlap ("Y53G8AR") ;   
# [000911 dl] Modified to reflect new YAC structure LRM nomenclature 
&add ("SUPERLINK_CB_IIIL") ; &overlap ("C38D4") ;
&add ("SUPERLINK_RW3B") ; &overlap ("PAR3") ;
&add ("SUPERLINK_CB_IIIR") ;

print "\nSequence CHROMOSOME_IV\nMethod Link\n" ; $pos = 1 ;
&add ("cTel4X") ; &overlap ("cTel4X") ;
# [010221 dl] Unchanged but modified tag in cTel4X should work 
&add ("SUPERLINK_RW4") ; &overlap ("H23L24") ;
&add ("SUPERLINK_CB_IV") ;

print "\nSequence CHROMOSOME_V\nMethod Link\n" ; $pos = 1 ;
&add ("cTel3X") ; &overlap ("cTel3X") ; 
# [000509 dl] No longer needed, the Overlap_right tag in cTel3X is sufficient
#$pos -= 147 ;	# because cTel3X -> B0348 which is at 148
&add ("SUPERLINK_RW5") ; &overlap ("H24G06") ;
&add ("SUPERLINK_CB_V") ;

print "\nSequence CHROMOSOME_X\nMethod Link\n" ; $pos = 1 ;
&add ("cTel7X") ; &overlap ("cTel7X") ;
&add ("SUPERLINK_RWXL") ; &overlap ("C23F12") ;
&add ("SUPERLINK_CB_X") ; &overlap ("C11G6") ;
&add ("SUPERLINK_RWXR") ; &overlap ("H11L12") ;
&add ("LINK_6R55") ;

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

#---------------------------
# Prints help and disappears
#
sub PrintHelp {
   exec ('perldoc',$0);
}


__END__

=pod

=head1 NAME - makeChromLinks.pl

=head2 DESCRIPTION

makeChromLinks.pl will build and display an ace file for CHROMOSOME* 
objects from autoace or cgcace, depending from the chosen switch.

=head2 MANDATORY arguments (one of the following):

=over 4

=item -a reads from autoace

=item -c reads from cgcace

=item -l database directory

=back

The database directory can be written in one of the following manners:

=over 2

=item ~username/physical/path,

=item /physical/path 

=item directoryname, when it is a subdirectory of the current user directory

=cut




