#!/usr/local/bin/perl
#
# script to read in a new version of a sequence from a new (finished) data directory
#
# RD 990818
#
#
# 001018 dl : modified $camdir path to /wormsrv2/camace
# 001018 dl : removed check and dbwrite for _1.ace (no longer used)

use Ace;

############################################
# set up some variables

$camdir = "/wormsrv2/camace" ;
$tmpdir = glob "~wormpub/tmp" ;
$cosdir = glob "~wormpub/analysis/cosmids" ;

use Sys::Hostname ; $host = hostname() ;

############################################

$seq = shift @ARGV || die "Usage: updateCam.pl <seqname>" ;

$version = `grep $seq $cosdir/current.versions` || die "can't find $seq in current.versions\n" ;
chomp $version ;
($s,$date) = split(/\//, $version) ;

$seqdir = "$cosdir/$seq" ;
$datedir = "$seqdir/$date" ;

if (! -d $datedir) { die "date directory $datedir missing\n" ; }

# check Date_directory contains required .ace files

$acefile = "$datedir/$seq" . "_0.ace" ; if (! -e $acefile) { die "file $acefile missing\n" ; }
#$acefile = "$datedir/$seq" . "_1.ace" ; if (! -e $acefile) { die "file $acefile missing\n" ; }
$acefile = "$datedir/$seq" . "_2.ace" ; if (! -e $acefile) { die "file $acefile missing\n" ; }

###############
# open database

$db = Ace->connect (-path => $camdir) || die "can't open camace\n" ;

$response = $db->raw_query("save") ;
if ($response =~ /(\w+) \(session (\d+)\) already has write access/) {
    die "can't get write access - $1 locking session $2\n" ;
}

# $status = $db->status; $status{"write"} || die "can't get write access to camace\n" ;
# above doesn't seem to work with direct path rather than server/client

###############################
# get object and do some checks

$seqobj = $db->fetch(Sequence => $seq) || die "can't find $seq in database\n" ;

$acedate = $seqobj->Date_directory(1) ;
if ($acedate) {
    if ($acedate eq $date) {
	print "camace $seq Date_directory $acedate matches current.versions - aborting\n" ;
	exit 1 ;
    }
    print "updating $seq from $acedate to $date\n" ;
}

if ($seqobj->Flipped(0)) { print "sequence $seq is currently flipped\n" ; }

$oright = $seqobj->Overlap_right(1) ;
$oright_num = $seqobj->Overlap_right(2) ;
if ($oright_num) {
    print "Current Overlap_right offset is $oright_num\n" ;
}

@subseq = $seqobj->Subsequence ;
@checked = grep { /^$seq\.[1-9]/ } @subseq ;
if (@checked > 0) { 
    print "Annotated subsequences @checked\n" ;
    die "Sorry - do not deal with annotated subsequences yet\n" ;
}

@children = $seqobj->S_Child(2) ;
@checked = grep { /^$seq\.[1-9]/ } @children ;
if (@checked > 0) { 
    print "Annotated children @checked\n" ;
    die "Sorry - do not deal with annotated children yet\n" ;
}

#########################################################################
# save potentially important data that may be affected by coordinate changes

open (LOG, ">$seqdir/update.$acedate-$date") || die "can't open update file\n" ;

# Confirmed (and predicted?) subsequence objects

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time) ;
$now = "$mday"."/"."$mon"."/"."$year".":"."$hour".":"."$min" ;

print LOG "camUpdate.pl $now\n\n" ;

print LOG "updating sequence $seq from acedate $acedate to current version $date\n\n" ;

print LOG "current length " . $seqobj->DNA(2) ."\n" ;
$composition = `composition $datedir/$seq.seq` ;
if ($composition =~ /(\d+) total/) { print LOG "new length $1\n" ; }
print LOG "current Overlap_right $oright $oright_num\n" ;

foreach $obj (@subseq) {
    if ($obj =~ /^$seq\.[1-9a-z]/) {
	print LOG "\n$obj " . $obj->asAce . $obj->fetch->asAce . $obj->asDNA ;
    }
}
foreach $obj (@children) {
    if ($obj =~ /^$seq\.[1-9a-z]/) {
	print LOG "\n$obj " . $obj->asAce . $obj->fetch->asAce . $obj->asDNA ;
    }
}

# Confirmed_intron UTR data
# Allele data

close LOG ;

###############################################
# now delete old information and add in the new

# first remove data from the object itself

$db->parse(<<END) ;
Sequence $seq
-D DNA
-D Source
-D S_Parent
-D Clone_left_end
-D Clone_right_end
-D Remark "unfinished fragment"
-D Assembly_tags
-D Splices
-D Allele
-D Analysis_details
-D Homol
-D Feature
END

if ($oright_num) {
    $db->parse(<<END) ;
Sequence $seq
-D Overlap_right $oright $oright_num
Overlap_right $oright
END
}

# next kill any subsequences/children

foreach $obj (@subseq) { if (! defined $obj->kill) { die "failed to kill subsequence $obj\n" ; } }
foreach $obj (@children) { if (! defined $obj->kill) { die "failed to kill child $obj\n" ; } }

# finally read in new analysis files

if (! $db->parse_file ("$datedir/$seq" . "_0.ace")) { die $db->error() . "\n" ; }
#if (! $db->parse_file ("$datedir/$seq" . "_1.ace")) { die $db->error() . "\n" ; }
if (! $db->parse_file ("$datedir/$seq" . "_2.ace")) { die $db->error() . "\n" ; }

$db->raw_query("save") ;

exit 0 ;

############### end of file
