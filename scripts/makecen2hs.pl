#!/usr/local/bin/perl
# Makecen2hs.pl makes cen2hs.ace querying Alan Coulson's
# cen2hs database and using Fred Wobus' contig_dumpace.
# Must run from an IRIX 6.5 machine (rathbin)
# [ag3 991126]

$|=1;

use File::Find;
$cemapdir = glob "~cemap";
$pmapdir = glob "~wormpub/acedb/ace4/autoace/physical_map";
$sylviadir = glob "~sylvia/homonyms";
$oldaces = glob "~wormpub/acedb/ace4/autoace/physical_map/*.ace";
$log = "$pmapdir/makecen2hs.log";
unlink $log;

open (LOG,">$log");
# Pick up all the old .gz files
# Erase them when they are older than 5 days
@DIRLIST=($pmapdir);
find sub {
  $filename = $File::Find::name;
  if (!-d $filename){
    if (($filename =~ /orig\.gz/)||($filename =~ /ace\.gz/)) {
    if (-M $filename > 5) {
      print LOG "$filename is older than 5 days - erased\n";
      unlink $filename;
    }
  }
  }
},@DIRLIST;

$realdate=&GetTime();
$date1="$pmapdir/$realdate.1";
$dateorig="$pmapdir/$realdate.orig";
$acefile="$pmapdir/$realdate.ace";
$cen2hs="$pmapdir/cen2hs.ace";

print LOG "Realdate: $realdate\n";
print LOG "Date1: $date1\n";
print LOG "Dateorig: $dateorig\n";
print LOG "Acefile: $acefile\n";
print LOG "Cen2hs: $cen2hs\n";

$dmp=system("$pmapdir/contig_dumpace $cemapdir/cen2hs $dateorig");
printf LOG ("Dump:: return value %d - signal %x \n", $dmp >> 8, $dmp & 0xff) ;
if (($dmp >> 8)!=0) {die ("Error in contig_dumpace");}
open (OUT,">$date1");
open (GREP,"grep -v \\\"\\\" $dateorig |") or die ("makecen2hs.pl: could not grep $dateorig\n");
while (<GREP>) {
	print OUT $_;
}
close GREP;
close OUT;

open (OUT,">$acefile") or die ("makecen2hs.pl: could not open $acefile\n");
open (OMONYM,"$pmapdir/homonym.p -T Positive_locus $sylviadir/homonym.gene $date1 |") or die ("makecen2hs.pl: could not do homonym.p\n");
while (<OMONYM>) {
	print OUT $_;
}
close OMONYM;
close OUT;

unlink $date1;
$gzp=system("gzip -f $dateorig");
printf LOG ("Gzip (1): return value %d - signal %x \n", $gzp >> 8, $gzp & 0xff) ;
if (($gzp >> 8)!=0) {warn ("Error in gzipping $dateorig:");}

$gzp=system("gzip -f $oldaces");
if (($gzp >> 8)!=0) {warn ("Error in gzipping $dateorig:");}
printf LOG ("Gzip (2): return value %d - signal %x \n", $gzp >> 8, $gzp & 0xff) ;

unlink $cen2hs;
`ln -s $acefile $cen2hs`;
print LOG ("makecen2hs.pl correctly executed at $realdate\n");

close LOG;

&email_log;

exit;

#------------------------
# Get time coordinates
#
sub GetTime {
  my ($SECS,$MINS,$HOURS,$DAY,$MONTH,$YEAR)=(localtime)[0,1,2,3,4,5];
  if ($MINS=~/^\d{1,1}$/) {
    $MINS="0"."$MINS";
  }
  my $REALMONTH=$MONTH+1;
  my $REALYEAR=$YEAR+1900;
  my $FILEDATE = "$REALMONTH"."_$DAY"."_$REALYEAR";
  return $FILEDATE;
} # end GetTime


#------------------------------
# Email the log file to wormpub
#
sub email_log {
  system ("mailx -s \"makecen2hs.pl log \"  wormpub\@sanger.ac.uk < $log");
  return;
}

