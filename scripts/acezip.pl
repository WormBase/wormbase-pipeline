#!/usr/local/bin/perl5.6.1 -w
#
# acezip.pl : an improved acecompress.pl
#
# dl
#
# Usage : acezip.pl [-options]
#
# Last edited by: $Author: gw3 $
# Last edited on: $Date: 2009-02-11 10:36:55 $
use lib $ENV{'CVS_DIR'};

use strict;
use Getopt::Long;
use Storable;
use Wormbase;
use Log_files;

my ($test, $debug, $store, $file, $bk, $build);

GetOptions (
	    "test"     => \$test,
	    "store:s"  => \$store,
	    "debug:s"  => \$debug,
	    "bk"       => \$bk,       #keep original as .bk
	    "file:s"   => \$file,
	    "build"    => \$build
	    )
		    or die("invalid commandline option\n");
	    ;

my $wormbase;
my $outfile = "/tmp/acezip$$.new";

if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}

my $log = Log_files->make_build_log($wormbase);

# Change input separator to paragraph mode, but store old mode in $oldlinesep
my $oldlinesep = $/;
$/ = "";

open (FILE, "<$file") or $log->log_and_die("cant open $file : $!\n");
open (FILE2, ">$file$$.presort") or $log->log_and_die("cant open $file$$.presort : $!\n");
while (my $record = <FILE>) {
  $record =~ s/\t/ /g;
  $record =~ s/\n/\t/g;
  print FILE2 $record, "\n";
}
close FILE;
  
# reset input line separator
$/= $oldlinesep;

# sort the file
$wormbase->run_command("sort -S 4G $file$$.presort -o $file$$.sort", $log);

# print output

open( ACE, "<$file$$.sort") or $log->log_and_die("cant read from $file$$.sort : $!\n");
open( ACE2, ">$outfile") or $log->log_and_die("cant write to $outfile : $!\n");
my $prev="";
while (my $line = <ACE>) {
  my ($first, $rest) = ($line =~ /(.+?)\t(.+)/);
  if ($first ne $prev) {print ACE2 "\n\n$first\n"; $prev = $first;}
  $rest =~ s/\t/\n/g;
  $rest =~ s/\n$//;
  print ACE2 "$rest";
}
close ACE2;
close ACE;

#replace original, retaining if bk option set
$wormbase->run_command("mv $file $file.bk", $log) if ($bk);
$file=~s/_uncompressed// if $build;
$wormbase->run_command("rm -f $file $file$$.sort $file$$.presort", $log);
$wormbase->run_command("mv -f $outfile $file", $log);

$log->mail;
exit(0);
