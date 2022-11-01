#!/usr/local/bin/perl5.6.1 -w
#
# acezip.pl : an improved acecompress.pl
#
# dl
#
# Usage : acezip.pl [-options]
#
# Last edited by: $Author: gw3 $
# Last edited on: $Date: 2009-09-17 09:42:46 $
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
# use the time and the process ID to make a unique file extension
my $time = time();
my $pid = "$$";
my $tmpDir= $ENV {'WB_SCRATCH'};
my $outfile = "/$tmpDir/acezip.$pid.$time.new";
my $presort = "$file.$pid.$time.presort";
my $sort = "$file.$pid.$time.sort";

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
open (FILE2, ">$presort") or $log->log_and_die("cant open $presort : $!\n");
while (my $record = <FILE>) {
  while ($record =~ /^\s*\/\//) {$record =~ s/^\s*\/\/.*?\n//} # strip out any comments at the start of the record
  $record =~ s/\t/ /g;
  $record =~ s/\n/\t/g;
  print FILE2 $record, "\n";
}
close FILE;
  
# reset input line separator
$/= $oldlinesep;

# sort the file - sometimes the sorted file doesn't appear - this is very odd - try a few times to make it
my $tries = 5;
while ($tries-- && ! -e "$sort") {
  $wormbase->run_command("sort -S 4G $presort -o $sort", $log);
  system('sleep 5'); # wait a few seconds for NFS to realise that there really is a file there
}

# print output

open( ACE, "<$sort") or $log->log_and_die("cant read from $sort : $!\n");
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
$wormbase->run_command("rm -f $file $sort $presort", $log);
$wormbase->run_command("mv -f $outfile $file", $log);

$log->mail;
exit(0);
