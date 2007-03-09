#!/usr/local/bin/perl5.6.1 -w
#
# acecompress.pl
#
# dl
#
# Usage : acecompress.pl [-options]
#
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2007-03-09 15:15:58 $
use lib $ENV{'CVS_DIR'};

use strict;
use Getopt::Long;
use Storable;
use Wormbase;
use Log_files;

my $homol;         # homol data   BLAT files
my $feature;       # feature data Good_intron files    
my %objects;       # hash for number of homol lines per object
my $ace_object;  
my %acedata;

my ($test, $debug, $store, $file, $bk,$build);

GetOptions (
	    "homol"    => \$homol,
	    "feature"  => \$feature,
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

if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}

my $log = Log_files->make_build_log($wormbase);

open (FILE, "<$file") or $log->log_and_die("cant open $file : $!\n");
while (<FILE>) {

    # BLAT homology data
    if ($homol) {
		if (/Homol_data : \"(\S+)\"/) {
	    $ace_object = $1;
	    $objects{$ace_object}++;
	    next;
	}
	if (/DNA_homol/) {
	    push (@{$acedata{$ace_object}{DNA_homol}},$_);
	}
    }
    # Confirmed_intron feature data
    elsif ($feature) {
	if (/Feature_data : \"(\S+)\"/) {
	    $ace_object = $1;
	    $objects{$ace_object}++;
	    next;
	}
	if (/Confirmed_intron/) {
	    push (@{$acedata{$ace_object}{Confirmed_intron}},$_);
	}
    }
}
close FILE;


# print output

my $outfile = "/tmp/acecompress.$$";

open( ACE,">$outfile") or $log->log_and_die("cant write to $outfile");
	foreach my $obj (keys %objects) {
   	print ACE "\n// $obj\n\n";

   	# homol (BLAT homology)
   	my $line;
   	if ($homol) {
			print ACE "Homol_data : \"$obj\"\n";
			foreach $line (@{$acedata{$obj}{DNA_homol}}) {
	    		print ACE "$line";
			}
    	}
    	# feature (confirmed_intron)
    	elsif ($feature) {
			print ACE "Feature_data : \"$obj\"\n";
			foreach $line (@{$acedata{$obj}{Confirmed_intron}}) {
	    	print ACE "$line";
		}
	}  
}
close ACE;

#replace original, retaining if bk option set
$wormbase->run_command("mv $file $file.bk", $log) if ($bk);
$file=~s/_uncompressed// if $build;
$wormbase->run_command("mv -f $outfile $file",$log);

$log->mail;
exit(0);
