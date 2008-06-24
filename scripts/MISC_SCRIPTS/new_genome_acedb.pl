#!/usr/local/bin/perl -w
use lib $ENV{'CVS_DIR'};

use strict;
use Getopt::Long;
use Wormbase;
use Log_files;
use File::Copy;

my($species,$sequence_file,$dbpath,$init);
my ($test, $debug);
GetOptions (
	    'species:s' => \$species,
	    'fasta:s'   => \$sequence_file,
	    'init'		=> \$init,
	    'dbpath:s'  => \$dbpath,
	    'test' 		=> \$test,
	    'debug:s'	=> \$debug
	   );

my $wormbase = new Wormbase ( -test => $test,
							  -debug=> $debug,
							  -organism => $species,
							 );
							 
my $log = Log_files->make_build_log($wormbase);
$dbpath = $wormbase->wormpub."/DATABASES/".$species unless $dbpath;
my $full_species = $wormbase->full_name or die "cant identify\n";

&initialise if $init;
&load_sequence;

$log->mail;

sub initialise {
	# create directory, get latest CVS of wspec dir and reinit database.
	$wormbase->run_command("mkdir $dbpath", $log) unless -e $dbpath;
	$wormbase->run_command("mkdir $dbpath/database", $log) unless -e "$dbpath/database";
	$wormbase->run_command("mkdir $dbpath/acefiles", $log) unless -e "$dbpath/acefiles";	
	$wormbase->run_command("cd $dbpath".';cvs -d :pserver:cvsuser@cvsro.sanger.ac.uk:/cvsroot/ensembl checkout -d wspec wormbase/wspec', $log);
	$wormbase->run_command("cd $dbpath".';cvs -d :pserver:cvsuser@cvsro.sanger.ac.uk:/cvsroot/ensembl checkout -d autoace_config wormbase/autoace_config', $log);
	&DbWrite('y',$wormbase->tace,$dbpath,"ReInitDB $species");
	
	#copy the genefinder table file
	$wormbase->run_command("cp -r ".$wormbase->database('current')."/wgf/ $dbpath/", $log);
	
	#rename database
	my $database_wrm = "$dbpath/wspec/database.wrm";
	$wormbase->run_command('sed s/WS0/'.$wormbase->full_name('-short' => 1).'/ < '.$database_wrm.' > '.${database_wrm}.'.new',$log);
	my $status = move("${database_wrm}.new", "$database_wrm");
}

sub load_sequence {
	$wormbase->load_to_database($dbpath,$sequence_file,'sequences', $log);
	$wormbase->load_to_database($dbpath,"$dbpath/autoace_config/misc_autoace_methods.ace",'methods',$log);
	$wormbase->load_to_database($dbpath,&create_seq_obj( $sequence_file ),'seq_objs', $log);
}

sub DbWrite {
    my ($command,$exec,$dir,$name)=@_;
    open (WRITEDB,"| $exec $dir") or die "$name DbWrite failed\n";
    print WRITEDB $command;
    close WRITEDB;
}

sub create_seq_obj {
	my $fasta = shift;
	my $acefile = "$dbpath/acefiles/genomic_canonical.ace";
	open( FASTA,"<$sequence_file") or $log->log_and_die("cant open $sequence_file: $!\n");
	open( ACE,">$acefile") or $log->log_and_die("cant write $acefile: $!\n");
	while(<FASTA>) {
		chomp;
		if(/>(\w+)\s?/){
			print ACE "\nSequence : \"$1\"\n";
			print ACE "method Genomic_canonical\n";
			print ACE "Genomic_canonical\n";
			print ACE "Species \"$full_species\"\n";
		}
	}
	return $acefile;
	close FASTA;
	close ACE;
}
	
	
	

