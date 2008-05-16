#/software/bin/perl -w
#
# make_acefiles.pl 
#
# dl1
#
# Generates the .acefiles from the primary databases as a prelim for building
# autoace.
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2008-05-16 12:52:47 $

#################################################################################
# Variables                                                                     #
#################################################################################

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use File::Copy;
use File::Spec;
use File::Path;
use Log_files;
use Storable;

##############################
# command-line options       #
##############################

our $help;       # Help perdoc
our $debug;      # Debug mode, verbose output to runner only
my $test;        # If set, script will use TEST_BUILD directory under ~wormpub
my $db;
my $database;
my $basedir;
my $store;
my $species;
my ($syntax, $merge);

GetOptions (	"debug=s"    => \$debug,
		"help"       => \$help,
		"database=s" => \$database,
		"db:s"       => \$db,
		"test"       => \$test,
		"store:s"    => \$store,
		"species:s"	 => \$species,
		"syntax"	 => \$syntax,
		"merge"		 => \$merge
	   	);

my $wormbase;
if( $store ) {
    $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
    $wormbase = Wormbase->new( -debug   => $debug,
			       -test    => $test,
			       -organism => $species
			       );
}

my $log = Log_files->make_build_log($wormbase);
my $elegans = $wormbase->build_accessor;

my $config = $wormbase->basedir."/autoace_config/".$wormbase->species;
$config .= ".merge" if $merge;
$config .=".config";

my $tace = $wormbase->tace;
my $path = $wormbase->acefiles.($merge ? "/MERGE" :"/primaries");
mkpath $path unless -e $path;

unless (-e $config) {
	$log->write_to("merge config file absent - database being skipped\n");
}
else {
	open (CFG,"<$config") or $log->log_and_die("cant open $config :$!\n");

	my $dbpath = "";
	while(<CFG>) {
		next if /#/;
		next unless /\w/;
    	my %makefile;
    	foreach my $pair (split(/\t+/,$_)) {
			my($tag,$value) = ($pair =~ /(\w+)=(.*)/);
			if( $tag and $value) {
				if($tag eq 'delete') {
					push(@{$makefile{$tag}},$value);
				}else {
		    		$makefile{$tag} = $value;
		    	}
			}
			else {
		    	$log->log_and_die("Ill formed config line $_ :\n");
		    }
    	}
    	next if $syntax;
    	my $query = "nosave\nquery ";
    	if( $makefile{'class'} and  $makefile{'db'} and  $makefile{'file'} ) {
    		if ($db) {
		    	next unless ($makefile{'db'} eq $db);
    		}
    		mkpath("$path/".$makefile{'db'}) unless -e "$path/".$makefile{'db'};
    		my $file = $path."/".$makefile{'db'}."/".$makefile{'file'};
    		open(ACE,">$file") or $log->log_and_die("cant open file $file : $!\n");
    		
    		if($makefile{'class'} eq 'DNA') {
 				$query .= "find Sequence";
 			}else {
		    	$query .= "find ".$makefile{'class'};
		    }
		    if( $makefile{'query'} ) {
		    	$query .= " WHERE ".$makefile{'query'};
			}
			$query .= "\n";    
 			if( $makefile{'delete'} ) {	 
 				foreach my $del (@{$makefile{'delete'}}) {   
		 			$query .= "eedit -D $del\n";
		 		}
 			}
 			if($makefile{'class'} eq 'DNA') {
 				$query .= "FOLLOW DNA\n";
 			}
			if( $makefile{'follow'} ) {
				$query .= "show -T ".$makefile{'follow'}." -a\n"; #output parent + followed tag
				$query .= "FOLLOW ".$makefile{'follow'}."\n";
			}
			$query .= "show -T -a ";
			if( $makefile{'tag'} ) {
				$query .= " -t ".$makefile{'tag'};
			}
			$query .= "\n";
			my $acedb = $dbpath."/".$makefile{'db'};
			$log->write_to("dumping $makefile{'class'} from $acedb\n");
			open(TACE,"echo '$query' | $tace $acedb | ") or $log->log_and_die("cant do query : $!\n");
			LINE: while(<TACE>) {
				next if (/acedb>/ or /\/\//);
				if( $makefile{'regex'} ) {	
					unless (/[^\w]/ or /$makefile{'class'}\s+\:\s+/ or /$makefile{'follow'}\s+\:\s+/) {
						next LINE unless /$makefile{'regex'}/;
					}			
    			}
    			print ACE;
    		}
    		close TACE;
    		close ACE;
		}
		elsif ($makefile{'path'}) {
			my $sub = $makefile{'path'};
			$dbpath = $wormbase->$sub;
		}
	}
}	

$log->write_to("syntax ok\n") if $syntax; #will have died before here if not ok
$log->mail;
	
	
	
	
	