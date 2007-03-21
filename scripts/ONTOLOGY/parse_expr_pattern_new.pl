#!/nfs/disk100/wormpub/bin/perl -w
          
use lib $ENV{'CVS_DIR'};
use strict;            
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

my ($ep, $ref, %genes, %at, %auth, $date);

$|=9;

my $program_name=$0=~/([^\/]+)$/ ? $1 : '';
my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($output, $acedbpath);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	  	  	"test"       => \$test,
	    	"verbose"    => \$verbose,
	    	"store:s"    => \$store,
	    	"database:s" => \$acedbpath,
	    	"output:s" 	 => \$output
	    );

if ($help) {
    print "usage: $program_name -o output -d database -p tace\n";
    print "       -help              help - print this message\n";
    print "       -output <output>     output file\n";
    print "       -database <database>   path to AceDB\n";
    exit;
}

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$acedbpath=$wormbase->autoace unless $acedbpath;
my $tace= $wormbase->tace;


warn "connecting to database... ";
my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or $log->log_and_die("Connection failure: ". Ace->error);
warn "done\n";

my %names=();
my @aql_results=$db->aql("select a, a->public_name from a in class gene");
foreach (@aql_results) {
    $names{$_->[0]}=$_->[1];
}
warn scalar keys %names , " genes read\n";

my $out;
$output = $wormbase->ontology."/anatomy_association.".$wormbase->get_wormbase_version_name.".wb" unless $output;
open($out, ">$output") or $log->log_and_die("cannot open $output : $!\n");	

my $query="find Expr_pattern Anatomy_term";
my $it=$db->fetch_many(-query=>$query);

my $count=0;
while (my $obj=$it->next) {
    $count++;
    if ($count % 1000 == 0) {
	warn "$count objects processed\n";
    }
    my @lines=split("\n", $obj->asAce);
    foreach (@lines) {
	if (/^Expr_pattern/) {
	    $ep=$_=~/\"(.+)\"/ ? $1 : '';
	}
	elsif (/^Reference/) {
	    $ref=$_=~/\"(WBPaper\d+)\"/ ? $1 : '';
	}
	elsif (/^Expression_of\s+Gene/) {
	    my $tmp=$_=~/\"(WBGene\d+)\"/ ? $1 : '';
	    next unless $tmp;
	    $genes{$tmp}++;
	}
	elsif (/^Expressed_in\s+Anatomy_term/) {
	    my $tmp=$_=~/\"(WBbt:\d+)\"/ ? $1 : '';
	    next unless $tmp;
	    $at{$tmp}[0]++;
	    if (/\s+Certain/i) {
		$at{$tmp}[1]='Certain';
	    }
	    elsif (/\s+Uncertain/i) {
		$at{$tmp}[1]='Uncertain';
	    }
	    elsif (/\s+Partial/i) {
		$at{$tmp}[1]='Partial';
	    }
	    
	}
	elsif (/^Experiment\s+Author/) {
	    my $tmp=$_=~/\"(.+)\"/ ? $1 : '';
	    next unless $tmp;
	    $auth{$tmp}++;
	}
	elsif (/^Experiment\s+Date/) {

	}
	
    }

    if ($ep) {
	foreach my $g (keys %genes) {
	    foreach my $a (keys %at) {
		my $q='';
		if ($at{$a}[1]) {
		    $q=$at{$a}[1];
		}
		print $out "WB\t$g\t$names{$g}\t$q\t$a\t$ref\tExpr_pattern\t$ep\t\t\t\t\t\t$date\tWB\n";
	    }
	}
    }
    $ep='';
    $ref='';
    %genes=();
    %at=();
    %auth=();
}

$log->mail;
exit;


