#!/nfs/disk100/wormpub/bin/perl -w
          
use lib $ENV{'CVS_DIR'};
use strict;
use Ace;         
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;	

my ($ep, %at, $var);
my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($output, $acedbpath);
my ($rnai, $ref, %genes, %pheno, %auth);

$|=9;


my %opts=();
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	  	  	"test"       => \$test,
	    	"verbose"    => \$verbose,
	    	"store:s"    => \$store,
	    	"database:s" => \$acedbpath,
	    	"output:s"   => \$output
	    );

my $program_name=$0=~/([^\/]+)$/ ? $1 : '';


if ($opts{h}) {
    print "usage: $program_name -output output -database database/n";
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

my $year  = (1900 + (localtime)[5]);
my $month = (1 + (localtime)[4]);
my $day   = (localtime)[3];

my $date=sprintf("%04d%02d%02d", $year, $month, $day);


$acedbpath=$wormbase->autoace unless $acedbpath;
my $tace=$wormbase->tace;


warn "connecting to database... ";
my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or $log->log_and_die("Connection failure: ". Ace->error);
warn "done\n";

my %names=();
my @aql_results=$db->aql("select a, a->public_name from a in class gene");
foreach (@aql_results) {
    $names{$_->[0]}=$_->[1];
}
warn scalar keys %names , " genes read\n";

my %variations=();
@aql_results=$db->aql("select a, a->gene from a in class variation where exists a->gene");
foreach (@aql_results) {
    push @{$variations{$_->[0]}}, $_->[1];
}
warn scalar keys %variations , " variations read\n";





my $out;
if ($output) {
    open($out, ">$output") or $log->log_and_die("cannot open $output : $!\n");
}
else {
    $out=*STDOUT;
}

my $query="find Variation Phenotype";

my $it=$db->fetch_many(-query=>$query);


my $count=0;
while (my $obj=$it->next) {
    $count++;
    if ($count % 1000 == 0) {
	warn "$count objects processed\n";
    }
    my @lines=split("\n", $obj->asAce);
    foreach (@lines) {
	if (/^Variation\s+/) {
	    $var=$_=~/\"(.+)\"/ ? $1 : '';
	}
	elsif (/Phenotype/) {
	    my $tmp=$_=~/\"(WBPhenotype\d+)\"/ ? $1 : '';
	    next unless $tmp;
	    $pheno{$tmp}[0]++;
	    my @tmp1=split('\s+');
	    if ($tmp1[3]) {
		$tmp1[3]=~s/\s+//g;
		if ($tmp1[3] eq 'Not') {
		    $pheno{$tmp}[1]='NOT';
		}
	    }
	    if (/Paper_evidence/) {
		my $paper=$_=~/\"(WBPaper\d+)\"/ ? $1 : '';
		if ($paper) {
		    $pheno{$tmp}[2]=$paper;
		}
	    }
	    
	}
    }
    if ($var) {
	foreach my $g (@{$variations{$var}}) {
	    foreach my $p (keys %pheno) {
		my $q=$pheno{$p}[1] ? 'NOT' : '';
		my $paper=$pheno{$p}[2] ? $pheno{$p}[2] : '';
		print $output "WB\t$g\t$names{$g}\t$q\t$p\t$paper\tVariation\t$var\t\t\t\t\t\t$date\tWB\n";
	    }
	}
    }
    
    $var='';
    %pheno=();
}

$query="find RNAi Phenotype";

$it=$db->fetch_many(-query=>$query);
$count=0;
while (my $obj=$it->next) {
    next unless $obj->isObject();
    $count++;
    if ($count % 1000 == 0) {
	warn "$count RNAi objects processed\n";
    }
    
    my @genes_tmp=$obj->Gene;
    my @genes=();
    my %pheno=();
    foreach (@genes_tmp) {
	if ($_->right(2) eq "RNAi_primary") {
	    push @genes, $_;
	}
    }
    my $ref=$obj->Reference;
    my @phen_array_tmp=$obj->Phenotype;
    foreach (@phen_array_tmp) {
	my $not=grep {/Not/} $_->tags();
	if ($not) {
	    $pheno{$_}[1]='NOT';
	}
	$pheno{$_}[0]++;
    }
    foreach my $g (@genes) {
	foreach my $p (keys %pheno) {
	    my $q=$pheno{$p}[1] ? 'NOT' : '';
	    print $output "WB\t$g\t$names{$g}\t$q\t$p\t$ref\tRNAi\t$obj\t\t\t\t\t\t$date\tWB\n";
	}
    }
}

$log->mail;
exit();

	

