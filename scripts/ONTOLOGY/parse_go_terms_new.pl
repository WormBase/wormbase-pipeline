 #!/nfs/disk100/wormpub/bin/perl
          
use lib $ENV{'CVS_DIR'};
use strict;
use Ace;         
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;	

my ($ep, $ref, %genes, %at, %auth, $date);
my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($output, $acedbpath, $rnai, $gene);

$|=9;


my %opts=();
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	  	  	"test"       => \$test,
	    	"verbose"    => \$verbose,
	    	"store:s"    => \$store,
	    	"database:s" => \$acedbpath,
	    	"rnai"   	 => \$rnai,
	    	"gene"		 => \$gene,
	    	"output:s"   => \$output,
	    );

my $program_name=$0=~/([^\/]+)$/ ? $1 : '';


if ($help) {
    print "usage: $program_name [options] -output output -database database\n";
    print "       -help            help - print this message\n";
    print "       -output <output>   output file\n";
    print "       -database <database> path to database\n";
    print "       -rnai            generate RNAi2GO mapping; default no\n";
    print "       -gene            generate gene mapping; default no\n";
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
my %aspect=(molecular_function=>'F',
	    biological_process=>'P',
	    cellular_component=>'C'
	    );
	    
my $year  = (1900 + (localtime)[5]);
my $month = (1 + (localtime)[4]);
my $day   = (localtime)[3];
$date=sprintf("%04d%02d%02d", $year, $month, $day);

$acedbpath = $wormbase->autoace unless $acedbpath;
warn "connecting to database... ";
my $db = Ace->connect(-path => $acedbpath,  -program => $wormbase->tace) or $log->log_and_die("Connection failure: ". Ace->error);
warn "done\n";


my %name_hash=();
my @aql_results=$db->aql("select a, a->public_name from a in class gene");
foreach (@aql_results) {
    $name_hash{$_->[0]}=$_->[1];
}
warn scalar keys %name_hash , " gene public names read\n";


my %seq_name_hash=();
@aql_results=$db->aql("select a, a->sequence_name from a in class gene");
foreach (@aql_results) {
    $seq_name_hash{$_->[0]}=$_->[1];
}
warn scalar keys %seq_name_hash , " gene sequence names read\n";


my %papers=();
@aql_results=$db->aql("select a, a->pmid from a in class paper");
foreach (@aql_results) {
    $papers{$_->[0]}=$_->[1];
}
warn scalar keys %papers , " papers read\n";


my %names=();

my $out;
$output = $wormbase->ontology."/gene_association.".$wormbase->get_wormbase_version_name.".wb" unless $output;
open($out, ">$output") or $log->log_and_die("cannot open $output : $!\n");

my $count=0;
my $line_count=0;

if ($gene) {   
    my $query="find gene go_term";
    my $it=$db->fetch_many(-query=>$query);
    
   
    while (my $obj=$it->next) {
	next unless $obj->isObject();
	$count++;
	if ($count % 1000 == 0) {
	    warn "$count gene objects processed\n";
	}
	
	my @lines=();
	eval {
	    @lines=split("\n", $obj->asAce);
	};
	if ($@) {
	    warn "$@\n";
	    next;
	}
	my $public_name='';
	my $species='';
	foreach (@lines) {
	    if (/Public_name\s+\"(.+)\"/) {
		$public_name=$1;
	    }
	    elsif (/Species\s+\"(.+)\"/) {
		$species=$1;
	    }
	    elsif (/GO_term\s+\"(GO:\d+)\"/) {
		my $term=$1;
		my $go_type=$db->fetch('GO_term', $term)->Type;
		s/\"//g;#"
		my @tmp=split("\t");
		
		my $ref='';
		if ($tmp[5]=~/WBPaper/) {
		    $ref="WB:$tmp[5]";
		    if ($papers{$tmp[5]}) {
			$ref.="|PMID:$papers{$tmp[5]}";
		    }
		}
		if (!$ref and $tmp[3] eq "IEA") {
		    $ref='PMID:12520011|PMID:12654719';
		}
		
		my $with='';
		if ($tmp[5]=~/WBPerson/ && $tmp[4] eq 'Person_evidence') {
		    $with="WB:$tmp[5]";
		}
		if ($tmp[4] eq 'Curator_confirmed') {
		    next;
		}
		if ($tmp[4] eq 'Inferred_automatically') {
		    $with=$tmp[5];
		}
		if ($tmp[5]=~/WBRNAi/) {     # do not parse RNAi GO terms via genes - it's done separately via phenotypes
		    next;
		}
		my $syn="";
		if ($public_name ne $seq_name_hash{$obj}) {
		    $syn=$seq_name_hash{$obj};
		}
		my $taxon="taxon:6239";
		if ($species=~/briggsae/) {
		    $taxon="taxon:6238";
		}
		my $a=$aspect{lc $go_type};
		my $type="gene";
		print $out "WB\t$obj\t$public_name\t\t$term\t$ref\t$tmp[3]\t$with\t$a\t\t$syn\t$type\t$taxon\t$date\tWB\n";
		$line_count++;
	    }
	}
	
    }
    
    warn "$count genes processed\n";
    warn "$line_count lines generated\n";
}


if ($rnai) {

    $count=0;

    my %phen2go=();
    @aql_results = $db->aql("select p, p->GO_term from p in class Phenotype where exists p->GO_term");
    foreach(@aql_results) { 
	if (defined $_->[1]) { 
	    $phen2go{$_->[0]}{$_->[1]}=1;
	}
    }
    warn scalar keys %phen2go, " phenotypes read\n";
    
    my $query="find phenotype go_term; follow rnai";
    
    my $it=$db->fetch_many(-query=>$query);
    
    
    while (my $obj=$it->next) {
	next unless $obj->isObject();
	$count++;
	if ($count % 1000 == 0) {
	    warn "$count RNAi objects processed\n";
	}
	
	my @genes_tmp=$obj->Gene;
	my @genes=();
	foreach (@genes_tmp) {
	    if ($_->right(2) eq "RNAi_primary") {
		push @genes, $_;
	    }
	    else {
	    }
	}
	my $species=$obj->Species;
	my $ref=$obj->Reference;
	my @phen_array_tmp=$obj->Phenotype;
	my @phen_array=();
	foreach (@phen_array_tmp) {
	    my $not=grep {/Not/} $_->tags();
	    if ($not) {
	    }
	    push @phen_array, $_ unless $not;
	}
	
	foreach my $gene (@genes) {
	    foreach my $phen (@phen_array) {
		if (! ($phen2go{$phen})) {
		    next;
		}
		my $taxon="taxon:6239";
		if (defined $species and $species=~/briggsae/) { #future problems!
		    $taxon="taxon:6238";
		}
		my $type="gene";
		my $public_name='';
		if ($name_hash{$gene}) {
		    $public_name=$name_hash{$gene};
		}
		my $ref_field='';
		if ($papers{$ref}) {
		    $ref_field="WB:$ref|PMID:$papers{$ref}";
		}
		else {
		    $ref_field="WB:$ref";
		}
		my $with="WB:$obj|WB:$phen";
		my $syn="";
		if ($public_name ne $seq_name_hash{$gene}) {
		    $syn=$seq_name_hash{$gene};
		}
		
		foreach my $term (keys %{$phen2go{$phen}}) {
		    my $go_type=$db->fetch('GO_term', $term)->Type;
		    my $a=$aspect{lc $go_type};
		    print $out "WB\t$gene\t$public_name\t\t$term\t$ref_field\tIMP\t$with\t$a\t\t$syn\t$type\t$taxon\t$date\tWB\n";
		    $line_count++;
		}
	    }
	}
    }
    
    
    warn "$count RNAi objects processed\n";
    warn "$line_count lines generated\n";
}

#separate species for gene associations
$wormbase->run_command("grep 'taxon:6239' $output > $output.ce", $log);
$wormbase->run_command("grep 'taxon:6238' $output > $output.cb", $log);

$log->mail;
exit();

