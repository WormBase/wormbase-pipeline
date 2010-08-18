#!/software/bin/perl -w
          
use lib $ENV{'CVS_DIR'};
use strict;
use Ace;         
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;	

my ($ep, $ref, %genes, %at, %auth, $date);
my ($help, $debug, $test, $store, $wormbase);
my ($output, $acedbpath, $rnai, $gene, $variation, $skiplist);

$|=9;


my %opts=();
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "store:s"    => \$store,
	    "database:s" => \$acedbpath,
	    "rnai"   	 => \$rnai,
	    "gene"     	 => \$gene,
	    "variation"  => \$variation,
	    "output:s"   => \$output,
	    "skiplist:s" => \$skiplist, 
	    );

my $program_name=$0=~/([^\/]+)$/ ? $1 : '';


if ($help) {
    print "usage: $program_name [options] -output output -database database\n";
    print "       -help            help - print this message\n";
    print "       -output <output>   output file\n";
    print "       -database <database> path to database\n";
    print "       -rnai            generate RNAi2GO mapping; default no\n";
    print "       -gene            generate gene mapping; default no\n";
    print "       -variation       generate variation2GO mapping; default no\n";
    print "       -skiplist <papers_to_skip> list of papers to skip when generating RNAi2GO mappings\n";
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

# check that the GO_term objects all have a Type tag set
&check_go_term();

warn "connecting to database... $acedbpath";
my $db = Ace->connect(-path => $acedbpath,  -program => $wormbase->tace) or $log->log_and_die("Connection failure: ". Ace->error);
warn "... done\n";

my %species_taxon_hash=("Caenorhabditis elegans"  => 6239,
			"Caenorhabditis briggsae" => 6238,
			"Caenorhabditis remanei"  => 31234,
			"Caenorhabditis brenneri" => 135651,
			"Caenorhabditis vulgaris" => 31233,
			"Caenorhabditis trinidad" => "N/A",
			"Brugia pahangi"          => 6280,
			"Onchocerca volvulus"     => 6282,
			"Pristionchus pacificus"  => 54126,
			'Caenorhabditis japonica' => 281687,
			);

my %name_hash=();
my @aql_results=$db->aql('select a, a->public_name from a in class gene');
foreach (@aql_results) {
    $name_hash{$_->[0]}=$_->[1];
}
warn scalar keys %name_hash , " gene public names read\n";


my %seq_name_hash=();
@aql_results=$db->aql('select a, a->sequence_name from a in class gene');
foreach (@aql_results) {
    $seq_name_hash{$_->[0]}=$_->[1];
}
warn scalar keys %seq_name_hash , " gene sequence names read\n";


my %paper_fields=();
my %paper_accno=();
@aql_results=$db->aql('select a, a->Database[2], a->Database[3]  from a in class paper');
foreach (@aql_results) {
  $paper_fields{$_->[0]} = $_->[1]; # database field e.g. 'PMID'
  $paper_accno{$_->[0]} = $_->[2];  # database accession number e.g. '12393910'
}
warn scalar keys %paper_fields , " papers read\n";

my %papers_to_skip=();
if ($skiplist) {
    open (SKIP, "<$skiplist") or $log->log_and_die("cannot open $skiplist : $!\n");
    while (<SKIP>) {
	chomp;
	s/\s+//g;
	next unless /WBPaper/;
	$papers_to_skip{$_}=1;
    }
    close (SKIP);
    warn scalar keys %papers_to_skip, " papers will be skipped\n";
}

my %names=();

my $out;
$output = $wormbase->ontology."/gene_association.".$wormbase->get_wormbase_version_name.".wb" unless $output;
open($out, ">$output") or $log->log_and_die("cannot open $output : $!\n");

my $count=0;
my $line_count=0;

if ($gene) {   
    my $it=$db->fetch_many(-query=>'find gene go_term');
    
   
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
		    if ($paper_fields{$tmp[5]}) {
			$ref.="|$paper_fields{$tmp[5]}:$paper_accno{$tmp[5]}"; # type of database field and the accession_number e.g. 'PMID:12393910'
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
		if ($tmp[5]=~/WBPhenotype/) {     # do not parse phenotype-based GO terms via genes - it's done separately via phenotypes
		    next;
		}
		my $syn="";
		if ($public_name ne $seq_name_hash{$obj}) {
		    $syn=$seq_name_hash{$obj};
		}
		my $taxon="taxon:$species_taxon_hash{$species}";
		if (!exists $species_taxon_hash{$species}) {print "ERROR: missing species = $species\n";}
		my $a=$aspect{lc $go_type};
		my $type="gene";
		print $out "WB\t$obj\t$public_name\t\t$term\t$ref\t$tmp[3]\t$with\t$a\t\t$syn\t$type\t$taxon\t$date\tWB\t\t\n";
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
    
    my $it=$db->fetch_many(-query=>'find phenotype go_term; follow rnai');
    
    while (my $obj=$it->next) {
	next unless $obj->isObject();
	$count++;
	if ($count % 1000 == 0) {
	    warn "$count RNAi objects processed\n";
	}
	
	my @genes_tmp=$obj->Gene;
	my @genes=();
	foreach (@genes_tmp) {
	    if ($_->right(2) eq 'RNAi_primary') {
		push @genes, $_;
	    }
	}
	my $ref=$obj->Reference;
	next if $papers_to_skip{$ref};
	my @phen_array_tmp=$obj->Phenotype;
	my @phen_array=();
	foreach (@phen_array_tmp) {
	    my $not=grep {/Not/} $_->tags();
	    push @phen_array, $_ unless $not;
	}
	
	foreach my $gene (@genes) {
	    my $species=$gene->Species;
	    foreach my $phen (@phen_array) {
		if (! ($phen2go{$phen})) {
		    next;
		}
		my $taxon="taxon:$species_taxon_hash{$species}";
		my $type="gene";
		my $public_name='';
		if ($name_hash{$gene}) {
		    $public_name=$name_hash{$gene};
		}
		my $ref_field='';
		if ($paper_fields{$ref}) {
		    $ref_field="WB:$ref|$paper_fields{$ref}:$paper_accno{$ref}"; # type of database field and the accession_number e.g. 'PMID:12393910'";
		}
		elsif ($ref) {
		    $ref_field="WB:$ref";
		} else {
		  # if the RNAi object does not have a Reference, give a warning and skip this RNAi because the $ref_field is a mandatory term
		  warn "No reference for $obj\n";
		  next;
		}
		my $with="WB:$obj|WB:$phen";
		my $syn="";
		if ($public_name ne $seq_name_hash{$gene}) {
		    $syn=$seq_name_hash{$gene};
		}
		
		foreach my $term (keys %{$phen2go{$phen}}) {
		    my $go_type=$db->fetch('GO_term', $term)->Type;
		    my $a=$aspect{lc $go_type};
		    print $out "WB\t$gene\t$public_name\t\t$term\t$ref_field\tIMP\t$with\t$a\t\t$syn\t$type\t$taxon\t$date\tWB\t\t\n";
		    $line_count++;
		}
	    }
	}
    }
    
    
    warn "$count RNAi objects processed\n";
    warn "$line_count lines generated\n";
}

if ($variation) {

    $count=0;

    my %phen2go=();
    @aql_results = $db->aql("select p, p->GO_term from p in class Phenotype where exists p->GO_term");
    foreach(@aql_results) { 
	if (defined $_->[1]) { 
	    $phen2go{$_->[0]}{$_->[1]}=1;
	}
    }
    warn scalar keys %phen2go, " phenotypes read\n";
    
    my $it=$db->fetch_many(-query=>'find phenotype go_term; follow variation');
    
    while (my $obj=$it->next) {
	next unless $obj->isObject();
	$count++;
	if ($count % 1000 == 0) {
	    warn "$count variation objects processed\n";
	}

	if ($obj->Gain_of_function) { # skip gain-of-function alleles
	    next;
	}

	my @genes=$obj->Gene;
	my @phen_array_tmp=$obj->Phenotype;
	my %phen_hash=();
	foreach (@phen_array_tmp) {
	    if (! ($phen2go{$_})) {
		next;
	    }
	    my $not=grep {/Not/} $_->tags();
	    if ($not) {
		next;
	    }

	    $phen_hash{$_}{count}++;
	    if ($_->at("Paper_evidence")) {
		$phen_hash{$_}{Paper_evidence}=$_->at("Paper_evidence[1]");
	    }
	    if ($_->at("Person_evidence")) {
		$phen_hash{$_}{Person_evidence}=$_->at("Person_evidence[1]")=~/(WBPerson\d+)/ ? $1 : '';
	    }
	    if ($_->at("Curator_confirmed")) {
		$phen_hash{$_}{Curator_confirmed}=$_->at("Curator_confirmed[1]")=~/(WBPerson\d+)/ ? $1 : '';
	    }

	}

	    
	
	foreach my $gene (@genes) {
	    my $species=$gene->Species;
	    foreach my $phen (sort {$a cmp $b} keys %phen_hash) {
		my $taxon="taxon:$species_taxon_hash{$species}";
		my $type="gene";
		my $public_name='';
		if ($name_hash{$gene}) {
		    $public_name=$name_hash{$gene};
		}
		my $ref_field='';
		if ($phen_hash{$phen}{Paper_evidence}) {
		    if ($paper_fields{$phen_hash{$phen}{Paper_evidence}}) {
			$ref_field="WB:$phen_hash{$phen}{Paper_evidence}|$paper_fields{$phen_hash{$phen}{Paper_evidence}}:$paper_accno{$phen_hash{$phen}{Paper_evidence}}"; # type of database field and the accession_number e.g. 'PMID:12393910'";
		    }
		    else {
			$ref_field="WB:$phen_hash{$phen}{Paper_evidence}";
		    }
		}
		elsif ($phen_hash{$phen}{Person_evidence}) {
		    $ref_field="WB:$phen_hash{$phen}{Person_evidence}";
		}
		elsif ($phen_hash{$phen}{Curator_confirmed}) {
		    $ref_field="WB:$phen_hash{$phen}{Curator_confirmed}";
		}
		my $with="WB:$obj|WB:$phen";
		my $syn="";
		if ($public_name ne $seq_name_hash{$gene}) {
		    $syn=$seq_name_hash{$gene};
		}
		
		foreach my $term (keys %{$phen2go{$phen}}) {
		    my $go_type=$db->fetch('GO_term', $term)->Type;
		    my $a=$aspect{lc $go_type};
		    print $out "WB\t$gene\t$public_name\t\t$term\t$ref_field\tIMP\t$with\t$a\t\t$syn\t$type\t$taxon\t$date\tWB\t\t\n";
		    $line_count++;
		}
	    }
	}
    }
    
    
    warn "$count variation objects processed\n";
    warn "$line_count lines generated\n";
}


close($out);

#separate species for gene associations
$wormbase->run_command("grep 'taxon:6239' $output > $output.ce", $log);
$wormbase->run_command("grep 'taxon:6238' $output > $output.cb", $log);
$wormbase->run_command("grep 'taxon:31234' $output > $output.rem", $log);
$wormbase->run_command("grep 'taxon:54126' $output > $output.ppa", $log);
$wormbase->run_command("grep 'taxon:281687' $output > $output.cjp", $log);
$wormbase->run_command("grep 'taxon:135651' $output > $output.cbn", $log);


##################
# Check the files
##################


$wormbase->check_file($output, $log,
minsize => 11000000,
maxsize => 70000000,
lines => ['^WB\tWBGene\d+\t\S+\t\tGO\:\d+'],
);

$wormbase->check_file("$output.ce", $log,
minsize => 5000000,
maxsize => 16000000,
lines => ['^WB\tWBGene\d+\t\S+\t\tGO\:\d+'],
);

$wormbase->check_file("$output.cb", $log,
minsize => 1500000,
maxsize => 8000000,
lines => ['^WB\tWBGene\d+\t\S+\t\tGO\:\d+'],
);

$wormbase->check_file("$output.cjp", $log,
minsize =>  1600000,
maxsize =>  8000000,
lines => ['^WB\tWBGene\d+\t\S+\t\tGO\:\d+'],
);

$wormbase->check_file("$output.cbn", $log,
minsize =>  7000000,
maxsize =>  9000000,
lines => ['^WB\tWBGene\d+\t\S+\t\tGO\:\d+'],
);


$wormbase->check_file("$output.ppa", $log,
minsize =>  1700000,
maxsize =>  8000000,
lines => ['^WB\tWBGene\d+\t\S+\t\tGO\:\d+'],
);

$wormbase->check_file("$output.rem", $log,
minsize =>  8000000,
maxsize =>  9000000,
lines => ['^WB\tWBGene\d+\t\S+\t\tGO\:\d+'],
);


$db->close;

$log->mail;
exit();



#####################################################################################




# Even though the GO term definitions are imported to autoace from
# citace, there may be some GO terms that are referenced by InterPro
# etc. that have not yet been defined by Caltech. We therefore do a
# check for those GO_term object that do not have a fully populated
# set of tags and we populate them here.

sub check_go_term {

  # get the set of GO_term objects missing a Type tag

  my $tace            = $wormbase->tace;        # TACE PATH
  my $cmd = "Query Find GO_term Where !Type\nlist -a\nquit";

  my @go;
  my %go;
  open (TACE, "echo '$cmd' | $tace $acedbpath |");
  while (<TACE>) {
    chomp;
    next if (/acedb\>/);
    next if (/\/\//);
    if (/GO_term\s+:\s+\"(\S+)\"/) {
      push @go, $1;
      $go{$1} = 1;
    }
  }
  close TACE;
  
  if (@go == 0) {
    $log->write_to("No GO_terms are missing a Type tag.\nAll appear OK.\n");
  } else {
    
    $log->write_to("There are " . scalar @go . " GO_terms that are incomplete - fixing them\n");
    
    $log->write_to("GO_terms objects are:\n@go\n\n");
    
    my $gocount = 0;
    # get the full GO file - this is in OBO format
    # http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo
    
    #[Term]
    #id: GO:0000003
    #name: reproduction
    #namespace: biological_process
    #alt_id: GO:0019952
    #alt_id: GO:0050876
    #def: "The production by an organism of new individuals that contain some portion of their genetic material inherited from that organism." [GOC:go_curators, GOC:isa_complete, ISBN:0198506732 "Oxford Dictionary of Biochemistry and Molecular Biology"]
    #subset: goslim_generic
    #subset: goslim_pir
    #subset: goslim_plant
    #subset: gosubset_prok
    #synonym: "reproductive physiological process" EXACT []
    #xref: Wikipedia:Reproduction
    #is_a: GO:0008150 ! biological_process
    
    # We only really want the 'namespace' field for the 'Type' tag, but
    # grap a few other things as well
    
    
    my $go_obo_file = "/tmp/GO_file_.".$wormbase->species;
    $log->write_to("Attempting to wget the latest GO data\n");
    `wget -q -O $go_obo_file http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo` and die "$0 Couldnt get goslim_generic.obo\n";
    
    
    $log->write_to("Opening file $go_obo_file\n");
    open (GOIN,"<$go_obo_file") or die "cant open $go_obo_file\n";
    
    
    my $acefile = "$acedbpath/acefiles/go_defs.ace";
    
    open (GOOUT,">$acefile") or die "cant write to $acefile\n";
    
    
    print "\treading data . . . \n";
    
    my $id="";
    my $namespace="";
    my $name="";
    my $def="";
    
    while (<GOIN>){
      chomp;
      if ($_ =~ /^id:\s+(\S+)/) {$id = $1}
      if ($_ =~ /name:\s+(\S+)/) {$name = $1}
      if ($_ =~ /namespace:\s+(\S+)/) {$namespace = $1}
      if ($_ =~ /def:\s+\"(.+)\"/) {$def = $1}
      if ($_ =~ /^\s*$/ && exists $go{$id}) {
	print GOOUT "\nGO_term : \"$id\"\n";
	print GOOUT "Definition \"$def\"\n";
	print GOOUT "Term \"$name\"\n";
	print GOOUT "Type \"".ucfirst($namespace)."\"\n";
	print GOOUT "\n";
	$gocount++;
      }
    }
    
    $log->write_to("added $gocount GO definitions\n");
    
    close GOIN;
    close GOOUT;
    
    # load file to autoace if -load specified
    $wormbase->load_to_database($wormbase->autoace, "$acefile", 'go_defs', $log);
    
    # tidy up and exit
    $wormbase->run_command("rm $go_obo_file",$log);
  }
}





