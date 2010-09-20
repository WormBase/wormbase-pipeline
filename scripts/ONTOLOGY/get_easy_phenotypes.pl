#!/software/bin/perl
          
#
# This was written at the suggestion of Todd for Megan Senchuk, msenchuk@mcb.harvard.edu
# who wanted a quick way of getting the observed phenotypes with RNAi evidence for all genes
#

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


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "database:s" => \$acedbpath,
	    "output:s"   => \$output
	    );

my $program_name=$0=~/([^\/]+)$/ ? $1 : '';


if ($help) {
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

$acedbpath=$wormbase->autoace unless $acedbpath;
my $tace=$wormbase->tace;

warn "connecting to database... ";
my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or $log->log_and_die("Connection failure: ". Ace->error);
warn "done\n";

my %names=();
my @aql_results=$db->aql('select a, a->sequence_name from a in class gene where a->Species = "Caenorhabditis elegans"');
foreach (@aql_results) {
  $names{$_->[0]} = $_->[1];
}
warn scalar keys %names , " genes read\n";

my %description=();
@aql_results=$db->aql('select a, a->Primary_name from a in class phenotype');
foreach (@aql_results) {
  $description{$_->[0]} = $_->[1];
}
warn scalar keys %description , " phenotypes read\n";


my $out;
$output = $wormbase->ontology."/rnai_phenotypes.".$wormbase->get_wormbase_version_name.".wb" unless $output;
open($out, ">$output") or $log->log_and_die("cannot open $output : $!\n");



my $it=$db->fetch_many(-query=>'find RNAi Phenotype');
my $count=0;

my %result;

while (my $RNAi_obj=$it->next) {
  next unless $RNAi_obj->isObject();
  $count++;
  if ($count % 1000 == 0) {
    warn "$count RNAi objects processed\n";
  }

  my @genes_tmp=$RNAi_obj->Gene;
  foreach my $this_gene (@genes_tmp) {
    if ($this_gene->right(2) eq 'RNAi_primary') {
      my @phen_array_tmp=$RNAi_obj->Phenotype;
      foreach my $phenotype (@phen_array_tmp) {
	next if (grep {/Not/} $phenotype->tags());
	my $RNAi_name = $RNAi_obj->name;
	my $RNAi_ref = $RNAi_obj->Reference;
	push @{$result{$this_gene}{$phenotype}}, "$RNAi_name|$RNAi_ref";

      }
    }
  }
}


# print the result
$count=0;
foreach my $gene (keys %names) {
  foreach my $phenotype (keys %{$result{$gene}}) {
    $count++;
    if ($count % 1000 == 0) {
      warn "$count Gene Phenotypes output\n";
    }    
    print $out "$gene\t$names{$gene}\t$description{$phenotype}\t$phenotype\t@{$result{$gene}{$phenotype}}\n";
  }
}




#while (my $obj=$it->next) {
#    next unless $obj->isObject();
#    $count++;
#    if ($count % 1000 == 0) {
#	warn "$count RNAi objects processed\n";
#    }
#    
#    my @genes_tmp=$obj->Gene;
#    my @genes=();
#    my %pheno=();
#    foreach (@genes_tmp) {
#	if ($_->right(2) eq 'RNAi_primary') {
#	    push @genes, $_;
#	}
#    }
#    my $ref=$obj->Reference;
#    my @phen_array_tmp=$obj->Phenotype;
#    foreach (@phen_array_tmp) {
#	my $not=grep {/Not/} $_->tags();
#
#	$pheno{$_}[1]='NOT' if $not;
#	$pheno{$_}[0]++;
#    }
#    foreach my $g (@genes) {
#	foreach my $p (keys %pheno) {
#	  if ($pheno{$p}[1] eq 'NOT') {next}
#	  if (!exists $names{$g}) {next} # not C. elegans
#	  print $out "$g\t$names{$g}\t$description{$p}\t$p\t$ref\t$obj\n";
#	}
#    }



close($out);


##################
# Check the files
##################


$wormbase->check_file($output, $log,
minsize => 3000000,
maxsize => 6000000,
);

$db->close;

$log->mail;
exit();

	

