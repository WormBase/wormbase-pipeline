#!/software/bin/perl -w

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

######################################
# variables and command-line options #
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $user,);
my $database;

my $acefile = "/nfs/wormpub/DATABASES/geneace/CGC/orthos.ace";
my $batchfile  = "/nfs/wormpub/DATABASES/geneace/CGC/batch_load";
my $genelist;

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
            "verbose"    => \$verbose,
            "store:s"    => \$store,
	    "database:s" => \$database,
	    "batch:s"    => \$batchfile,
	    "ace:s"      => \$acefile,
	    "list:s"     => \$genelist,
	    "user:s"     => \$user, #manditory option as required for .ace output.
            );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} 
else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			   );
}

my $log = Log_files->make_build_log($wormbase);
$log->log_and_die("-list is compulsory. This is the list of elegans genes to transfer names from!\n") unless ($genelist);
if ($debug) {
  print "-list is compulsory. This is the list of elegans genes to transfer names from!\n" unless ($genelist);
}

my %CGC_species = ('briggsae' => 'Cbr',
		   'remanei'  => 'Cre',
		   'japonica' => 'Cjp',
		   'brenneri' => 'Cbn',
		   'pacificus'=> 'Ppa',
		   'brugia'   => 'Bma',
                   'ovolvulus'=> 'Ovo',
                   'sratti'   => 'Sra',
		  );

$database = $database or $wormbase->database('geneace');
$log->write_to("Database : $database\n\n");

my $ace;
my $namedb;
open($ace, ">$acefile")      or $log->log_and_die("cant write acefile - $acefile : $!\n");
open($namedb, ">$batchfile") or $log->log_and_die("cant write batch load file - $batchfile: $!\n");
open(GENES,"<$genelist")     or $log->log_and_die("cant read gene list - $genelist : $!\n");

my %person = ('mt3' => 'WBPerson2970',
	      'pad' => 'WBPerson1983',
	      'klh' => 'WBPerson3111',
	      'gw3' => 'WBPerson4025',
	      'mh6' => 'WBPerson4055',
	     );

unless (defined $user) {
  print "\nERROR: please specify a user ID - e.g. mt3\n\n[INPUT]:";
  my $tmp = <STDIN>;
  chomp $tmp;
  if (($tmp ne 'mt3') && ($tmp ne 'pad') && ($tmp ne 'ar2') && ($tmp ne 'gw3') && ($tmp ne 'mh6')){
    $log->log_and_die("UNKNOWN USER.....TERMINATE\n\n");
    print "UNKNOWN USER.....TERMINATE\n\n" if $debug;
  }
  else {
    $user = "$tmp";
  }
}


##################################################################
# Main Body
##################################################################

my $acedb = Ace->connect('-path' => $database) or Ace->error;


while(<GENES>){
    chomp;
    my $gene = $_;
    print "Processing:$gene\n";
    unless ($gene =~ /WBGene\d{8}$/) { warn "$_ bad gene format\n";next; }
    

    my $geneObj = $acedb->fetch('Gene',$gene);
    unless ($geneObj->Species->name  eq $wormbase->full_name) { warn "$gene is not a ".$wormbase->full_name('-gspecies' => 1). " gene\n";next;}
    if (!defined $geneObj) { warn "Could not retrieve $gene data, please check your input data is correct\n";next;}

    my $cgc;
    if($geneObj->CGC_name) {
	$cgc = $geneObj->CGC_name->name;
	my @orthos = $geneObj->Ortholog;
	my %store;
	foreach my $orth (@orthos){
	    my $id = $orth->name;
	    my ($spe) = $orth->right->name =~ /(\w+)$/;
	    push(@{$store{$spe}}, $id);
	}

	foreach my $species(keys %store) {
	    if( scalar @{$store{$species}} > 1){
		print STDERR "paralogs : cel $gene $species: ".@{$store{$species}}."\n";
	    }
	    else {
		print ${$store{$species}}[0]." is $cgc\n";
		my $new_name = $CGC_species{$species}."-$cgc";
		&write_new_orthology($gene, ${$store{$species}}[0], $new_name); #gene ortho cgc
	    }
	}
    }
    else {
	print "$gene has no CGC_name\n";
    }
}
close $ace;
close $namedb;
close GENES;
$log->write_to("\nCreated orthos.ace and batch_load under /nfs/wormpub/DATABASES/geneace/CGC/\n\n
Dont forget to load the ace file in to Geneace.(default is orthos.ace)\n\nFinished\n");
$log->mail;
print "\nCreated orthos.ace and batch_load under /nfs/wormpub/DATABASES/geneace/CGC/\n\n
Dont forget to load the ace file in to Geneace.(default is orthos.ace)\n\nFinished\n" if $debug;
exit;

sub write_new_orthology {
    my $gene = shift;
    my $ortholog = shift;
    my $new_name = shift;
    my $geneObj = $acedb->fetch('Gene',$ortholog);
    my $version = $geneObj->Version->name;
    $version++;

    print $ace "\nGene : $ortholog\nVersion $version\nCGC_name $new_name From_analysis Inferred_from_orthology\n";
    print $ace "Public_name $new_name\nVersion_change $version now $person{$user} Name_change CGC_name $new_name\n";
    my ($class) = $new_name =~ /-(\w+)-/;
    print $ace "Gene_class $class\n";
    $log->write_to("Transfering CGC_name: $new_name from $gene to $ortholog");

    if($geneObj->CGC_name) {
	#get existing evidence
	my $evidence = $geneObj->CGC_name->right(2);
	my $old_name = $geneObj->CGC_name->name;
	my $old_class = $geneObj->Gene_class->name;
	print $ace "Version_change ",$version," now $person{$user} Name_change Other_name $old_name\n";
	print $ace "Other_name $old_name\n";
	#print old name as Other_name with evidence transferred.
	foreach ($geneObj->CGC_name(2)){
	    print $ace "Other_name $old_name ", $_->name."\t".$_->right->name."\n";
	}

	#gene class update
	print $ace "\nGene_class : $old_class\nOld_member $old_name\n";
	$log->write_to(" and replacing $old_name");
    }

    print $namedb "$ortholog\t$new_name\n";
    $log->write_to("\n");
}
