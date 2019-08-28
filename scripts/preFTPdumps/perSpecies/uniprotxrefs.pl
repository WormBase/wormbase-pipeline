#!/usr/bin/env perl
# generate a uniprot xref report on a per gene paper.
# lines are:
# gene paper pmid type
# WBGene1 WBPaper1 pmed123 Variation;RNAi
# types are: GO_annotation Interaction RNAi Variation Expr_pattern Gene Variation

use Ace;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use IO::File;
use strict;

my ($debug,$test,$store,$database,$outfile,$species);
GetOptions (
            'debug=s'    => \$debug,
	    'test'       => \$test,
	    'store=s'    => \$store,
	    'database=s' => \$database,
	    'outfile=s'  => \$outfile,
	    'species=s'  => \$species,
	  )||die($!);

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug => $debug, -test => $test, -organism => $species);
}

$outfile||=$wormbase->reports.'/uniprot_xrefs.txt';

my $dbpath = $database||$wormbase->autoace;
my $db = Ace->connect(-path => $dbpath);

my $log = Log_files->make_build_log($wormbase);

# will be $g{WBGeneID}->{WBPaperId}->{Type}=1
my %geneHash;

my $t = localtime;
$log->write_to("$t - Processing GO Annotations\n") if $debug;
# GO_annotations
my $it = $db->fetch_many(-query => 'find GO_annotation;Reference');
while (my $go = $it->next){
   my @papers = $go->Reference;
   my @genes  = $go->Gene;
   foreach my $g (@genes){
	   next unless $g->Species eq $wormbase->long_name;
	   map {$geneHash{"$g"}->{"$_"}->{GO}=1} @papers;
   }
}

my $t = localtime;
$log->write_to("$t - Processing Interactions\n") if $debug;
# Interaction
my $it = $db->fetch_many(-query => 'find Interaction;Physical;Paper');
while (my $interaction = $it->next){
   my @papers = $interaction->Paper;
   my @genes  = $interaction->Interactor_overlapping_gene;
   foreach my $g (@genes){
	   next unless $g->Species eq $wormbase->long_name;
	   map {$geneHash{"$g"}->{"$_"}->{PPI}=1} @papers;
   }
}

my $t = localtime;
$log->write_to("$t - Processing RNAi\n") if $debug;
# RNAi
my $it = $db->fetch_many(-query => 'find RNAi;Phenotype;Reference');
while (my $rnai = $it->next){
   my @papers = $rnai->Reference;
   my @genes  = $rnai->Gene;
   foreach my $g (@genes){
	   next unless $g->Species eq $wormbase->long_name;
	   map {$geneHash{"$g"}->{"$_"}->{Phenotype}=1} @papers;
   }
}

my $t = localtime;
$log->write_to("$t - Processing Variations\n") if $debug;
# Variations - has a hack to prevent giface leakage
my @vars = grep {/WBVar\d+/} map {"$_"} $db->fetch(-query => 'find Variation;Gene;Reference');
while ( my $v = shift @vars){
   my $var = $db->fetch(Variation => $v);
   next unless $var->Species eq $wormbase->long_name;
   my @papers = $var->Reference;
   my @genes  = $var->Gene;
   foreach my $g (@genes){
	   map {$geneHash{"$g"}->{"$_"}->{Phenotype}=1} @papers;

           if ($var->Sequence && $var->Evidence && $var->Type_of_mutation){
	     map {$geneHash{"$g"}->{"$_"}->{Sequence}=1} grep {$_->class eq 'Paper'} $var->get('Type_of_mutation',5);
           }
   }
}

my $t = localtime;
$log->write_to("$t - Processing Expression Patterns\n") if $debug;
# Expression pattern
my $it = $db->fetch_many(-query => 'find Expr_pattern;Gene;Reference');
while (my $ep = $it->next){
   my @papers = $ep->Reference;
   my @genes  = $ep->Gene;
   foreach my $g (@genes){
	   next unless $g->Species eq $wormbase->long_name;
	   map {$geneHash{"$g"}->{"$_"}->{Expression}=1} @papers;
   }
}

my $t = localtime;
$log->write_to("$t - Processing Genes\n") if $debug;
# Gene / Diseases
my $it = $db->fetch_many( -query => 
	'find Gene WBGene*;Species="'.$wormbase->long_name.'";Disease_relevance OR Experimental_model'
);
while (my $g = $it->next){
   my @papers;
   if ($g->Disease_relevance){
	   push @papers, grep {$_->class eq 'Paper'} $g->get('Disease_relevance',4);
   }
   if ($g->Experimental_model){
	   push @papers, grep {$_->class eq 'Paper'} $g->get('Experimental_model',4);
   }
   map {$geneHash{"$g"}->{"$_"}->{Disease}=1} @papers;
}

# print the whole thing
my $fh = IO::File->new(">$outfile");
if (defined $fh) {
  while (my($geneid,$v)=each %geneHash){
	while(my ($wbpaper,$types)= each %$v){
		print $fh join("\t",$geneid,$wbpaper,get_pmid($wbpaper),join(';',keys %$types),"\n")
	}
  }
  $fh->close
}else{
	$log->write_to("ERROR: cannot open $outfile for writing\n");
}

$log->mail();

###########################
sub get_pmid{
	my ($wbpaperId)=@_;
	my $paper = $db->fetch(Paper => $wbpaperId);
	my $pmid = 'nopmid';
	foreach my $pref ($paper->Database) {
         if ($pref->name eq 'MEDLINE') {
          $pmid = $pref->right->right->name;
         }
        }
	return $pmid;
}
