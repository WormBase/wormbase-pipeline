#!/usr/local/bin/perl5.6.0 -w
# 
# current_DB_check.pl
#
# by Keith Bradnam
#
# Script to run consistency checks on the current_DB database
# to look for bogus sequence entries
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2002-07-18 15:52:22 $

use Ace;
use lib "/wormsrv2/scripts/"; 
use Wormbase;
use strict;
use Getopt::Std;
 
##############################
# command-line options       #
##############################
our $opt_d = "";      # Help/Usage page
our $opt_v = "";
getopts ('d:v');

# set which database to check
my $database;

if ($opt_d){
  $database = $opt_d;
}
else{
  $database = "/wormsrv2/current_DB/";
}


our $log;
&create_log_files($database);


# toggle verbose mode
my $verbose;
if ($opt_v){
  $verbose = "ON";
}

############################################
# Initialise variables
############################################

my $errors;
my %problems; # store problems in a double hash, first key being timestamp name
my @other; # store uncategorised problems

print "Checking $database\n\n";

our $tace = glob("~wormpub/ACEDB/bin.ALPHA_4/tace");   # tace executable path
my $db = Ace->connect(-path  => "$database",
		      -program =>$tace) || do { print LOG "Connection failure: ",Ace->error; die();};

my @genome_seqs = $db->fetch(-class => 'Genome_sequence',
		      -name  => '*');

print "\nLooking for spurious sequences\n";

foreach my $seq (@genome_seqs){
  my @subseqs = $db->fetch(-class=>'Sequence',-name=>"$seq.*");
  foreach my $subseq (@subseqs){
    if(!defined($subseq->at('Structure.From.Source'))){  
      $errors++;     
      my $category = 0;
      
      if(defined($subseq->at('DB_info.Database'))){
	my $tag = "Database";
	my $timestamp = &aql_query($subseq, $tag);
	$category = 1;
      }
      if(defined($subseq->at('Visible.RNAi_result'))){
	my $tag = "RNAi_result";
	my $tag_pair = "Predicted_gene";
	my $pair_class = "RNAi";
	my $timestamp = &aql_query($subseq, $tag, $tag_pair, $pair_class);
	$category = 1;
      }
      if(defined($subseq->at('Visible.Locus_genomic_seq'))){
	my $tag = "Locus_genomic_seq";
	my $tag_pair = "Genomic_sequence";
	my $pair_class = "Locus";
	my $timestamp = &aql_query($subseq, $tag, $tag_pair, $pair_class);
	$category = 1;
      }
      if(defined($subseq->at('Visible.Contained_in_operon'))){
	my $tag = "Contained_in_operon";
	my $tag_pair = "Contains_CDS";
	my $pair_class = "Operon";
	my $timestamp = &aql_query($subseq, $tag, $tag_pair, $pair_class);
	$category = 1;
      }
      if(defined($subseq->at('Allele'))){
	my $tag = "Allele";
	my $tag_pair = "Sequence";
	my $pair_class = "Allele";
	my $timestamp = &aql_query($subseq, $tag, $tag_pair, $pair_class);
	$category = 1;
      }
      if(defined($subseq->at('Visible.Matching_cDNA'))){
	my $tag = "Matching_cDNA";
	my $tag_pair ="Matching_Genomic";
	my $pair_class = "Sequence";
	my $timestamp = &aql_query($subseq, $tag, $tag_pair, $pair_class);
	$category = 1;
      }
      if(defined($subseq->at('Visible.Reference'))){
	my $tag = "Reference";
	my $tag_pair ="Sequence";
	my $pair_class = "Paper";
	my $timestamp = &aql_query($subseq, $tag, $tag_pair, $pair_class);
	$category = 1;
      }
      if(defined($subseq->at('Visible.Paired_read'))){
	my $tag = "Paired_read";
	my $tag_pair = "Paired_read";
	my $pair_class = "Sequence";
	my $timestamp = &aql_query($subseq, $tag, $tag_pair, $pair_class);
	$category = 1;
      }
      if(defined($subseq->at('DNA'))){
	my $tag = "DNA";
	my $timestamp = &aql_query($subseq, $tag);
	$category = 1;
      }
      if(defined($subseq->at('Visible.Remark'))){
	my $tag = "Remark";
	my $timestamp = &aql_query($subseq, $tag);
	$category = 1;
      }
      if(defined($subseq->at('Visible.Expr_pattern'))){
	my $tag = "Expr_pattern";
	my $tag_pair = "Sequence";
	my $pair_class = "Expr_pattern";
	my $timestamp = &aql_query($subseq, $tag, $tag_pair, $pair_class);
	$category = 1;
      }
      if(defined($subseq->at('Visible.Provisional_description'))){
	my $tag = "Provisional_description";
	my $timestamp = &aql_query($subseq, $tag);
	$category = 1;
      }
      if(defined($subseq->at('Visible.Detailed_description'))){
	my $tag = "Detailed_description";
	my $timestamp = &aql_query($subseq, $tag);
	$category = 1;
      }
      if(defined($subseq->at('Visible.Concise_description'))){
	my $tag = "Concise_description";
	my $timestamp = &aql_query($subseq, $tag);
	$category = 1;
      }
      if(defined($subseq->at('Visible.Corresponding_protein'))){
	my $tag = "Corresponding_protein";
	my $tag_pair = "Corresponding_DNA";
	my $pair_class = "Protein";
	my $timestamp = &aql_query($subseq, $tag, $tag_pair, $pair_class);
	$category = 1;
      }
      if(defined($subseq->at('Visible.GO_term'))){
	my $tag = "GO_term";
	my $tag_pair = "Sequence";
	my $pair_class = "GO_term";
	my $timestamp = &aql_query($subseq, $tag, $tag_pair, $pair_class);
	$category = 1;
      }
      if(defined($subseq->at('Properties'))){
	my $tag = "Properties";
	my $timestamp = &aql_query($subseq, $tag);
	$category = 1;
      }
      if ($category == 0){
	push(@other,$subseq);
	print "$subseq - Other problem\n" if $verbose;
      }
    }    
  }
}

print "\n$errors errors found\n" if $verbose;
print LOG "\n$errors errors found\n\n";

foreach my $i (@other){
  print LOG "Other - $i\n";
}

foreach my $j (keys(%problems)){
  foreach my $k (keys(%{$problems{$j}})){
    print LOG "$j $k @{${$problems{$j}}{$k}}\n";
  }
}


##########################################
# Tidy up and mail relevant log files
##########################################

$db->close;
print LOG "\n$0 ended at ",`date`,"\n";
close(LOG);

my $maintainer = "wormbase-dev\@wormbase.org";
&mail_maintainer($0,$maintainer,$log);
exit(0);

################################################


sub create_log_files{
  my $database = shift;
  my $rundate    = `date +%y%m%d`; chomp $rundate;
  $log = "/wormsrv2/logs/current_DB_check.log.$rundate.$$";
  open(LOG,">$log") || die "cant open $log";

  print LOG "$0 started at ",`date`,"\n";
  print LOG "Checking $database\n";
  print LOG "=============================================\n";


}

################################################
# Grab timestamps from each object using aql
################################################
sub aql_query{
  my $subseq = shift;
  my $tag = shift;
  my $value = ""; 
  ($value) = $subseq->$tag;
  my $tag_pair = shift;
  my $pair_class = shift;

  my $aql_query = "select s,s->$tag.node_session from s in object(\"Sequence\",\"$subseq\")";
#  print "$aql_query\n";
  my @aql = $db->aql($aql_query);
  my $source = $aql[0]->[1];
  $source =~ s/\d{4}\-\d{2}\-\d{2}_\d{2}:\d{2}:\d{2}_//;

  # check for splice variants
  my $splice_check ="";
  my $splice_variant = $db->fetch(-class => 'Sequence',-name  => "${subseq}a");
  $splice_check =  " - splice variants exist for this sequence" if(defined($splice_variant));

  
  
  # if timestamp is wormpub, need to look at the corresponding object to see if that has 
  # more informative timestamp
  if ($source =~ m/wormpub/){
    if(defined($tag_pair) && defined($pair_class)){
      my $aql_query = "select s,s->${tag_pair}.node_session from s in object(\"$pair_class\",\"$value\")";
      my @aql = $db->aql($aql_query);
      my $source = $aql[0]->[1];
      $source =~ s/\d{4}\-\d{2}\-\d{2}_\d{2}:\d{2}:\d{2}_//;
      $problems{$source}{$subseq} = [$tag,$splice_check];
    }
  }
  else{
    $problems{$source}{$subseq} = [$tag,$splice_check];
#    print "@{$problems{$source}{$subseq}}\n";
  }

  print "$subseq \t$tag \t$source\n" if $verbose;   

  return($source);
}
