#!/usr/local/bin/perl5.6.1 -w
# 
# current_DB_check.pl
#
# by Keith Bradnam
#
# Script to run consistency checks on the current_DB database
# to look for bogus sequence entries
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2003-03-21 15:48:37 $

use strict;
use lib "/wormsrv2/scripts/"; 
use Wormbase;
use Ace;
use Getopt::Std;
 
##############################
# command-line options       #
##############################
our $opt_d = "";      # specify database to run against
our $opt_v = "";      # verbose mode
our $opt_h = "";      # display help
our $opt_t = "";      # test mode (run script but don't email anyone)
getopts ('d:vht');

&usage if ($opt_h);

# toggle verbose mode (turns on reporting of errors, but doesn't email people)
our $verbose;
$verbose = "ON" if ($opt_v);
my $test;
$test = "ON" if ($opt_t);

########################################################
# Set database details, default is /wormsrv2/current_DB
########################################################
our $database = "/wormsrv2/current_DB/";
$database = $opt_d if ($opt_d);
our $tace = &tace;   # tace executable path
our $db = Ace->connect(-path  => "$database",
		      -program =>$tace) || do { print ALL_LOG "Connection failure: ",Ace->error; die();};

print "Checking $database\n\n";



##############################################
# Initialise variables, set up log files
###############################################
my %problems; # store problems in a double hash, first key being timestamp name
my @other; # store uncategorised problems

# make one global log file but also split for different groups
our $all_log;
our $sanger_log;
our $cshl_log;
our $stlouis_log;
our $caltech_log;
our $WS_version  = &get_wormbase_version_name;

&create_log_files($database);


###################################
# Report errors to log file
###################################

# Count all the problems!
my $sanger_counter = 0;
my $caltech_counter = 0;
my $cshl_counter = 0;
my $stlouis_counter = 0;
my $all_counter = 0;

############################################################
# Check ?Sequence class (including predicted genes)
############################################################

print "\nLooking for spurious sequences\n";


# Checks sequences connected to multiple loci
&find_multiple_loci;

# Search everything else, one sequence at a time
&process_sequences;


# Count problems and print output

foreach my $i (@other){
  print ALL_LOG "Other - $i\n";
  print SANGER_LOG "Other - $i\n";
  print "Other - $i\n" if $verbose;
  $all_counter++;
  $sanger_counter++;
}

  

foreach my $j (keys(%problems)){
  foreach my $k (keys(%{$problems{$j}})){
    print ALL_LOG "$j $k @{${$problems{$j}}{$k}}\n";
    if ($j ne "caltech" && $j ne "csh" && $j ne "stlace" && $j ne "brigace"){
      print SANGER_LOG "$j $k @{${$problems{$j}}{$k}}\n";      
      $sanger_counter++;
      $all_counter++;
    }
    
    if ($j eq "csh"){
      print CSHL_LOG "$k @{${$problems{$j}}{$k}}\n";
      $cshl_counter++;
      $all_counter++;
    }
    if ($j eq "caltech"){
      print CALTECH_LOG "$k @{${$problems{$j}}{$k}}\n";
      $caltech_counter++;
      $all_counter++;
    }
    if ($j eq "stlace" && $j eq "brigace"){
      print STLOUIS_LOG "$j $k @{${$problems{$j}}{$k}}\n";
      $stlouis_counter++;
      $all_counter++;
    }
    
    #    print "$j $k @{${$problems{$j}}{$k}}\n";
  }
}

print "\n$all_counter problems found\n" if $verbose;



##########################################
# Tidy up and mail relevant log files
##########################################

$db->close;

print ALL_LOG "\n$all_counter problems found\n\n$0 ended at ",`date`,"\n";
print SANGER_LOG "\n$sanger_counter problems found\n\n$0 ended at ",`date`,"\n";
print CSHL_LOG "\n$cshl_counter problems found\n\n$0 ended at ",`date`,"\n";
print CALTECH_LOG "\n$caltech_counter problems found\n\n$0 ended at ",`date`,"\n";
print STLOUIS_LOG "\n$stlouis_counter problems found\n\n$0 ended at ",`date`,"\n";

close(ALL_LOG);
close(SANGER_LOG);
close(CSHL_LOG);
close(CALTECH_LOG);
close(STLOUIS_LOG);

my $all = "wormbase-dev\@wormbase.org";

my $caltech = "wormbase\@its.caltech.edu, krb\@sanger.ac.uk";
my $sanger = "wormbase\@sanger.ac.uk";
my $cshl = "stein\@cshl.org, cunningh\@cshl.edu, harris\@cshl.org, krb\@sanger.ac.uk";
my $stlouis = "dblasiar\@watson.wustl.edu, jspieth\@watson.wustl.edu, krb\@sanger.ac.uk";

&mail_maintainer("$WS_version integrity checks: Sanger","$sanger",$sanger_log) unless ($test || ($sanger_counter == 0));
&mail_maintainer("$WS_version integrity checks: CSHL","$cshl",$cshl_log) unless ($test || ($cshl_counter ==0));
&mail_maintainer("$WS_version integrity checks: Caltech","$caltech",$caltech_log) unless ($test || ($caltech_counter == 0));
&mail_maintainer("$WS_version integrity checks: St. Louis","$stlouis",$stlouis_log) unless ($test || ($stlouis_counter == 0));

exit(0);





####################################################################
#
#  T H E    S U B R O U T I N E S
#
#####################################################################


##############################################################
# Main loop to step through each sequence object which is a
# child of a genomic clone sequence
##############################################################

sub process_sequences{

  my @genome_seqs = $db->fetch(-class => 'Genome_sequence',
			       -name  => '*');
  my $class = "Sequence";

  # loop through each subseq-style name
  
  foreach my $seq (@genome_seqs){

    my @subseqs = $db->fetch(-class=>'Sequence',-name=>"$seq.*");
    
    foreach my $subseq (@subseqs){
      my $category = 0;
      
      if(!defined($subseq->at('Visible.Corresponding_protein'))){
	if(defined($subseq->at('Properties.Coding.CDS'))){  
	  my $tag = "Origin";	
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"No Corresponding_protein set for this sequence\""];
	  print "$timestamp $class:$subseq \"No Corresponding_protein tag\"\n" if $verbose;
	  $category = 1;
	}
      }    
      # only look for some specific errors if no Source tag present
      if(!defined($subseq->at('Structure.From.Source'))){  
	if(defined($subseq->at('DB_info.Database'))){
	  my $tag = "Database";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Visible.RNAi_result'))){
	  my $tag = "RNAi_result";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);	
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Visible.Locus_genomic_seq'))){
	  my $tag = "Locus_genomic_seq";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Visible.Contained_in_operon'))){
	  my $tag = "Contained_in_operon";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Allele'))){
	  my $tag = "Allele";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Visible.Matching_cDNA'))){
	  my $tag = "Matching_cDNA";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Visible.Reference'))){
	  my $tag = "Reference";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Visible.Paired_read'))){
	  my $tag = "Paired_read";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	$category = 1;
	}
	if(defined($subseq->at('DNA'))){
	  my $tag = "DNA";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Visible.Remark'))){
	  my $tag = "Remark";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Visible.Expr_pattern'))){
	  my $tag = "Expr_pattern";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Visible.Provisional_description'))){
	  my $tag = "Provisional_description";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Visible.Detailed_description'))){
	  my $tag = "Detailed_description";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Visible.Concise_description'))){
	  my $tag = "Concise_description";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Visible.Corresponding_protein'))){
	  my $tag = "Corresponding_protein";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Visible.GO_term'))){
	  my $tag = "GO_term";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Properties'))){
	  my $tag = "Properties";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Visible.Clone'))){
	  my $tag = "Clone";
	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($subseq->at('Has_allele'))){
	  my $tag = "Has_allele";
#	  my $timestamp = &get_timestamp($class, $subseq, $tag);
	  # added because AcePerl doesn't like Has_allele tag for some reason?!?
	  my $timestamp = "geneace_allele";
	  my $comment = &splice_variant_check($subseq);
	  $problems{$timestamp}{$class.":".$subseq} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$subseq \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if ($category == 0){
	  push(@other,$subseq);
	  print "$subseq - Other problem\n" if $verbose;
	}
      }    
      undef($subseq);
    }
    undef($seq);
  }



  
}

##########################################################################

sub create_log_files{
  my $database = shift;
  my $rundate    = `date +%y%m%d`; chomp $rundate;
  $all_log = "/wormsrv2/logs/current_DB_check.all.log.$rundate.$$";
  open(ALL_LOG,">$all_log") || die "can't open $all_log";
  $sanger_log = "/wormsrv2/logs/current_DB_check.sanger.log.$rundate.$$";
  open(SANGER_LOG,">$sanger_log") || die "can't open $sanger_log";
  $cshl_log = "/wormsrv2/logs/current_DB_check.cshl.log.$rundate.$$";
  open(CSHL_LOG,">$cshl_log") || die "can't open $cshl_log";
  $stlouis_log = "/wormsrv2/logs/current_DB_check.stlouis.log.$rundate.$$";
  open(STLOUIS_LOG,">$stlouis_log") || die "can't open $stlouis_log";
  $caltech_log = "/wormsrv2/logs/current_DB_check.caltech.log.$rundate.$$";
  open(CALTECH_LOG,">$caltech_log") || die "can't open $caltech_log";

  print ALL_LOG "$0 started at ",`date`,"\n";
  print ALL_LOG "This file contains information on possible errors in the current_DB database\n";
  print ALL_LOG "==========================================================================\n";

  print SANGER_LOG "$0 started at ",`date`,"\n";
  print SANGER_LOG "This file contains information on possible errors in the latest $WS_version release\n";
  print SANGER_LOG "which have been traced to the camace or geneace databases. This list is generated\n";
  print SANGER_LOG "automatically at the end of the $WS_version build process. Most items on this list will be\n";
  print SANGER_LOG "sequences created by your database which should now be replaced by splice variants\n";
  print SANGER_LOG "or removed altogether. Email wormbase\@sanger.ac.uk if you have any questions.\n";
  print SANGER_LOG "==========================================================================\n\n";

  print CSHL_LOG "$0 started at ",`date`,"\n";
  print CSHL_LOG "This file contains information on possible errors in the latest $WS_version release\n";
  print CSHL_LOG "which have been traced to the cshace database.  This list is generated\n";
  print CSHL_LOG "automatically at the end of the $WS_version build process. Most items on this list will be\n";
  print CSHL_LOG "sequences created by your database which should now be replaced by splice variants\n";
  print CSHL_LOG "or removed altogether. Email wormbase\@sanger.ac.uk if you have any questions.\n";
  print CSHL_LOG "==========================================================================\n\n";

  print STLOUIS_LOG "$0 started at ",`date`,"\n";
  print STLOUIS_LOG "This file contains information on possible errors in the latest $WS_version release\n";
  print STLOUIS_LOG "which have been traced to the stlace or brigace databases. This list is generated\n";
  print STLOUIS_LOG "automatically at the end of the $WS_version build process. Most items on this list will be\n";
  print STLOUIS_LOG "sequences created by your database which should now be replaced by splice variants\n";
  print STLOUIS_LOG "or removed altogether. Email wormbase\@sanger.ac.uk if you have any questions.\n";
  print STLOUIS_LOG "==========================================================================\n\n";

  print CALTECH_LOG "$0 started at ",`date`,"\n";
  print CALTECH_LOG "This file contains information on possible errors in the latest $WS_version release\n";
  print CALTECH_LOG "which have been traced to the citace database. This list is generated\n";
  print CALTECH_LOG "automatically at the end of the $WS_version build process. Most items on this list will be\n";
  print CALTECH_LOG "sequences created by your database which should now be replaced by splice variants\n";
  print CALTECH_LOG "or removed altogether. Email wormbase\@sanger.ac.uk if you have any questions.\n";
  print CALTECH_LOG "==========================================================================\n\n";
}


################################################
# Grab timestamps from each object using aql
################################################
sub get_timestamp{
  my $class = shift;
  my $object = shift;
  my $tag = shift;
  my $value = ""; 
  
  ($value) = $object->$tag;


  my $aql_query = "select s,s->$tag.node_session from s in object(\"$class\",\"$object\")";
  my @aql = $db->aql($aql_query);
  my $timestamp = $aql[0]->[1];
  $timestamp =~ s/\d{4}\-\d{2}\-\d{2}_\d{2}:\d{2}:\d{2}_//;

  # if timestamp is wormpub, need to look at the corresponding object to see if that has 
  # more informative timestamp
  if (($timestamp =~ m/wormpub/)&&($class eq "Sequence")){
    ($class = "Protein", $tag = "Sequence")               if ($tag eq "GO_term");
    ($class = "Protein", $tag = "Corresponding_protein")  if ($tag eq "Corresponding_DNA");
    ($class = "Expr_pattern", $tag = "Sequence")          if ($tag eq "Expr_pattern");
    ($class = "Paper", $tag = "Sequence")                 if ($tag eq "Reference");
    ($class = "Sequence", $tag = "Matching_genomic")      if ($tag eq "Matching_cDNA");
    ($class = "Allele", $tag = "Sequence")                if ($tag eq "Allele");
    ($class = "Operon", $tag = "Contains_CDS")            if ($tag eq "Contained_in_operon");
    ($class = "Locus", $tag = "Genomic_sequence")         if ($tag eq "Locus_genomic_seq");
    ($class = "RNAi", $tag = "Predicted_gene")            if ($tag eq "RNAi_result");
    ($class = "Clone", $tag = "Sequence")                 if ($tag eq "Clone");
    ($class = "Allele", $tag = "Predicted_gene")          if ($tag eq "Has_allele");

    my $aql_query = "select s,s->${tag}.node_session from s in object(\"$class\",\"$value\")";
    my @aql = $db->aql($aql_query);
    $timestamp = $aql[0]->[1];
    $timestamp =~ s/\d{4}\-\d{2}\-\d{2}_\d{2}:\d{2}:\d{2}_// unless (!defined($timestamp));

  }

  $timestamp = "unknown" if (!defined($timestamp));

  return($timestamp);

}


#################################################################

sub splice_variant_check{
  my $object = shift;
  my $splice_variant = $db->fetch(-class => 'Sequence',-name  => "${object}a");
  my $comment =  "";
  $comment =  "Splice variants now exist for this sequence" if(defined($splice_variant));
  return($comment);

}

############################################
# Check sequences connected to multiple loci
############################################

sub find_multiple_loci {


  print "\nLooking for sequences attached to multiple loci...\n" if $verbose;
  my $get_seqs_with_multiple_loci=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/get_seq_has_multiple_loci.def" quit 
EOF

  my ($seq, $locus);

  open (FH, "echo '$get_seqs_with_multiple_loci' | tace $database | ") || die "Couldn't access $database\n";
  while (<FH>){
    chomp($_);
    if ($_ =~ /^\"/){
      $_ =~ s/\"//g;
      ($seq, $locus)=split(/\s+/, $_);

      # Don't know how to get timestamps for each value next to Locus_genomic_seq tag but can
      # get the timestamp from the Locus_genomic_seq tag itself
      my $aql_query = "select s,s->Locus_genomic_seq.node_session from s in object(\"Sequence\",\"$seq\")";
      my @aql = $db->aql($aql_query);
      my $timestamp = "$aql[0]->[1]";
      $timestamp =~ s/\d{4}\-\d{2}\-\d{2}_\d{2}:\d{2}:\d{2}_//;
      $problems{$timestamp}{"Sequence".":".$seq} = ["\"connected to multiple loci\""];
    }
  }  
}
   
###########################################
sub usage {
    system ('perldoc',$0);
    exit;       
}
############################################




__END__

=pod

=head2 NAME - current_DB_check.pl

=head1 USAGE

=over 4

=item current_DB_check.pl  [-options]

=back

This script looks for errant sequence objects in current_DB database on 
wormsrv2.  These are sequence objects which have a cosmid.number name format
but do not have a Source tag.  From studying the timestamps and tags in the
offending sequence object, the script sorts the problems based on timestamp
information (i.e. from camace, cshace, build script etc.).  The script then
sends separate emails to each WormBase group listing the problems stemming
from their database(s).  If no source database can be traced, problems will
be sent to Sanger.

This script also checks valid looking sequence objects to see if they have
a Corresponding_protein tag.

current_DB_check.pl MANDATORY arguments:

=over 4

=item none

=back

current_DB_check.pl  OPTIONAL arguments:

=over 4

=item -d, database

By default this script will check /wormsrv2/current_DB, but you can use the 
-d flag to checkn another database (i.e. /wormsrv2/autoace)

=item -t, test

Will run all of the script as normal but will not send any emails.  Separate 
log files are still written to /wormsrv2/logs

=item -v, verbose mode

Turning on verbose mode will output errors to the screen as it finds them, it still
writes a log file, but it doesn't email wormbase-dev

=item -h, Help

This help.

=back

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Keith Bradnam (krb@sanger.ac.uk)

=back

=cut
