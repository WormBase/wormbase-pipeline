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
# Last updated on: $Date: 2002-10-15 13:34:56 $

use Ace;
use lib "/wormsrv2/scripts/"; 
use Wormbase;
use strict;
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
my $verbose;
$verbose = "ON" if ($opt_v);
my $test;
$test = "ON" if ($opt_t);

########################################################
# Set database details, default is /wormsrv2/current_DB
########################################################
my $database = "/wormsrv2/current_DB/";
$database = $opt_d if ($opt_d);
our $tace = glob("~wormpub/ACEDB/bin.ALPHA_4/tace");   # tace executable path
my $db = Ace->connect(-path  => "$database",
		      -program =>$tace) || do { print ALL_LOG "Connection failure: ",Ace->error; die();};

print "Checking $database\n\n";



##############################################
# Initialise variables, set up log files
###############################################
my $errors;
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


############################################################
# Check ?Sequence class (including predicted genes)
############################################################

print "\nLooking for spurious sequences\n";

my @genome_seqs = $db->fetch(-class => 'Genome_sequence',
		      -name  => '*');
my $class = "Sequence";

# loop through each subseq-style name

foreach my $seq (@genome_seqs){
  my @subseqs = $db->fetch(-class=>'Sequence',-name=>"$seq.*");

  foreach my $subseq (@subseqs){
    my $category = 0;

    if(!defined($subseq->at('Visible.Corresponding_protein'))&&defined($subseq->at('Properties.Coding.CDS'))){  
      $errors++;     
      my $tag = "Origin";
      my $timestamp = &get_timestamp($class, $subseq, $tag);
      $problems{$timestamp}{$class.":".$subseq} = ["\"No Corresponding_protein set for this sequence\""];
      print "$timestamp $class:$subseq \"No Corresponding_protein tag\"\n" if $verbose;
      $category = 1;
    }
    
    # only look for some specific errors if no Source tag present
    if(!defined($subseq->at('Structure.From.Source'))){  
      $errors++;     
      
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
      if ($category == 0){
	push(@other,$subseq);
	print "$subseq - Other problem\n" if $verbose;
      }
    }    
  }
}

###################################
# Report errors to log file
###################################
print "\n\n$errors errors found\n" if $verbose;
print ALL_LOG "\n$errors errors found\n\n\n";



foreach my $i (@other){
  print ALL_LOG "Other - $i\n";
  print SANGER_LOG "Other - $i\n";
  print "Other - $i\n" if $verbose;
}

foreach my $j (keys(%problems)){
  foreach my $k (keys(%{$problems{$j}})){
    print ALL_LOG "$j $k @{${$problems{$j}}{$k}}\n";
    print SANGER_LOG "$j $k @{${$problems{$j}}{$k}}\n" if ($j ne "caltech" && $j ne "csh" && $j ne "stlace" && $j ne "brigace");
    print CSHL_LOG "$k @{${$problems{$j}}{$k}}\n" if ($j eq "csh");
    print CALTECH_LOG "$k @{${$problems{$j}}{$k}}\n" if ($j eq "caltech");
    print STLOUIS_LOG "$j $k @{${$problems{$j}}{$k}}\n" if ($j eq "stlace" && $j eq "brigace");
							  
#    print "$j $k @{${$problems{$j}}{$k}}\n";
  }
}


##########################################
# Tidy up and mail relevant log files
##########################################

$db->close;
print ALL_LOG "\n$0 ended at ",`date`,"\n";
print SANGER_LOG "\n$0 ended at ",`date`,"\n";
print CSHL_LOG "\n$0 ended at ",`date`,"\n";
print CALTECH_LOG "\n$0 ended at ",`date`,"\n";
print STLOUIS_LOG "\n$0 ended at ",`date`,"\n";

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

&mail_maintainer("End of build database integrity checks: Sanger","$sanger",$sanger_log) unless $test;
&mail_maintainer("End of build database integrity checks: CSHL","$cshl",$cshl_log) unless $test;
&mail_maintainer("End of build database integrity checks: Caltech","$caltech",$caltech_log) unless $test;
&mail_maintainer("End of build database integrity checks: St. Louis","$stlouis",$stlouis_log) unless $test;

exit(0);

################################################


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
    ($class = "Protein")              if ($tag eq "GO_term");
    ($tag = "Sequence")               if ($tag eq "GO_term");
    ($class = "Protein")              if ($tag eq "Corresponding_DNA");
    ($tag = "Corresponding_protein")  if ($tag eq "Corresponding_DNA");
    ($class = "Expr_pattern")         if ($tag eq "Expr_pattern");
    ($tag = "Sequence")               if ($tag eq "Expr_pattern");
    ($class = "Paper")                if ($tag eq "Reference");
    ($tag = "Sequence")               if ($tag eq "Reference");
    ($tag = "Matching_genomic")       if ($tag eq "Matching_cDNA");
    ($class = "Allele")               if ($tag eq "Allele");
    ($tag = "Sequence")               if ($tag eq "Allele");
    ($class = "Operon")               if ($tag eq "Contained_in_operon");
    ($tag = "Contains_CDS")           if ($tag eq "Contained_in_operon");
    ($class = "Locus")                if ($tag eq "Genomic_sequence");
    ($tag = "Locus_genomic_seq")      if ($tag eq "Genomic_sequence");
    ($class = "RNAi")                 if ($tag eq "RNAi_result");
    ($tag = "Predicted_gene")         if ($tag eq "RNAi_result");
   

    my $aql_query = "select s,s->${tag}.node_session from s in object(\"$class\",\"$value\")";
    my @aql = $db->aql($aql_query);
    $timestamp = $aql[0]->[1];
    $timestamp =~ s/\d{4}\-\d{2}\-\d{2}_\d{2}:\d{2}:\d{2}_// unless (!defined($timestamp));

  }

  $timestamp = "unknown" if (!defined($timestamp));

  return($timestamp);

}

sub splice_variant_check{
  my $object = shift;
  my $splice_variant = $db->fetch(-class => 'Sequence',-name  => "${object}a");
  my $comment =  "";
  $comment =  "Splice variants now exist for this sequence" if(defined($splice_variant));
  return($comment);

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

Will run all of the script as normal but will not send any emails.  Separate l
og files are still written to /wormsrv2/logs

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
