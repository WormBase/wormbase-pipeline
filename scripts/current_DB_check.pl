#!/usr/local/bin/perl5.8.0 -w
# 
# current_DB_check.pl
#
# by Keith Bradnam
#
# Script to run consistency checks on the current_DB database
# to look for bogus sequence entries
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2004-08-02 16:12:51 $


use strict;
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Getopt::Long;
use GENEACE::Geneace;



##############################
# command-line options       #
##############################

my $help;                # Help/Usage page
my $verbose;             # turn on extra output
my $database;            # path to database
my $debug;               # For sending output to just one person
my $test;                # for running in test mode
my $maintainers = "All"; # log file recipients


GetOptions ("database=s" => \$database,
            "verbose"    => \$verbose,
            "test"       => \$test,
            "debug=s"    => \$debug,
            "help"       => \$help);

&usage if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}



#############################################################################
# Set database details, default is /nfs/disk100/wormpub/DATABASES/current_DB
##############################################################################

$database = "/nfs/disk100/wormpub/DATABASES/current_DB/" if (!$database);

my $tace = &tace;   # tace executable path

my $db = Ace->connect(-path  => "$database",
		      -program =>$tace) || do { print ALL_LOG "Connection failure: ",Ace->error; die();};

print "Checking $database\n\n";



##############################################
# Initialise variables, set up log files
###############################################

my %problems;         # store problems in a double hash, first key being timestamp name
my @other;            # store uncategorised problems
my @gene_id_error;    # store out-of-date gene id in acefiles of each group for the build: subdir of /wormsrv2/wormbase

# make one global log file but also split for different groups
my $all_log;
my $sanger_log;
my $cshl_log;
my $stlouis_log;
my $caltech_log;
my $WS_version  = &get_wormbase_version_name;

&create_log_files;


###################################
# Report errors to log file
###################################

# Count all the problems!
my $sanger_counter = 0;
my $caltech_counter = 0;
my $cshl_counter = 0;
my $stlouis_counter = 0;
my $all_counter = 0;



# Checks for CDS connected to multiple loci
print "\nChecking for CDSs connected to multiple loci\n" if ($verbose);
&find_multiple_loci;

# Search everything else, one sequence at a time
print "\nChecking CDS objects\n" if ($verbose);
&process_sequences;

# Check out-of-date gene id
print "\nChecking for connections to non-live Gene objects\n" if ($verbose);
&find_out_of_date_gene_id;

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

foreach (@gene_id_error){
  print ALL_LOG $_;
  print $_ if $verbose;
  $all_counter++;
  if ($_ =~ /^csh/ ){ $cshl_counter++; print CSHL_LOG $_ }
  if ($_ =~ /^camace|^misc/ ){ $sanger_counter++; print SANGER_LOG $_ }
  if ($_ =~ /^caltech/ ){ $caltech_counter++; print CALTECH_LOG $_ }
  if ($_ =~ /^stlace/ ){ $stlouis_counter++; print STLOUIS_LOG $_ }
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
my $cshl = "stein\@cshl.org, harris\@cshl.org, chenn\@cshl.edu, krb\@sanger.ac.uk";
my $stlouis = "dblasiar\@watson.wustl.edu, jspieth\@watson.wustl.edu, krb\@sanger.ac.uk";

if($debug){
  $caltech = $sanger = $cshl = $stlouis = $debug;
}

&mail_maintainer("$WS_version database checks: Sanger","$sanger",$sanger_log)    unless ($test || ($sanger_counter == 0));
&mail_maintainer("$WS_version database checks: CSHL","$cshl",$cshl_log)          unless ($test || ($cshl_counter ==0));
&mail_maintainer("$WS_version database checks: Caltech","$caltech",$caltech_log) unless ($test || ($caltech_counter == 0));
&mail_maintainer("$WS_version database checks: WashU","$stlouis",$stlouis_log)   unless ($test || ($stlouis_counter == 0));

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
  my $class = "CDS";

  # loop through each subseq-style name

  foreach my $seq (@genome_seqs){

    my @CDSs = $db->fetch(-class=>'CDS',-name=>"$seq.*");

    foreach my $CDS (@CDSs){
      our $category = 0;

      # only look for some errors for CDS objects with no Method (most likely to be errant objects)
      if(!defined($CDS->at('Method'))){  
	if(defined($CDS->at('DB_info.Database'))){
	  my $tag = "Database";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.RNAi_result'))){
	  my $tag = "RNAi_result";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);	
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.Gene'))){
	  my $tag = "Gene";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.Contained_in_operon'))){
	  my $tag = "Contained_in_operon";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Allele'))){
	  my $tag = "Allele";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.Matching_cDNA'))){
	  my $tag = "Matching_cDNA";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.Reference'))){
	  my $tag = "Reference";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.Paired_read'))){
	  my $tag = "Paired_read";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	$category = 1;
	}
	if(defined($CDS->at('DNA'))){
	  my $tag = "DNA";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.Remark'))){
	  my $tag = "Remark";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.Expr_pattern'))){
	  my $tag = "Expr_pattern";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.Provisional_description'))){
	  my $tag = "Provisional_description";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.Detailed_description'))){
	  my $tag = "Detailed_description";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.Concise_description'))){
	  my $tag = "Concise_description";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.Corresponding_protein'))){
	  my $tag = "Corresponding_protein";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.GO_term'))){
	  my $tag = "GO_term";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Properties'))){
	  my $tag = "Properties";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.Clone'))){
	  my $tag = "Clone";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.Alleles'))){
	  my $tag = "Alleles";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.Y2H_target'))){
	  my $tag = "Y2H_target";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	
	if(defined($CDS->at('Visible.Y2H_bait'))){
	  my $tag = "Y2H_bait";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.SAGE_transcript'))){
	  my $tag = "SAGE_transcript";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}

	if(defined($CDS->at('Visible.Drives_Transgene'))){
	  my $tag = "Drives_Transgene";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}

	if(defined($CDS->at('Visible.Transgene_product'))){
	  my $tag = "Transgene_product";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}
	if(defined($CDS->at('Visible.Microarray_results'))){
	  my $tag = "Microarray_results";
	  my $timestamp = &get_timestamp($class, $CDS, $tag);
	  my $comment = &splice_variant_check($CDS);
	  $problems{$timestamp}{$class.":".$CDS} = ["\"$tag tag is creating this sequence. $comment\""];
	  print "$timestamp $class:$CDS \"$tag tag is creating this sequence. $comment\"\n" if $verbose;
	  $category = 1;
	}


	if ($category == 0){
	  push(@other,$CDS);
	  print "$CDS - Other problem\n" if $verbose;
	}
      }
      $CDS->DESTROY();
    }
    $seq->DESTROY();
  }
}

sub find_out_of_date_gene_id {

  my @models = `cat /wormsrv2/autoace/wspec/models.wrm`;
  my ($class, %class_has_Gene);

  # ----- checks which ?class used ?Gene in model.wrm file and build a hash
  foreach my $line (@models) {
    chomp $line;
    if ($line =~ /^\?(\w+)\s+.+/) {
      $class = $1;
    }
    if ($line =~ /\?Gene.+/) {
      $class_has_Gene{$class}++;
    }
  }

  # ----- get latest live gene id / non_live gene id from geneace
  my $ga = init Geneace();
  my ($Live_gene_id, $Non_Live_gene_id) = $ga -> gene_id_status(); # hash ref

  # ----- loop through each acefile of all groups subdir of at /wormsrv2/wormbase and
  #       see if there is WBGxxxxxxx. If yes, check if gene id is up-to-date

  foreach my $group ("stlace", "camace", "csh", "caltech", "misc") {
    my @config = `cat /wormsrv2/autoace_config/autoace.config`;
    my $error =();
    my ($acefile, $exist, $obj);
    foreach my $line (@config){
      if ($line =~ /^$group\s+(.+\.ace)\s+(\w+)\s+/){
	$acefile = $1;
	$acefile =~ s/\s+//g;
	$acefile = "/wormsrv2/wormbase/$group/$acefile";
	$exist = `grep "WBGene" $acefile`; # check only acefiles that have WBGxxxxxxx info
	
	if ( exists $class_has_Gene{$2} && $exist ){  # $2 is ?class
	  $class = ();
	  my @acefile = `cat $acefile`;
	  foreach my $line (@acefile){
	    chomp $line;
	    if ( $line =~ /(\w+) : \"(.+)\" -O.+/ ){
	      $class = $1;
	      $obj = $2;
	    }
	    if ( $line =~ /(WBGene\d+)/ ){
	      if ( !exists $Live_gene_id->{$1} && exists $Non_Live_gene_id->{$1} ){
		push (@gene_id_error, "$group $class: $obj \"attached gene id $1 has been merged into $Non_Live_gene_id->{$1}\"\n");
	      }
	      if ( !exists $Live_gene_id->{$1} && !exists $Non_Live_gene_id->{$1} ){
		push (@gene_id_error, "$group $class: $obj \"attached gene id $1 is invalid\"\n");
	      }
	    }
	  }
	}
      }
    }
  }
}

##########################################################################

sub create_log_files{

  my $rundate = &rundate;
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
  if (($timestamp =~ m/wormpub/)&&($class eq "CDS")){
    ($class = "Protein", $tag = "Sequence")               if ($tag eq "GO_term");
    ($class = "Protein", $tag = "Corresponding_protein")  if ($tag eq "Corresponding_CDS");
    ($class = "Expr_pattern", $tag = "CDS")               if ($tag eq "Expr_pattern");
    ($class = "Paper", $tag = "CDS")                      if ($tag eq "Reference");
    ($class = "Sequence", $tag = "Matching_CDS")          if ($tag eq "Matching_cDNA");
    ($class = "Allele", $tag = "Predicted_gene")          if ($tag eq "Alleles");
    ($class = "Operon", $tag = "Contains_CDS")            if ($tag eq "Contained_in_operon");
    ($class = "Gene", $tag = "Corresponding_CDS")         if ($tag eq "Gene");
    ($class = "RNAi", $tag = "CDS")                       if ($tag eq "RNAi_result");
    ($class = "Clone", $tag = "Sequence")                 if ($tag eq "Clone");
    ($class = "Allele", $tag = "CDS")                     if ($tag eq "Allele");
    ($class = "Y2H", $tag = "Target_overlapping_CDS")     if ($tag eq "Y2H_target");
    ($class = "Y2H", $tag = "Bait_overlapping_CDS")       if ($tag eq "Y2H_bait");
    ($class = "SAGE_transcript", $tag = "Predicted_gene") if ($tag eq "SAGE_transcript");
    ($class = "Transgene", $tag = "Driven_by_CDS_promoter") if ($tag eq "Drives_transgene");
    ($class = "Transgene", $tag = "CDS")                  if ($tag eq "Transgene_product");
    ($class = "Microarray_results", $tag = "CDS")                  if ($tag eq "Microarray_results");

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
  my $splice_variant = $db->fetch(-class => 'elegans_CDS',-name  => "${object}a");
  my $comment =  "";
  $comment =  "Splice variants now exist for this sequence" if(defined($splice_variant));
  return($comment);

}

############################################
# Check sequences connected to multiple loci
############################################

sub find_multiple_loci {

  my $get_seqs_with_multiple_loci=<<EOF;
  Table-maker -p "/wormsrv2/autoace/wquery/CDSs_with_multiple_genes.def" quit 
EOF

  my ($cds, $gene);

  open (FH, "echo '$get_seqs_with_multiple_loci' | tace $database | ") || die "Couldn't access $database\n";
  while (<FH>){
    chomp($_);
    if ($_ =~ /^\"/){
      $_ =~ s/\"//g;
      ($cds, $gene)=split(/\s+/, $_);

      # Don't know how to get timestamps for each value next to Gene tag but can
      # get the timestamp from the Corresponding_CDS tag itself
      my $aql_query = "select s,s->Gene.node_session from s in object(\"Corresponding_CDS\",\"$cds\")";

      my @aql = $db->aql($aql_query);
      my $timestamp = "$aql[0]->[1]";
      $timestamp =~ s/\d{4}\-\d{2}\-\d{2}_\d{2}:\d{2}:\d{2}_//;
      $problems{$timestamp}{"CDS".":".$cds} = ["\"connected to multiple genes\""];
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

This script looks for errant CDS objects in current_DB database on wormsrv2.  
These are ?CDS objects which have a cosmid.number name format but do not have 
a Sequence tag to connect them to their parent.  From studying the timestamps 
and tags in the offending sequence object, the script sorts the problems based 
on timestamp information (i.e. from camace, cshace, build script etc.).  The 
script then sends separate emails to each WormBase group listing the problems 
stemming from their database(s).  If no source database can be traced, problems 
will be sent to Sanger.

This script also checks valid looking sequence objects to see if they have
a Corresponding_protein tag.

current_DB_check.pl MANDATORY arguments:

=over 4

=item none

=back

current_DB_check.pl  OPTIONAL arguments:

=over 4

=item -d, database

By default this script will check ~wormpub/DATABASES/current_DB, but you can use the 
-d flag to check another database (i.e. /wormsrv2/autoace)

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
