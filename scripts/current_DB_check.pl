#!/usr/local/bin/perl5.8.0 -w
# 
# current_DB_check.pl
#
# by Keith Bradnam
#
# Script to run consistency checks on the current_DB database
# to look for bogus sequence entries
#
# Last updated by: $Author: gw3 $

# Last updated on: $Date: 2006-01-03 17:35:31 $

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
use GENEACE::Geneace;



##############################
# command-line options       #
##############################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my $database;            # path to database

GetOptions ("database=s" => \$database,
            "verbose"    => \$verbose,
            "test"       => \$test,
            "debug=s"    => \$debug,
            "help"       => \$help,
	    "store"      => \$store,
);


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);



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
my $WS_version  = $wormbase->get_wormbase_version_name;

&create_log_files;


#############################################################################
# Set database details, default is current_DB
##############################################################################

$database = $wormbase->databases('current') if (!$database);
my $tace = $wormbase->tace;        # TACE PATH

my $db = Ace->connect(-path  => "$database",
		      -program =>$tace) || $all_log->log_and_die("Connection failure: " . Ace->error);
print "Checking $database\n\n";



###################################
# Report errors to log file
###################################

# Count all the problems!
my $sanger_counter = 0;
my $caltech_counter = 0;
my $cshl_counter = 0;
my $stlouis_counter = 0;
my $all_counter = 0;


# Check out-of-date gene id
print "\nChecking for connections to non-live Gene objects\n" if ($verbose);
&find_out_of_date_gene_id;

# look for CDS, Pseudogene, and Transcript objects created by XREF
&process_genes;



# Count problems and print output

foreach (@gene_id_error){
  $all_log->write_to( $_ );
  print $_ if $verbose;
  $all_counter++;
  if ($_ =~ /^csh/ )         { $cshl_counter++;    $cshl_log->write_to( $_ )}
  if ($_ =~ /^camace|^misc/ ){ $sanger_counter++;  $sanger_log->writeto( $_ )}
  if ($_ =~ /^caltech/ )     { $caltech_counter++; $caltech_log->write_to( $_ )}
  if ($_ =~ /^stlace/ )      { $stlouis_counter++; $stlouis_log->write_to( $_ )}
}


foreach my $j (keys(%problems)){
  foreach my $k (keys(%{$problems{$j}})){
    $all_log->write_to("ERROR: $j $k @{${$problems{$j}}{$k}}\n");
    if ($j ne "caltech" && $j ne "csh" && $j ne "stlace" && $j ne "brigace"){
      $sanger_log->writeto( "ERROR: $j $k @{${$problems{$j}}{$k}}\n");
      $sanger_counter++;
      $all_counter++;
    }
    if ($j eq "csh"){
      $cshl_log->write_to( "ERROR: $k @{${$problems{$j}}{$k}}\n");
      $cshl_counter++;
      $all_counter++;
    }
    if ($j eq "caltech"){
      $caltech_log->write_to( "ERROR: $k @{${$problems{$j}}{$k}}\n");
      $caltech_counter++;
      $all_counter++;
    }
    if ($j eq "stlace" && $j eq "brigace"){
      $stlouis_log->write_to( "ERROR: $j $k @{${$problems{$j}}{$k}}\n");
      $stlouis_counter++;
      $all_counter++;
    }
    #    print "$j $k @{${$problems{$j}}{$k}}\n";
  }
}

foreach my $i (@other){
  $all_log->write_to( "ERROR: Other - $i\n");
  $sanger_log->writeto( "ERROR: Other - $i\n");
  print "Other - $i\n" if $verbose;
  $all_counter++;
  $sanger_counter++;
}



print "\n$all_counter problems found\n" if $verbose;



##########################################
# Tidy up and mail relevant log files
##########################################

$db->close;

$all_log->write_to( "\n$all_counter problems found\n\n$0 ended at ",`date`,"\n");
$sanger_log->writeto( "\n$sanger_counter problems found\n\n$0 ended at ",`date`,"\n");
$cshl_log->write_to( "\n$cshl_counter problems found\n\n$0 ended at ",`date`,"\n");
$caltech_log->write_to( "\n$caltech_counter problems found\n\n$0 ended at ",`date`,"\n");
$stlouis_log->write_to( "\n$stlouis_counter problems found\n\n$0 ended at ",`date`,"\n");

my $all = "wormbase-dev\@wormbase.org";

my $caltech = "wormbase\@its.caltech.edu";
my $sanger  = "wormbase\@sanger.ac.uk";
my $cshl    = "stein\@cshl.org, harris\@cshl.org, chenn\@cshl.edu";
my $stlouis = "dblasiar\@watson.wustl.edu, jspieth\@watson.wustl.edu";

if($debug){
  $caltech = $sanger = $cshl = $stlouis = $debug;
}

# +++ check this! do we really want to mail wormbase-dev with these problems?
#$all_log->mail($all, "$WS_version database checks: All")             unless ($test || ($all_counter == 0));

$sanger_log->mail($sanger, "$WS_version database checks: Sanger")    unless ($test || ($sanger_counter == 0));
$cshl_log->mail($cshl, "$WS_version database checks: CSHL")          unless ($test || ($cshl_counter == 0));
$caltech_log->mail($caltech, "$WS_version database checks: Caltech") unless ($test || ($caltech_counter == 0));
$stlouis_log->mail($stlouis, "$WS_version database checks: WashU")   unless ($test || ($stlouis_counter == 0));

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);



####################################################################
#
#  T H E    S U B R O U T I N E S
#
#####################################################################


sub process_genes{

  # loop through each class  
  my @classes_to_check = ("CDS","Transcript","Pseudogene");

  foreach my $class (@classes_to_check){
    print "\n\n\nCHECKING $class CLASS\n\n\n" if ($verbose);

    my @genes = $db->fetch(-query=>"Find $class WHERE !Method");


    foreach my $gene (@genes){
      my $category = 0;

      print "\n\n*$gene*\n" if ($verbose);
      # only look for some errors for CDS objects with no Method (most likely to be errant objects)
      if(!defined($gene->at('Method'))){  
	if(defined($gene->at('DB_info.Database'))){
	  my $tag = "Database";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.RNAi_result'))){
	  my $tag = "RNAi_result";
	  my $timestamp = &get_timestamp($class, $gene, $tag);	
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Gene'))){
	  my $tag = "Gene";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Contained_in_operon'))){
	  my $tag = "Contained_in_operon";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Allele'))){
	  my $tag = "Allele";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Matching_cDNA'))){
	  my $tag = "Matching_cDNA";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Reference'))){
	  my $tag = "Reference";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Paired_read'))){
	  my $tag = "Paired_read";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('DNA'))){
	  my $tag = "DNA";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Remark'))){
	  my $tag = "Remark";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Expr_pattern'))){
	  my $tag = "Expr_pattern";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Provisional_description'))){
	  my $tag = "Provisional_description";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Detailed_description'))){
	  my $tag = "Detailed_description";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Concise_description'))){
	  my $tag = "Concise_description";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Corresponding_protein'))){
	  my $tag = "Corresponding_protein";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.GO_term'))){
	  my $tag = "GO_term";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Properties'))){
	  my $tag = "Properties";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Clone'))){
	  my $tag = "Clone";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Alleles'))){
	  my $tag = "Alleles";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Y2H_target'))){
	  my $tag = "Y2H_target";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	
	if(defined($gene->at('Visible.Y2H_bait'))){
	  my $tag = "Y2H_bait";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.SAGE_transcript'))){
	  my $tag = "SAGE_transcript";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}

	if(defined($gene->at('Visible.Drives_Transgene'))){
	  my $tag = "Drives_Transgene";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}

	if(defined($gene->at('Visible.Transgene_product'))){
	  my $tag = "Transgene_product";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Microarray_results'))){
	  my $tag = "Microarray_results";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}
	if(defined($gene->at('Visible.Corresponding_oligo_set'))){
	  my $tag = "Corresponding_oligo_set";
	  my $timestamp = &get_timestamp($class, $gene, $tag);
	  my $comment = &check_fate_of_gene($gene,$class);
	  $problems{$timestamp}{$class.":".$gene} = ["\'$tag\' tag is creating this sequence. $comment"];
	  print "$timestamp $class:$gene \'$tag\' tag is creating this sequence. $comment\n" if $verbose;
	  $category = 1;
	}



	if ($category == 0){
	  push(@other,$gene);
	  print "$gene - Other problem\n" if $verbose;
	}
      }
      $gene->DESTROY();
    }
  }
}

sub find_out_of_date_gene_id {
  
  my $ace_dir = $wormbase->autoace;     # AUTOACE DATABASE DIR
  my $misc_static_dir = $wormbase->misc_static;

  my @models = `cat $ace_dir/wspec/models.wrm`;
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
  
  my $db = Ace->connect(-path  => $database,
			-program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

  
  # get list of live genes and add to hash
  my @genes = $db->fetch(-query=>"Find Gene WHERE !Live");
  my %genes;
  foreach my $gene (@genes){
    my $fate;
    if($gene->Merged_into){
      $fate = $gene->Merged_into;
      $genes{$gene->name} = "$fate";
    }
    else{
      $genes{$gene->name} = "Other fate";
    }
  }
  $db->close;

  foreach my $group ("stlace", "camace", "csh", "caltech", "misc") {
    print "\n\nProcessing $group data\n\n" if ($verbose);
    my @config = `cat /wormsrv2/autoace_config/autoace.config`;
    my $error =();
    my ($acefile, $exist, $obj);
    foreach my $config_line (@config){
      if ($config_line =~ /^$group\s+(.+\.ace)\s+(\w+)\s+/){
	$acefile = $1;
	$acefile =~ s/\s+//g;

	# fix path if group is misc
	if($group eq "misc"){
	  $acefile = "$misc_static_dir/$1";
	} else {
	  $acefile = $wormbase->primary($group) . "/$acefile";
	}

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
	      # skip gene history lines which will by definition be non-Live genes
	      next if ($line =~ m/Gene_history/);

	      if (exists $genes{$1}){
		# ignore CDS objects if there is no Merged_into tag, these are probably all
		# cases of CDSs being made into Transposon CDSs
		unless(($class eq "CDS") && ($genes{$1} eq "Other fate")){
		  print "$class : $obj\n$line\nGene $1 is no longer live, gene was merged into $genes{$1}\n\n" if ($verbose); 
		  if($genes{$1} eq "Other fate"){		    
		    push (@gene_id_error, "$group $class:$obj connects to gene $1 which is no longer live and has not been merged into any other gene.  Please investigate\n");
		  }
		  else{
		    push (@gene_id_error, "$group $class:$obj connects to gene $1 which is no longer live, gene was merged into $genes{$1}\n");
		  }
		}
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

  my $rundate = $wormbase->rundate;
  my $logs_dir = $wormbase->logs;

  $all_log = "$logs_dir/current_DB_check.all.log.$rundate.$$";
  #open(ALL_LOG,">$all_log") || die "can't open $all_log";
  $all_log = Log_files->make_log($all_log);

  $sanger_log = "$logs_dir/current_DB_check.sanger.log.$rundate.$$";
  #open(SANGER_LOG,">$sanger_log") || die "can't open $sanger_log";
  $sanger_log = Log_files->make_log($sanger_log);

  $cshl_log = "$logs_dir/current_DB_check.cshl.log.$rundate.$$";
  #open(CSHL_LOG,">$cshl_log") || die "can't open $cshl_log";
  $cshl_log = Log_files->make_log($cshl_log);

  $stlouis_log = "$logs_dir/current_DB_check.stlouis.log.$rundate.$$";
  #open(STLOUIS_LOG,">$stlouis_log") || die "can't open $stlouis_log";
  $stlouis_log = Log_files->make_log($stlouis_log);

  $caltech_log = "$logs_dir/current_DB_check.caltech.log.$rundate.$$";
  #open(CALTECH_LOG,">$caltech_log") || die "can't open $caltech_log";
  $caltech_log = Log_files->make_log($caltech_log);

  $all_log->write_to( "$0 started at ",`date`,"\n");
  $all_log->write_to( "This file contains information on possible errors in the current_DB database\n");
  $all_log->write_to( "==========================================================================\n");

  $sanger_log->writeto( "$0 started at ",`date`,"\n");
  $sanger_log->writeto( "This file contains information on possible errors in the latest $WS_version release\n");
  $sanger_log->writeto( "which have been traced to the camace or geneace databases. This list is generated\n");
  $sanger_log->writeto( "automatically at the end of the $WS_version build process. Most items on this list will be\n");
  $sanger_log->writeto( "sequences created by your database which should now be replaced by splice variants\n");
  $sanger_log->writeto( "or removed altogether. Email wormbase\@sanger.ac.uk if you have any questions.\n");
  $sanger_log->writeto( "==========================================================================\n\n");

  $cshl_log->write_to( "$0 started at ",`date`,"\n");
  $cshl_log->write_to( "This file contains information on possible errors in the latest $WS_version release\n");
  $cshl_log->write_to( "which have been traced to the cshace database.  This list is generated\n");
  $cshl_log->write_to( "automatically at the end of the $WS_version build process. Most items on this list will be\n");
  $cshl_log->write_to( "sequences created by your database which should now be replaced by splice variants\n");
  $cshl_log->write_to( "or removed altogether. Email wormbase\@sanger.ac.uk if you have any questions.\n");
  $cshl_log->write_to( "==========================================================================\n\n");

  $stlouis_log->write_to( "$0 started at ",`date`,"\n");
  $stlouis_log->write_to( "This file contains information on possible errors in the latest $WS_version release\n");
  $stlouis_log->write_to( "which have been traced to the stlace or brigace databases. This list is generated\n");
  $stlouis_log->write_to( "automatically at the end of the $WS_version build process. Most items on this list will be\n");
  $stlouis_log->write_to( "sequences created by your database which should now be replaced by splice variants\n");
  $stlouis_log->write_to( "or removed altogether. Email wormbase\@sanger.ac.uk if you have any questions.\n");
  $stlouis_log->write_to( "==========================================================================\n\n");

  $caltech_log->write_to( "$0 started at ",`date`,"\n");
  $caltech_log->write_to( "This file contains information on possible errors in the latest $WS_version release\n");
  $caltech_log->write_to( "which have been traced to the citace database. This list is generated\n");
  $caltech_log->write_to( "automatically at the end of the $WS_version build process. Most items on this list will be\n");
  $caltech_log->write_to( "sequences created by your database which should now be replaced by splice variants\n");
  $caltech_log->write_to( "or removed altogether. Email wormbase\@sanger.ac.uk if you have any questions.\n");
  $caltech_log->write_to( "==========================================================================\n\n");
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
    ($class = "Protein", $tag = "Sequence")                 if ($tag eq "GO_term");
    ($class = "Protein", $tag = "Corresponding_protein")    if ($tag eq "Corresponding_CDS");
    ($class = "Expr_pattern", $tag = "CDS")                 if ($tag eq "Expr_pattern");
    ($class = "Paper", $tag = "CDS")                        if ($tag eq "Reference");
    ($class = "Sequence", $tag = "Matching_CDS")            if ($tag eq "Matching_cDNA");
    ($class = "Allele", $tag = "Predicted_gene")            if ($tag eq "Alleles");
    ($class = "Operon", $tag = "Contains_CDS")              if ($tag eq "Contained_in_operon");
    ($class = "Gene", $tag = "Corresponding_CDS")           if ($tag eq "Gene");
    ($class = "RNAi", $tag = "CDS")                         if ($tag eq "RNAi_result");
    ($class = "Clone", $tag = "Sequence")                   if ($tag eq "Clone");
    ($class = "Allele", $tag = "CDS")                       if ($tag eq "Allele");
    ($class = "Y2H", $tag = "Target_overlapping_CDS")       if ($tag eq "Y2H_target");
    ($class = "Y2H", $tag = "Bait_overlapping_CDS")         if ($tag eq "Y2H_bait");
    ($class = "SAGE_transcript", $tag = "Predicted_gene")   if ($tag eq "SAGE_transcript");
    ($class = "Transgene", $tag = "Driven_by_CDS_promoter") if ($tag eq "Drives_transgene");
    ($class = "Transgene", $tag = "CDS")                    if ($tag eq "Transgene_product");
    ($class = "Microarray_results", $tag = "CDS")           if ($tag eq "Microarray_results");
    ($class = "Oligo_set", $tag = "Overlaps_CDS")           if ($tag eq "Corresponding_oligo_set");

    my $aql_query = "select s,s->${tag}.node_session from s in object(\"$class\",\"$value\")";
    my @aql = $db->aql($aql_query);
    $timestamp = $aql[0]->[1];
    $timestamp =~ s/\d{4}\-\d{2}\-\d{2}_\d{2}:\d{2}:\d{2}_// unless (!defined($timestamp));

  }

  $timestamp = "unknown" if (!defined($timestamp));

  return($timestamp);

}

#################################################################
# Check whether corresponding gene is dead or there are isoforms 
# or data should now be connected to a different class etc.
#################################################################

sub check_fate_of_gene{
  my $object = shift;
  my $class = shift;

  # treat object name to remove splice variant suffixes
  my $new_name = $object;
  if($new_name =~ m/[a-z]$/){
    $new_name =~ s/[a-z]$//;
  }

  # store comment to return to parent subroutine
  my $comment =  "";

  # get gene object via Gene_name class
  my $gene_name = $db->fetch(-class => 'Gene_name',-name  => "$new_name");

  # exit early if can't get valid gene name
  if(!defined($gene_name)){
    $comment = "$object apppears to have no corresponding ?Gene_name object...possible typo in $object name?";
    return($comment);
  }
  else{
    # now query corresponding gene object
    my $gene = $gene_name->Sequence_name_for;
    
    # test Public_name_for if not defined as gene may be dead (no Sequence_name_for field set)
    if(!defined($gene)){
      $gene = $gene_name->Public_name_for;
    }
    
    if(defined($gene)){
      # test if gene is Live
      if(!defined($gene->at("Identity.Live"))){
	$comment = "$object connects to a dead gene (no 'Live' tag)...maybe gene has been merged into another?";
	return($comment);
      }
      # else remind what valid objects should connect to the gene
      else{
	my @CDSs        = $gene->Corresponding_CDS;
	my @transcripts = $gene->Corresponding_transcript;
	my @pseudogenes = $gene->Corresponding_pseudogene;

	if(@CDSs || @pseudogenes || @transcripts){
	  $comment .= "$class : $object should be replaced by: ";
	}
	if(@CDSs){
	  $comment .= "CDSs: @CDSs,";
	}
	if(@transcripts){
	  $comment .= "Transcripts: @transcripts,";
	}
	if(@pseudogenes){
	  $comment .= "Pseudogenes: @pseudogenes";
	}

	return($comment);
      }
    }
  }
}



##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################




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
log files are still written to $logs_dir

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
