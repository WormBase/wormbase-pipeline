#!/usr/local/bin/perl5.8.0 -w
#
# make_wormpep
# 
# Usage : make_wormpep.pl 
#
# Builds a wormpep data set from the current autoace database
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2006-02-23 11:51:14 $

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use IO::Handle;
use Ace;
use Socket;
use File::Copy;

##############################
# command-line options       #
##############################

my ($help, $debug, $test, $verbose, $store, $wormbase);

my $initial;             # run script in initial mode at start of build
my $final;               # for full run at end of build

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
            "initial"    => \$initial,
	    "final"      => \$final,
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


#########################################
# sanity checks on command-line options #
#########################################

$log->log_and_die("=> You are running initial and full at same time\n\n") if ($initial and $final);


#######################################
# misc. variables                     #
#######################################

my $tace = $wormbase->tace; 
my $release = $wormbase->get_wormbase_version;

my %peptide2number;   # peptide sequence is key, CE number (and just the number) is value
my @number2peptide;   # stores peptide sequence at array positions corresponding to CE numbers
my @number2accession; # stores full Wormpep accessions (e.g. CE04323) in array, indexed as above
my @CDSs;             # stores list of CDS object names from autoace
my %cds2number;       # stores cds name (e.g. AH6.1) as key, Wormpep number (e.g. 4323) as value
my %cds2cgc_name;     # cds name and corresponding CGC name name
my %cds2gene;         # cds name and corresponding Gene ID (e.g. WBGene00012312)
my %cds2id;           # 'Brief_identification' field for each CDS
my %cds_status;       # 'Confirmed', 'Partially_confirmed', 'Predicted' status for each CDS
my %cds2protein_id;   # protein ID for each CDS
my %cds2protein_ac;   # Trembl/Swissprot protein accession for each CDS
my %cds2protein_db;   # TREMBL/TREMBLNEW/SWISSPROT
my $wpmax = 0;        # holds highest CE number in Wormpep (old release proteins + new)
my $old_wpmax;        # holds highest CE number in Wormpep (just in old release)


##########################################
# Set up database paths                  #
##########################################

# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir   = $wormbase->basedir;
my $dbdir     = $wormbase->autoace;
my $new_wpdir = $wormbase->wormpep;

# need to get previous build WORMPEP
my ($stem, $old_release ) = $new_wpdir =~ /(.*wormpep)(\d+)/;
my $wpdir     = "$stem".--$old_release;



#####################################################################################################
#
#
#                   M  A  I  N      B  O  D  Y      O  F      S  C  R  I  P  T
#
#
#####################################################################################################


# create directory structure and read wp.fasta file from previous release
&setup;

# now query autoace to get details of latest set of CDS names, writes wormpep.dnaXXX file
# also fills in hashes used elsewhere
&write_wormpep_dna;

# grab protein information for each CDS from autoace, using Table-Maker
# in -initial mode, this is just getting CGC name and Brief_identification information
# in -final mode, also gets protein ID, and protein Database accessions, and confirmation status
&retrieve_cds_data;

# write main wormpepXXX file and wormpep.tableXXX file 
# just writes wormpepXXX file if in initial mode
&write_main_wormpep_and_table;

if ($initial) {

  # write new wp.fastaXXX file
  &write_wormpep_fasta;

  #generate file to ad new peptides to mySQL database.
  $wormbase->run_script("new_wormpep_entries.pl", $log) if $initial;

  # write the wormep.historyXXX and the wormpep.diffXXX files
  &write_wormpep_history_and_diff;
}

if ($final) {

  # get all of the other info ( PFAM, InterPro, proteinId )
  &get_additional_data;

  # count the isoforms of each CDS (stats for release letter)
  &count_isoforms;

  # write wormpep accession file
  &write_wormpep_accession;


  # update common data
  $wormbase->run_script("update_Common_data.pl --build --cds2wormpep", $log);
}

##############################
# Tidy up and finish stuff   #
##############################

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);





#################################################################################
# Subroutines                                                                   #
#################################################################################

#################################################################################################
#
# Process existing wp.fasta file to populate data structures with existing information
#
# fills in %peptide2number, @number2peptide, and @number2accession
#
#################################################################################################

sub setup{

  # Make new directory for current release         
  if (-e $new_wpdir){
    $log->write_to("$new_wpdir already exists\n");
  }
  else {
    mkdir ("$new_wpdir" , 0755) || $log->log_and_die("=> Failed to create a new directory for wormpep release wormpep$release\n\n");
    $log->write_to($wormbase->runtime.": making wormpep$release\n\n");
  }
  

  # read in the wp.fasta file, contains all protein sequences ever assigned
  $log->write_to($wormbase->runtime." : initiate parsing of wp.fasta file\n");
  print $wormbase->runtime, ": initiate parsing of wp.fasta file\n" if ($verbose);

  undef (my $id) ;
  my $peptide    = ""; # will store peptide sequence for each Wormpep protein
  my $duplicates = 0;  # for tracking duplicate protein sequences

  open (WP , "$wpdir/wp.fasta$old_release") || $log->log_and_die("=> Failed to open the old wp.fasta for wormpep release wormpep$old_release\n\n");
  while (<WP>) {
    chomp;
    # is it the FASTA header line?
    if (/^>(\S+)/) {
      my $new_id = $1; # e.g. CE04232

      # Only enter this if loop if you have already grabbed peptide sequence
      if ($id) {
	# set $number to be just the actual numerical part of $id e.g. 4232
	$id =~ /CE0*([1-9]\d*)/ ; my $number = $1; 
        # bewarned - lurking 100000 protein bug - uses this pattern to cope with numbers beginning with 0.
	
	# raise maximum number?
	($wpmax = $number) if ($number > $wpmax);

	# add peptide sequnce into array, CE number is index in array
	$number2peptide[$number] = $peptide;

	# If this is a new peptide sequence in hash...
	unless (exists ($peptide2number{$peptide})) { 
	  # add to hash: where peptide sequence is key, and CE number is value
	  $peptide2number{$peptide} = $number;
	  # also store CE accession in array, indexed by number part of accession
	  $number2accession[$number] = $id;
	} 
	# else treat as duplicate protein sequence
	else {
	  ++$duplicates;
	  my $old_id = $peptide2number{$peptide};
	  my $old_id_pad = sprintf "%05d" , $old_id;
	  # append info about old accession number
	  $number2accession[$number] = "$id -> CE$old_id_pad";
	}
      }
      $id = $new_id ; $peptide = "" ;
    }
    # are you at end of file?
    elsif (eof) {
      if ($id) {
	$id =~ /CE0*([1-9]\d*)/ ; my $number = $1 ;
	($wpmax = $number) if ($number > $wpmax);

	$peptide .= $_ ;
	$number2peptide[$number] = $peptide;
	unless (exists ($peptide2number{$peptide})) { 
	  $peptide2number{$peptide} = $number;
	  $number2accession[$number] = $id;
	} else {
	  ++$duplicates;
	  my $old_id = $peptide2number{$peptide};
	  my $old_id_pad = sprintf "%05d" , $old_id;
	  $number2accession[$number] = "$id -> CE$old_id_pad";
	}
      }
    } 
    # read a line of amino acid sequence, add to $peptide
    else {
      $peptide .= $_ ;
    }
  }
  
  # set $old_wpmax to be highest CE number from release x-1
  # $wpmax can then be incremented with info from any new gene predictions
  $old_wpmax = $wpmax;

  $log->write_to("=> wp.fasta file contains $duplicates duplicate sequences\n");
  print     "=> wp.fasta file contains $duplicates duplicate sequences\n" if ($verbose);
  $log->write_to("=> wpmax of wp.fasta$old_release equals $old_wpmax\n");
  print     "=> wpmax of wp.fasta$old_release equals $old_wpmax\n" if ($verbose);
  $log->write_to($wormbase->runtime." : completed parsing of wp.fasta file\n\n");
  print     $wormbase->runtime, ": completed parsing of wp.fasta file\n\n" if ($verbose);
  close (WP);

}


##########################################################################################
#
# write_wormpep_dna
# 
# 1) retrieves list of valid ?CDS objects
# 2) gets dna and peptide sequence for each ?CDS, writes dna.fasta file,
# 3) maps the peptide sequences onto wormpep, deleting Peptides containing an X or a *,
# 4) creates new Wormpep IDs if necessary, creating %dot2num
#
##########################################################################################

sub write_wormpep_dna{

  $log->write_to($wormbase->runtime." : connecting to $dbdir\n");
  print     $wormbase->runtime, ": connecting to $dbdir\n" if ($verbose);

  # grab list of valid CDS names, ignore any temp genes
  my $db = Ace->connect (-path => $dbdir, -program => $tace) || $log->log_and_die("=> Failed to connect to primary database $dbdir\n\n");
  @CDSs = $db->fetch (-query => 'FIND elegans_CDS NOT *temp*');
  @CDSs = sort @CDSs;

  $log->write_to("=> ".scalar(@CDSs)." CDSs\n");
  print     "=> " . scalar(@CDSs) . " CDSs\n" if ($verbose);


  # write DNA for each CDS to separate file
  $log->write_to($wormbase->runtime.": creating wormpep.dna file\n\n");
  open (DNA , ">$new_wpdir/wormpep.dna$release") || $log->log_and_die("=> Failed to create a new wormpep.dna for wormpep release wormpep$release\n\n");


  foreach my $cds (@CDSs) {
    
    # get dna 
    my $dna = $cds->asDNA();
    $log->write_to("cannot extract dna for CDS $cds\n") if ((!defined ($dna)) || ($dna eq ""));
    $dna =~ /^\n>(\S+)\s+(\w.*)/s ; my $dna_seq = $2 ; $dna_seq =~ tr/a-z/A-Z/ ; $dna_seq =~ s/\s//g;
    if ($dna_seq =~ /[^ACGT]/) {
      if ($dna_seq =~ /\-/) {                                       # - seems to indicate that e.g the subsequence
	$log->write_to("ERROR: $cds - DNA sequence contains a -\n"); # coordinates differ from the last exon coordinate
      } 
      elsif ($dna_seq =~ /N/) {
	$log->write_to("ERROR: $cds - DNA sequence contains an N\n");
      } 
      else {                         
	$log->write_to("ERROR: $cds - DNA sequence contains a non-ACGT character (which isn't a '-' or 'N')\n");
      }
    }
    print DNA "$dna";

    # grab peptide sequence
    my $peptide = $cds->asPeptide();
    if ((!defined ($peptide)) || ($peptide eq "")) {
      $log->write_to("cannot extract peptide sequence for CDS $cds\n");
      next;
    }
    $peptide =~ /^>(\S+)\s+([\w\*].*)/s ; 
    my $peptide_seq = $2 ; $peptide_seq =~ tr/a-z/A-Z/ ; $peptide_seq =~ s/\s//g;


    # various santity checks
    $log->write_to("ERROR: $cds has an incorrect dotname\n") unless ($cds =~ /^[A-Z0-9_]+\.[1-9]\d?[A-Za-z]?$/i);    
    $log->write_to("ERROR: $cds - peptide sequence does not start with a M\n") unless ($peptide_seq =~ /^M/);
    $log->write_to("ERROR: $cds - peptide sequence contains an X\n")           if ($peptide_seq =~ /X/);
    if ($peptide_seq =~ /\*/) {
      $log->write_to("ERROR: $cds - peptide sequence contains a *\n");
      next;
    }
    
    # check current peptide sequence against peptides loaded in from last releases wp.fasta file,
    # if peptide is new add to hash
    unless (exists ($peptide2number{$peptide_seq})){
      $number2peptide[++$wpmax] = $peptide_seq; # adds new sequence to array, increases $wpmax
      $peptide2number{$peptide_seq} = $wpmax;
      $cds2number{$cds} = $wpmax;
      my $pad = sprintf "%05d" , $wpmax;
      $number2accession[$wpmax] = "CE$pad\t$cds";
    } 
    else {
      $cds2number{$cds} = $peptide2number{$peptide_seq};
      my $number = $peptide2number{$peptide_seq};
      $number2accession[$number] .= "\t$cds";
    }             
    
    # tidy up
    $cds->DESTROY();
    
  }   
  close DNA;
  $log->write_to($wormbase->runtime." : finished writing wormpep.dna file\n\n");

  # close database connection
  $db->close;
  $log->write_to($wormbase->runtime." : finished connection to database\n\n");
  print     $wormbase->runtime, ": finished connection to database\n\n" if ($verbose);

}


###############################################################################
# retrieve various protein data from autoace using a tablemaker call          #
# in -intitial mode there will not be much data you can get                   #
###############################################################################

sub retrieve_cds_data{

  $log->write_to($wormbase->runtime." : retrieving data from autoace for each CDS\n\n");
  $log->log_and_die("Tablemaker def file wormpep.def doesn't exist\n") unless -e ("$dbdir/wquery/wormpep.def");
  my $command = "Table-maker -p $dbdir/wquery/wormpep.def\nquit\n";
  
  open (TACE, "echo '$command' | $tace $dbdir | ");
  while (<TACE>) {
    print if ($verbose && $final);
    chomp;
    s/acedb\> //g;      # only need this is using 4_9i code, bug fixed in 4_9k onward (should be redundant)
    next if ($_ eq "");
    next if (/\/\//);
    s/\"//g;
    (/^(\S+)\s/);
    my ($cds,$prot_id_parent,$prot_id,$prot_id_ver,$prot_db,$prot_ac,$gene,$cgc_name,$confirmed,$partial,$brief_id) = split /\t/;

    # load hashes with data
    ($cds2cgc_name{$cds} = $cgc_name) if (defined($cgc_name));
    ($cds2id{$cds} = $brief_id)       if (defined($brief_id));
    ($cds2gene{$cds} = $gene)         if (defined($gene));

    # can only work with some of this data in full wormpep mode, not in initial
    if($final){
      $cds2protein_ac{$cds} = $prot_ac;
      $cds2protein_id{$cds} = $prot_id . "." . $prot_id_ver;
      $cds2protein_id{$cds} = "" if $cds2protein_id{$cds} eq ".";
      $cds2protein_db{$cds} = $prot_db;
      
      # confirmed CDS - [EST|mRNA] data
      if ($confirmed) {
	$cds_status{$cds} = "Confirmed";
      }
      # supported CDS - [EST] data
      elsif ($partial) {
	$cds_status{$cds} = "Partially_confirmed";
      }
      # predicted CDS - no data
      else {
	$cds_status{$cds} = "Predicted";
      }
    }
    print "$gene $cds Locus:$cds2cgc_name{$cds} Protein_ID:$cds2protein_id{$cds} Protein_AC:$cds2protein_ac{$cds} Protein_DB:$cds2protein_db{$cds} Brief_ID:$cds2id{$cds}\n\n" if ($verbose);
    
  }
  close(TACE);
  
  $log->write_to($wormbase->runtime." : finished retrieving data from autoace for each CDS\n\n");
}



############################################################################
# write the new wp.fasta file, and wrap the sequences using rd's seqpress  #
############################################################################

sub write_wormpep_fasta{

  $log->write_to($wormbase->runtime." : writing wp.fasta$release file\n\n");

  open (WPFASTA , ">$new_wpdir/wp.fasta_unwrap$release") or $log->log_and_die("cant open $new_wpdir/wp.fasta_unwrap$release :$!\n");

  # loop through array of peptide sequences
  for (my $i = 1 ; $i <= $wpmax ; $i++) {
    if (defined ($number2peptide[$i])) {
      my $i_pad = sprintf "%05d" , $i;
      print WPFASTA ">CE$i_pad\n$number2peptide[$i]\n";
    }
  }
  close(WPFASTA);

  my $status = system ("rewrap $new_wpdir/wp.fasta_unwrap$release > $new_wpdir/wp.fasta$release");
  if(($status >>8) != 0){
    $log->write_to("ERROR: rewrap command failed. \$\? = $status\n");
  }
  unlink ("$new_wpdir/wp.fasta_unwrap$release") || $log->write_to("cannot delete $new_wpdir/wp.fasta_unwrap$release\n");

  $log->write_to($wormbase->runtime." : finished writing wormpep.fasta file\n\n");

}


##########################################################################################
# get from autoace the data required to write the wormpep file and the wormpep.table file
# use rd's seqpress to wrap the sequence lines in the wormpep file 
##########################################################################################

sub write_main_wormpep_and_table{

  $log->write_to($wormbase->runtime." : Build wormpep & wormpep.table files\n\n");
  
  open (FASTA , ">$new_wpdir/wormpep_unwrap$release") || $log->log_and_die("cannot create wormpep_unwrap$release\n");

  FASTA->autoflush();

  if ($initial) {
    my $ace_dir         = $wormbase->autoace;     # AUTOACE DATABASE DIR

    open (CONNECTIONS, ">".$wormbase->acefiles."/CDS2wormpep.ace") ||  $log->log_and_die("cannot create CDS2wormpep\n");
    CONNECTIONS->autoflush();
  } elsif ($final) {
    open (TABLE, ">".$wormbase->wormpep."/wormpep.table$release")  || $log->log_and_die("cannot create wormpep.table$release\n");
    TABLE->autoflush();
  }
  
  foreach my $cds (@CDSs) {
   
      # reset all fields
      
      my $wpid     = $cds2number{$cds};
      my $wpid_pad = sprintf "%05d" , $wpid;
      my $pepseq   = $number2peptide[$wpid];
      
      # set the fields to be printed depending on whether they exist

      my $output = ">$cds CE$wpid_pad";
      ($output .= " $cds2gene{$cds}")                    if ($cds2gene{$cds});
      ($output .= " locus\:".$cds2cgc_name{$cds})        if ($cds2cgc_name{$cds});
      ($output .= " $cds2id{$cds}")                      if ($cds2id{$cds});
      ($output .= " status\:".$cds_status{$cds})         if ($cds_status{$cds});
      
      # can only get the following information if running in final mode, else won't yet be in autoace
      if ($final) {
	  ($output .= " SW:$cds2protein_ac{$cds}")           if (($cds2protein_db{$cds} eq "SwissProt") && defined($cds2protein_ac{$cds}));
	  ($output .= " TR:$cds2protein_ac{$cds}")           if (($cds2protein_db{$cds} eq "TREMBL")    && defined($cds2protein_ac{$cds}));
	  ($output .= " TN:$cds2protein_ac{$cds}")           if (($cds2protein_db{$cds} eq "TREMBLNEW") && defined($cds2protein_ac{$cds}));
	  ($output .= " protein_id\:".$cds2protein_id{$cds}) if ($cds2protein_id{$cds});
      }
      $output .= "\n$pepseq\n";
      print FASTA "$output";

      # print out connections in -initial mode
      
      if ($initial) {

	  $output  = "CDS : \"$cds\"\n";
	  $output .= "Corresponding_protein WP:CE${wpid_pad}\n\n";
	  print CONNECTIONS "$output";

      }
      
      # print to table in -final mode
    
      if ($final) {
	  my $output = ">$cds\tCE$wpid_pad";
	  $output   .= "\t";
	  ($output  .= "$cds2cgc_name{$cds}")       if ($cds2cgc_name{$cds});
	  $output   .= "\t";
	  ($output  .= "$cds2id{$cds}")             if ($cds2id{$cds});
	  $output   .= "\t";
	  ($output  .= "$cds_status{$cds}")         if ($cds_status{$cds});
	  $output .= "\t";
	  ($output .= "SW:$cds2protein_ac{$cds}")  if (($cds2protein_db{$cds} eq "SwissProt") && ($cds2protein_ac{$cds}));
	  ($output .= "TR:$cds2protein_ac{$cds}")  if (($cds2protein_db{$cds} eq "TREMBL")    && ($cds2protein_ac{$cds}));
	  ($output .= "TN:$cds2protein_ac{$cds}")  if (($cds2protein_db{$cds} eq "TREMBLNEW") && ($cds2protein_ac{$cds}));
	  $output .= "\t";
	  ($output .= "$cds2protein_id{$cds}")     if ($cds2protein_id{$cds});
	  print TABLE "$output\n";
      } 
  } 
  
  close (FASTA);
  close (CONNECTIONS) if ($initial);
  close (TABLE) if ($final);
  
  my $status = system ("rewrap $new_wpdir/wormpep_unwrap$release > $new_wpdir/wormpep$release");
  if(($status >>8) != 0){
    $log->write_to("ERROR: rewrap command failed. \$\? = $status\n");
  }

  # create a blast'able database (indexing) using setdb for Wublast (not formatdb, which is  for blastall)
  if($final){
    $status = copy("$new_wpdir/wormpep$release", "$new_wpdir/wormpep_current");
    $log->write_to("ERROR: Couldn't copy file: $!\n") if ($status == 0);
    $status = system ("/usr/local/pubseq/bin/setdb $new_wpdir/wormpep_current > $new_wpdir/wormpep_current.log");
    if(($status >>8) != 0){
      $log->write_to("ERROR: setdb command failed. \$\? = $status\n");
    }
  }

  unlink ("$new_wpdir/wormpep_unwrap$release")  || $log->write_to("cannot delete $new_wpdir/wormpep_unwrap$release\n");
  $log->write_to($wormbase->runtime." : Finished building wormpep & wormpep.table files\n\n");

}


############################################################################
# read in the current wormpep.history file, update it, and read it back out,
# wormpep.history contains tab-separated:  dotname, wpid, start (release), end (release)

sub write_wormpep_history_and_diff{
  open (OLDHISTORY, "$wpdir/wormpep.history$old_release")  || $log->log_and_die("cannot open $wpdir/wormpep.history$old_release\n");
  open (HISTORY,    ">$new_wpdir/wormpep.history$release") || $log->log_and_die("cannot create $wpdir/wormpep.history$release\n");
  open (DIFF,       ">$new_wpdir/wormpep.diff$release")    || $log->log_and_die("cannot create $wpdir/wormpep.diff$release\n");
  
  my %line;
  while (my $line = <OLDHISTORY>) {
    chomp $line;
    my $cds = ""; my $wpid = ""; my $start = ""; my $end = "";
    ($cds , $wpid , $start , $end) = split (/\t/ , $line);
    $wpid =~ /CE0*([1-9]\d*)/ ; my $num = $1;
    $line{$cds} = $line;

    if (!exists($cds2number{$cds})){
      if ($end eq "") {
	print HISTORY "$cds\t$wpid\t$start\t$release\n";
	print DIFF "lost:\t$cds\t$wpid\n";
      } 
      else{
	print HISTORY "$cds\t$wpid\t$start\t$end\n";
      }
    }
    elsif (($cds2number{$cds} ne $num) && ($end eq "")) {
      print HISTORY "$cds\t$wpid\t$start\t$release\n";
      my $new_num = $cds2number{$cds};
      my $new_pad = sprintf "%05d" , $new_num;
      print HISTORY "$cds\tCE$new_pad\t$release\t$end\n";
      print DIFF "changed:\t$cds\t$wpid --> CE$new_pad\n";
    } 
    else {
      print HISTORY "$cds\t$wpid\t$start\t$end\n";
    }
  }
  
  foreach my $cds (keys (%cds2number)) {
    my $empty = "";
    if (!exists ($line{$cds})) {
      my $num = $cds2number{$cds};
      my $pad = sprintf "%05d" , $num;
      print HISTORY "$cds\tCE$pad\t$release\t$empty\n";
      print DIFF "new:\t$cds\tCE$pad\n";
      next;
    }
    my $cdsname = ""; my $wpid = ""; my $start = ""; my $end = "";
    ($cdsname , $wpid , $start , $end) = split (/\t/ , $line{$cds});
    if ($end ne "") {   
      my $new_num = $cds2number{$cds};
      my $new_pad = sprintf "%05d" , $new_num;
      print HISTORY "$cds\tCE$new_pad\t$release\t$empty\n";
      print DIFF "reappeared:\t$cds\tCE$new_pad\n";
    }
  }
  
  close OLDHISTORY;
  close HISTORY;
  close DIFF;
  my $wpdiff = $wpmax - $old_wpmax;
  $log->write_to("\n\nnew wpmax of wp.fasta$release equals $wpmax\n$wpdiff new sequences have been added\n");
}



##########################################################################################
# count the CDS (with and without alternate splice forms) based on the wormpep file
##########################################################################################

sub count_isoforms{

  my @wormpep_CDSs;
  open (FILE , "$new_wpdir/wormpep$release") || $log->log_and_die("cannot open the wormpep$release file\n");
  while (<FILE>) {
    if (/^\>(\S+)\s+\S+/) {
      my $cds = $1;
      push (@wormpep_CDSs , $cds);
    }
  }
  close (FILE);
  
  my $total_cds_count  = 0;
  my $no_isoform_count = 0;
  my $isoform_count    = 0;
  my %new_name2x;
  
  foreach (@wormpep_CDSs) {
    $total_cds_count++;
    /^([A-Z0-9_]+)\.(.*)$/i;
    my $cds_prefix = $1;
    my $cds_suffix = $2;
    if ($cds_suffix =~ /^[1-9][0-9]?$/) {
      $no_isoform_count++;
    }
    elsif ($cds_suffix =~ /(^[1-9][0-9]?)([a-z])/) {
      my $number = $1;
      my $letter = $2;
      my $new_name = $cds_prefix."_".$number;
      unless (exists ($new_name2x{$new_name})) {
	$no_isoform_count++;
      } 
      else {
	$isoform_count++;
      }
      $new_name2x{$new_name} = "x";
    } 
    else {
      $log->write_to("$_ has a non\-acceptable name in wormpep$release \(has not been counted\)\n");
      next;
    }
  }
  $log->write_to("\n\nthere are $no_isoform_count CDS in autoace, $total_cds_count when counting \($isoform_count\) alternate splice_forms\n");

  # write the release letter (also does some checks)
  if($final) {
    $wormbase->release_wormpep($no_isoform_count,$total_cds_count,$isoform_count);
    chmod (0444 , "$new_wpdir/*") || $log->write_to("cannot chmod $new_wpdir/ files\n");
  }
  
}

###########################################################################
# create the new wormpep.accession,
# wp.table contains tab-separated:  wpid, list of associated CDSs if active, 
# empty if inactive, and link to active wpid in case of duplicate sequences in wp.fasta

sub write_wormpep_accession{
  open (WPTABLE , ">$new_wpdir/wormpep.accession$release") ||$log->log_and_die("cannot create wormpep.accession$release\n");
  
  for (my $i = 1 ; $i <= $wpmax ; $i++) {
    if (defined ($number2accession[$i])) {
      print WPTABLE "$number2accession[$i]\n";
    }
  }   
  close (WPTABLE);
  
}

sub get_additional_data {
  # get Pfam domains (this step loads resulting ace file)
  $wormbase->run_script("GetPFAM_motifs.pl -load");#

  # get interpro domains (this step loads resulting ace file )
  $wormbase->run_script("GetInterPro_motifs.pl -load");

  # make interpro2go connections (to be used by getProteinID)
  $wormbase->run_script("make_Interpro2GO_mapping.pl");

  # Get protein IDs (this step writes to ~wormpub/analysis/SWALL and loads wormpep info)
  $wormbase->run_script("getProteinID.pl -load");
}


__END__

=pod

=head2   NAME - make_wormpep


=head1 USAGE

=over 4

=item make_wormpep [-options]

=back

make_wormpep will generate a peptide data set from the autoace
database directory.

autoace_minder mandatory arguments:

=over 4

=item -i or -f

=back

autoace_minder OPTIONAL arguments:

=over 4

=item -h, Help page

=item -d, Verbose/Debug mode

=item -i ,initial
just writes a FASTA style file for use in the BLAST analysis wormpep88.

=item -f , full.
creates the version of Wormpep

=item -t, test
Writes to a test directory rather than the actual release one.

=back

=head1 EXAMPLES:

=over 4

=item make_wormpep -f

=back

Creates a new wormpep data set in the (new) $basedir/WORMPEP/wormpep40 directory

=cut




















