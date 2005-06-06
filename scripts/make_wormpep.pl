#!/usr/local/bin/perl5.8.0 -w
#
# make_wormpep
# 
# Usage : make_wormpep.pl 
#
# Builds a wormpep data set from the current autoace database
#
# Last updated by: $Author: dl1 $
# Last updated on: $Date: 2005-06-06 10:51:57 $


use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use IO::Handle;
use Ace;
use Socket;
use File::Copy;


##############################
# command-line options       #
##############################

my $help;                # Help/Usage page
my $verbose;             # turn on extra output
my $debug;               # For sending output to just one person
my $maintainers = "All"; # log file recipients
my $initial;             # run script in initial mode at start of build
my $final;               # for full run at end of build
my $test;                # for running in test environment ~wormpub/TEST_BUILD

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "verbose"     => \$verbose,
            "initial"     => \$initial,
	    "final"       => \$final,
            "test"        => \$test
           );

#########################################
# sanity checks on command-line options #
#########################################

&error(0) if ($help);
&error(1) if ($initial and $final);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}


#######################################
# misc. variables                     #
#######################################

my $tace = &tace; 
my $release; 
if($test){
  $release = "666";
}
else{
  $release = &get_wormbase_version; 
}
my $old_release = $release-1;


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
my $log;              # for log file


##########################################
# Set up database paths                  #
##########################################

# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir   = "/wormsrv2";
$basedir      = glob("~wormpub")."/TEST_BUILD" if ($test); 

my $dbdir     = "$basedir/autoace";
my $wpdir     = "$basedir/WORMPEP/wormpep$old_release";
my $new_wpdir = "$basedir/WORMPEP/wormpep$release";



#####################################################################################################
#
#
#                   M  A  I  N      B  O  D  Y      O  F      S  C  R  I  P  T
#
#
#####################################################################################################


&create_log_files;


# create directory structure and read wp.fasta file from previous release
&setup;


# now query autoace to get details of latest set of CDS names, writes wormpep.dnaXXX file
# also fills in hashes used elsewhere
&write_wormpep_dna;


# write new wp.fastaXXX file
&write_wormpep_fasta;


# grab protein information for each CDS from autoace, using Table-Maker
# in -initial mode, this is just getting CGC name and Brief_identification information
# in -final mode, also gets protein ID, and protein Database accessions, and confirmation status
&retrieve_cds_data;


# write main wormpepXXX file and wormpep.tableXXX file 
# just writes wormpepXXX file if in initial mode
&write_main_wormpep_and_table;


# write the wormep.historyXXX and the wormpep.diffXXX files
&write_wormpep_history_and_diff;


# count the isoforms of each CDS (stats for release letter)
&count_isoforms if ($final);

# write wormpep accession file
&write_wormpep_accession if ($final);


##############################
# Tidy up and finish stuff   #
##############################

close (LOG);

#Email log
my $subject = "BUILD REPORT: make_wormpep.pl";
$subject = "TEST MODE: WormBase Build Report: make_wormpep.pl" if ($test);

&mail_maintainer("$subject",$maintainers,$log);

exit(0);





#################################################################################
# Subroutines                                                                   #
#################################################################################

##########################
# run details            #
##########################

sub run_details {
  print "# make_wormpep\n";
  print "# run details    : ",&rundate," ",&runtime,"\n";
  print "\n";
  print "Wormpep version  : wormpep$release\n";
  print "Primary database : $dbdir\n\n";
  
} # end of sub 'run details'

##########################
# errors from the script #
##########################

sub error {
  my $error = shift;
  # Error 0 - help page
  if ($error == 0) {
    exec ('perldoc',$0);
    exit (0);
  }
  #   # Error  1 - initial and full at same time - WARN only
  elsif ($error == 1) {
    &run_details;
    print "=> Running initial and full at same time\n\n";
    print "are you sure you want to continue ( y / n )\n";
    my $yn = <STDIN>;
    if ("$yn" ne "y" ) {
      return;
    }
    else {
      print "Exiting\n";
      exit(0);
    }
  }
  # Error  2 - not chosen initial or full
  elsif ($error == 2) {
    &run_details;
    print "=> not chosen initial or full.\n\n";
    exit(0);
  }
  # Error  3 - cannot create new wormpep directory 
  elsif ($error == 3) {
    &run_details;
    print "=> Failed to create a new directory for wormpep release wormpep$release\n\n";
    print LOG "=> Failed to create a new directory for wormpep release wormpep$release\n\n";
    print LOG &runtime, ": Exiting script\n";
    close LOG;
    &mail_maintainer("WormBase Report: make_wormpep",$maintainers,$log);
  }
  # Error  4 - cannot open new wp.log file 
  elsif ($error == 4) {
    &run_details;
    print "=> Failed to create a new wp.log for wormpep release wormpep$release\n\n";
    print LOG &runtime, ": Exiting script\n";
    close LOG;
    &mail_maintainer("WormBase Report: make_wormpep",$maintainers,$log);
  }
  # Error  5 - cannot open old wp.fasta file 
  elsif ($error == 5) {
    &run_details;
    print "=> Failed to open the old wp.fasta for wormpep release wormpep$old_release\n\n";
    print LOG &runtime, ": Exiting script\n";
    close LOG;
    &mail_maintainer("WormBase Report: make_wormpep",$maintainers,$log);
  }
  # Error  6 - cannot connect to ACEDB database 
  elsif ($error == 6) {
    &run_details;
    print "=> Failed to connect to primary database 'dbdir'\n\n";
    print LOG &runtime, ": Exiting script\n";
    close LOG;
    &mail_maintainer("WormBase Report: make_wormpep",$maintainers,$log);
  }
  # Error  7 - cannot open new wp.log file 
  elsif ($error == 7) {
    &run_details;
    print "=> Failed to create a new wormpep.dna for wormpep release wormpep$release\n\n";
    print LOG &runtime, ": Exiting script\n";
    close LOG;
    &mail_maintainer("WormBase Report: make_wormpep",$maintainers,$log);
  }
  
  exit(1);
} # end of sub 'error'



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
    print LOG "$new_wpdir already exists\n";
  }
  else {
    mkdir ("$new_wpdir" , 0755) || &error(3);               # die "cannot create the $new_wpdir directory\n";
    print LOG &runtime, ": making wormpep$release\n\n";
  }
  

  # read in the wp.fasta file, contains all protein sequences ever assigned
  print LOG &runtime, ": initiate parsing of wp.fasta file\n";
  print     &runtime, ": initiate parsing of wp.fasta file\n" if ($verbose);

  undef (my $id) ;
  my $peptide    = ""; # will store peptide sequence for each Wormpep protein
  my $duplicates = 0;  # for tracking duplicate protein sequences

  open (WP , "$wpdir/wp.fasta$old_release") || &error(5); # die "couldn't open $wpdir/wp.fasta$old_release\n";
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

  print LOG "=> wp.fasta file contains $duplicates duplicate sequences\n";
  print     "=> wp.fasta file contains $duplicates duplicate sequences\n" if ($verbose);
  print LOG "=> wpmax of wp.fasta$old_release equals $old_wpmax\n";
  print     "=> wpmax of wp.fasta$old_release equals $old_wpmax\n" if ($verbose);
  print LOG &runtime, ": completed parsing of wp.fasta file\n\n";
  print     &runtime, ": completed parsing of wp.fasta file\n\n" if ($verbose);
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

  print LOG &runtime, ": connecting to $dbdir\n";
  print     &runtime, ": connecting to $dbdir\n" if ($verbose);

  # grab list of valid CDS names, ignore any temp genes
  my $db = Ace->connect (-path => $dbdir, -program => $tace) || &error(6); # die "cannot connect to autoace\n";
  @CDSs = $db->fetch (-query => 'FIND elegans_CDS NOT *temp*');
  @CDSs = sort @CDSs;

  print LOG "=> " . scalar(@CDSs) . " CDSs\n";
  print     "=> " . scalar(@CDSs) . " CDSs\n" if ($verbose);


  # write DNA for each CDS to separate file
  print LOG &runtime, ": creating wormpep.dna file\n\n";
  open (DNA , ">$new_wpdir/wormpep.dna$release") || &error(7); # die "cannot create $new_wpdir/wormpep.dna$release\n";


  foreach my $cds (@CDSs) {
    
    # get dna 
    my $dna = $cds->asDNA();
    print LOG "cannot extract dna for CDS $cds\n" if ((!defined ($dna)) || ($dna eq ""));
    $dna =~ /^>(\S+)\s+(\w.*)/s ; my $dna_seq = $2 ; $dna_seq =~ tr/a-z/A-Z/ ; $dna_seq =~ s/\s//g;
    if ($dna_seq =~ /[^ACGT]/) {
      if ($dna_seq =~ /\-/) {                                       # - seems to indicate that e.g the subsequence
	print LOG "ERROR: $cds - DNA sequence contains a -\n"; # coordinates differ from the last exon coordinate
      } 
      elsif ($dna_seq =~ /N/) {
	print LOG "ERROR: $cds - DNA sequence contains an N\n";
      } 
      else {                         
	print LOG "ERROR: $cds - DNA sequence contains a non-ACGT character (which isn't a '-' or 'N')\n";  
      }
    }
    print DNA "$dna";

    # grab peptide sequence
    my $peptide = $cds->asPeptide();
    if ((!defined ($peptide)) || ($peptide eq "")) {
      print LOG "cannot extract peptide sequence for CDS $cds\n";
      next;
    }
    $peptide =~ /^>(\S+)\s+([\w\*].*)/s ; 
    my $peptide_seq = $2 ; $peptide_seq =~ tr/a-z/A-Z/ ; $peptide_seq =~ s/\s//g;


    # various santity checks
    print LOG "ERROR: $cds has an incorrect dotname\n" unless ($cds =~ /^[A-Z0-9_]+\.[1-9]\d?[A-Za-z]?$/i);    
    print LOG "ERROR: $cds - peptide sequence does not start with a M\n" unless ($peptide_seq =~ /^M/);
    print LOG "ERROR: $cds - peptide sequence contains an X\n"           if ($peptide_seq =~ /X/);
    if ($peptide_seq =~ /\*/) {
      print LOG "ERROR: $cds - peptide sequence contains a *\n";              
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
  print LOG &runtime, ": finished writing wormpep.dna file\n\n";

  # close database connection
  $db->close;
  print LOG &runtime, ": finished connection to database\n\n";
  print     &runtime, ": finished connection to database\n\n" if ($verbose);

}


###############################################################################
# retrieve various protein data from autoace using a tablemaker call          #
# in -intitial mode there will not be much data you can get                   #
###############################################################################

sub retrieve_cds_data{

  print LOG &runtime, ": retrieving data from autoace for each CDS\n\n";
  
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
  
  print LOG &runtime, ": finished retrieving data from autoace for each CDS\n\n";
}



############################################################################
# write the new wp.fasta file, and wrap the sequences using rd's seqpress  #
############################################################################

sub write_wormpep_fasta{

  print LOG &runtime, ": writing wp.fasta$release file\n\n";

  open (WPFASTA , ">$new_wpdir/wp.fasta_unwrap$release");

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
    print LOG "ERROR: rewrap command failed. \$\? = $status\n";
  }
  unlink ("$new_wpdir/wp.fasta_unwrap$release") || print LOG "cannot delete $new_wpdir/wp.fasta_unwrap$release\n";

  print LOG &runtime, ": finished writing wormpep.fasta file\n\n";

}


##########################################################################################
# get from autoace the data required to write the wormpep file and the wormpep.table file
# use rd's seqpress to wrap the sequence lines in the wormpep file 
##########################################################################################

sub write_main_wormpep_and_table{

  print LOG &runtime, ": Build wormpep & wormpep.table files\n\n";
  
  # open filehandles for output files
  
  open (FASTA, ">$new_wpdir/wormpep_unwrap$release") || die "cannot create wormpep_unwrap$release\n";
  FASTA->autoflush();

  if ($initial) {
      open (CONNECTIONS, ">/wormsrv2/autoace/acefiles/CDS2wormpep.ace") || die "cannot create CDS2wormpep\n";
      CONNECTIONS->autoflush();
  }
  elsif ($final) {
      open (TABLE, ">$new_wpdir/wormpep.table$release")  || die "cannot create wormpep.table$release\n";
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
    print LOG "ERROR: rewrap command failed. \$\? = $status\n";
  }

  # create a blast'able database (indexing) using setdb for Wublast (not formatdb, which is  for blastall)
  if($final){
    $status = copy("$new_wpdir/wormpep$release", "$new_wpdir/wormpep_current");
    print LOG "ERROR: Couldn't copy file: $!\n" if ($status == 0);
    $status = system ("/usr/local/pubseq/bin/setdb $new_wpdir/wormpep_current > $new_wpdir/wormpep_current.log");
    if(($status >>8) != 0){
      print LOG "ERROR: setdb command failed. \$\? = $status\n";
    }
  }

  unlink ("$new_wpdir/wormpep_unwrap$release")  || print LOG "cannot delete $new_wpdir/wormpep_unwrap$release\n";
  print LOG &runtime, ": Finished building wormpep & wormpep.table files\n\n";

}


############################################################################
# read in the current wormpep.history file, update it, and read it back out,
# wormpep.history contains tab-separated:  dotname, wpid, start (release), end (release)

sub write_wormpep_history_and_diff{
  open (OLDHISTORY, "$wpdir/wormpep.history$old_release")  || die "cannot open $wpdir/wormpep.history$old_release\n";
  open (HISTORY,    ">$new_wpdir/wormpep.history$release") || die "cannot create $wpdir/wormpep.history$release\n";
  open (DIFF,       ">$new_wpdir/wormpep.diff$release")    || die "cannot create $wpdir/wormpep.diff$release\n";
  
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
  print LOG "\n\nnew wpmax of wp.fasta$release equals $wpmax\n$wpdiff new sequences have been added\n";
}



##########################################################################################
# count the CDS (with and without alternate splice forms) based on the wormpep file
##########################################################################################

sub count_isoforms{

  my @wormpep_CDSs;
  open (FILE , "$new_wpdir/wormpep$release") || print LOG  "cannot open the wormpep$release file\n";
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
      print LOG "$_ has a non\-acceptable name in wormpep$release \(has not been counted\)\n";
      next;
    }
  }
  print LOG "\n\nthere are $no_isoform_count CDS in autoace, $total_cds_count when counting \($isoform_count\) alternate splice_forms\n";

  # write the release letter (also does some checks)
  if($final) {
    &release_wormpep($no_isoform_count,$total_cds_count,$isoform_count);
    chmod (0444 , "$new_wpdir/*") || print LOG "cannot chmod $new_wpdir/ files\n";
  }
  
}

###########################################################################
# create the new wormpep.accession,
# wp.table contains tab-separated:  wpid, list of associated CDSs if active, 
# empty if inactive, and link to active wpid in case of duplicate sequences in wp.fasta

sub write_wormpep_accession{
  open (WPTABLE , ">$new_wpdir/wormpep.accession$release") || die "cannot create wormpep.accession$release\n";
  
  for (my $i = 1 ; $i <= $wpmax ; $i++) {
    if (defined ($number2accession[$i])) {
      print WPTABLE "$number2accession[$i]\n";
    }
  }   
  close (WPTABLE);
  
}


#################################################################################################
#
# Create logfile
#
#################################################################################################

sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch $basedir/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate = &rundate;
  # Log name depends on which mode script is being run in: initial/final
  $log        = "$basedir/logs/$script_name.initial.WS$release.$rundate.$$" if ($initial);
  $log        = "$basedir/logs/$script_name.final.WS$release.$rundate.$$"   if ($final);

  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name\n";
  print LOG "Wormpep version  : wormpep$release\n\n";
  print LOG "initial phase\n" if ( $initial );
  print LOG "full phase\n"    if ( $final );
  print LOG "started at ",&runtime,"\n";
  print LOG "=============================================\n";
  print LOG "\n";

  LOG->autoflush();
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




















