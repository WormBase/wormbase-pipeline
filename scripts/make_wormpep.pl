#!/usr/local/bin/perl5.8.0 -w
#
# make_wormpep
# 
# Usage : make_wormpep.pl 
#
# Builds a wormpep data set from the current autoace database
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2003-12-02 09:32:07 $

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use IO::Handle;
use Ace;
use Socket;


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
            "intitial"    => \$initial,
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

my $tace      = &tace; 
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
my %cds2num;          #
my @dotnames;         #
my @c_dotnames;       #
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
$ENV{'ACEDB'} = $dbdir;





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

# now query autoace to get details of latest set of CDSs
&process_cds_class;

# files in wormpep directory
#   wp.fasta            all sequences ever assigned (wpid and sequence)
#   wormpep.table       dotname (primary), wpid, locus, brief_id, lab origin
#   wormpep             current proteins with info from wormpep.table
#   wormpep.accession   wpid, if active all associated dotnames, if duplicate link to active wpid
#   wormpep.history     dotname, wpid, start (version Nr.), stop (version Nr.) if inactive
########################################################################

 ############################################################################
 # write the new wp.fasta file, and wrap the sequences using rd's seqpress  #
 ############################################################################


print LOG &runtime, ": writing wormpep.fasta file\n\n";

open (WPFASTA , ">$new_wpdir/wp.fasta_unwrap$release");
for (my $i = 1 ; $i <= $wpmax ; $i++) {
    if (defined ($number2peptide[$i])) {
       my $i_pad = sprintf "%05d" , $i;
       print WPFASTA ">CE$i_pad\n$number2peptide[$i]\n";
    }
}
close WPFASTA;

system ("rewrap $new_wpdir/wp.fasta_unwrap$release > $new_wpdir/wp.fasta$release");

#system ("/nfs/disk100/wormpub/bin.ALPHA/seqpress -a $new_wpdir/wp.fasta_unwrap$release > $new_wpdir/wp.fasta$release");
#chmod (0444 , "$new_wpdir/wp.fasta$release") || print LOG "cannot chmod $new_wpdir/wp.fasta$release\n";

print LOG &runtime, ": finished writing wormpep.fasta file\n\n";

 ##################################################
 # retrieve most of the data from the database    #
 # via a tablemaker call                          #
 ##################################################

print LOG &runtime, ": retrieving data from autoace for each CDS\n\n";


my %CDS_locus      = ""; 
my %CDS_id         = ""; 
my %CDS_source     = ""; 
my %CDS_lab        = "";
my %laboratory     = "";
my %CDS_status     = "";
my @CDS            = "";
my %EST            = "";
our %Protein_id     = "";
our %Protein_ac     = "";
our %Protein_db     = "";

$ENV{'ACEDB'} = $dbdir;
my $command=<<EOF;
Table-maker -p "$dbdir/wquery/wormpep.def"
quit
EOF
   
open (TACE, "echo '$command' | $tace | ");
while (<TACE>) {
    print;
    chomp;
    s/acedb\> //g;      # only need this is using 4_9i code, bug fixed in 4_9k onward (should be redundant)
    next if ($_ eq "");
    next if (/\/\//);
    s/\"//g;
    (/^(\S+)\s/);
    my ($cds,$parent,$prot_id_parent,$prot_id,$prot_id_ver,$prot_db,$prot_ac,$loci,$confirmed,$est,$brief_id) = split /\t/;

    $CDS_locus{$cds}  = $loci;
    $CDS_id{$cds}     = $brief_id;
    $CDS_lab{$cds}    = $parent;
    $EST{$cds}        = $est;
    $Protein_ac{$cds} = $prot_ac;
    $Protein_id{$cds} = $prot_id . "." . $prot_id_ver;
    $Protein_id{$cds} = "" if $Protein_id{$cds} eq ".";
    $Protein_db{$cds} = $prot_db;

    print "$cds Locus:$CDS_locus{$cds} Protein_ID:$Protein_id{$cds} Protein_AC:$Protein_ac{$cds} Protein_DB:$Protein_db{$cds} Brief_ID:$CDS_id{$cds}\n\n";
    
    # confirmed CDS - [EST|mRNA] data
    if ($confirmed) {
	$CDS_status{$cds} = "Confirmed";
    }
    else {
	# supported CDS - [EST] data
	if ($est > 0) {
	    $CDS_status{$cds} = "Partially_confirmed";
	}
	# predicted CDS - no data
	else {
	    $CDS_status{$cds} = "Predicted";
	}
    }
}
close TACE;

my $command1=<<EOF;
Table-maker -p "$dbdir/wquery/wormpep2.def"
quit
EOF
    
open (TACE, "echo '$command1' | $tace | ");
while (<TACE>) {
    chomp;
    next if ($_ eq "");
    next if (/\/\//);
    s/acedb\>//g;
    s/\"//g;
    my ($clone,$lab) = (/^(\S+)\s(\S+)$/);
    $laboratory{$clone} = $lab;
}
close TACE;

print LOG &runtime, ": finished retrieving data from autoace for each CDS\n\n";

 ###########################################################################
 # get from autoace the data required to write the wormpep file and the wormpep.table file
 # use rd's seqpress to wrap the sequence lines in the wormpep file 
 ###########################################################################

print LOG &runtime, ": Build wormpep & wormpep.table files\n\n";

open (USERFASTA , ">$new_wpdir/wormpep_unwrap$release") || die "cannot create wormpep_unwrap$release\n";
open (USERTABLE , ">$new_wpdir/wormpep.table$release")  || die "cannot create wormpep.table$release\n";
USERFASTA->autoflush();
USERTABLE->autoflush();

foreach (@dotnames) {
    undef (my $locus_1);
    undef (my $locus);
    undef (my $brief_id);
    undef (my $cosmid);

    my $cds = $_;
    my $c_dot = $cds;
#    $c_dot =~ tr/a-z/A-Z/;
    my $wpid = $cds2num{$c_dot};
    my $wpid_pad = sprintf "%05d" , $wpid;
    my $pepseq = $number2peptide[$wpid];

    if ($CDS_locus{$cds} ne "") {
	$locus_1 = "locus\:".$CDS_locus{$cds};
    }
    
    if ( (defined ($Protein_id{$cds})) && ($Protein_id{$cds} ne "") ) {
	if ($Protein_db{$cds} eq "SwissProt") {
	    print USERFASTA ">$c_dot CE$wpid_pad  $locus_1 $CDS_id{$cds} status\:$CDS_status{$cds} SW\:$Protein_ac{$cds} protein_id\:$Protein_id{$cds}\n$pepseq\n";
	    print USERTABLE "$c_dot\tCE$wpid_pad\t$CDS_locus{$cds}\t$CDS_id{$cds}\t$CDS_status{$cds}\tSW\:$Protein_ac{$cds}\t$Protein_id{$cds}\n";
	}
	elsif ($Protein_db{$cds} eq "TREMBL") {
	    print USERFASTA ">$c_dot CE$wpid_pad  $locus_1 $CDS_id{$cds} status\:$CDS_status{$cds} TR\:$Protein_ac{$cds} protein_id\:$Protein_id{$cds}\n$pepseq\n";
	    print USERTABLE "$c_dot\tCE$wpid_pad\t$CDS_locus{$cds}\t$CDS_id{$cds}\t$CDS_status{$cds}\tTR\:$Protein_ac{$cds}\t$Protein_id{$cds}\n";
	}
	elsif ($Protein_db{$cds} eq "TREMBLNEW") {
	    print USERFASTA ">$c_dot CE$wpid_pad  $locus_1 $CDS_id{$cds} status\:$CDS_status{$cds} TN\:$Protein_ac{$cds} protein_id\:$Protein_id{$cds}\n$pepseq\n";
	    print USERTABLE "$c_dot\tCE$wpid_pad\t$CDS_locus{$cds}\t$CDS_id{$cds}\t$CDS_status{$cds}\tTN\:$Protein_ac{$cds}\t$Protein_id{$cds}\n";
	}
	else {
	    print USERFASTA ">$c_dot CE$wpid_pad  $locus_1 $CDS_id{$cds} status\:$CDS_status{$cds} protein_id\:$Protein_id{$cds}\n$pepseq\n";
	    print USERTABLE "$c_dot\tCE$wpid_pad\t$CDS_locus{$cds}\t$CDS_id{$cds}\t$CDS_status{$cds}\t$Protein_id{$cds}\n";
	}
    }
    else {
	if ($Protein_db{$cds} eq "SwissProt") {
	    print USERFASTA     ">$c_dot CE$wpid_pad  $locus_1 $CDS_id{$cds} status\:$CDS_status{$cds} SW\:$Protein_ac{$cds}\n$pepseq\n";
	    print USERTABLE      "$c_dot\tCE$wpid_pad\t$CDS_locus{$cds}\t$CDS_id{$cds}\t$CDS_status{$cds}\tSW\:$Protein_ac{$cds}\n";
	} 
	elsif ($Protein_db{$cds} eq "TREMBL") {
	    print USERFASTA ">$c_dot CE$wpid_pad  $locus_1 $CDS_id{$cds} status\:$CDS_status{$cds} TR\:$Protein_ac{$cds}\n$pepseq\n";
	    print USERTABLE  "$c_dot\tCE$wpid_pad\t$CDS_locus{$cds}\t$CDS_id{$cds}\t$CDS_status{$cds}\tTR\:$Protein_ac{$cds}\n";
	}
	elsif ($Protein_db{$cds} eq "TREMBLNEW") {
	    print USERFASTA ">$c_dot CE$wpid_pad  $locus_1 $CDS_id{$cds} status\:$CDS_status{$cds} TN\:$Protein_ac{$cds}\n$pepseq\n";
	    print USERTABLE  "$c_dot\tCE$wpid_pad\t$CDS_locus{$cds}\t$CDS_id{$cds}\t$CDS_status{$cds}\tTR\:$Protein_ac{$cds}\n";
	}
	else {
	    print USERFASTA ">$c_dot CE$wpid_pad  $locus_1 $CDS_id{$cds} status\:$CDS_status{$cds}\n$pepseq\n";
	    print USERTABLE  "$c_dot\tCE$wpid_pad\t$CDS_locus{$cds}\t$CDS_id{$cds}\t$CDS_status{$cds}\n";
	}
    }
} # end of loop foreach (@dotnames)

close USERFASTA;
close USERTABLE;

system ("rewrap $new_wpdir/wormpep_unwrap$release > $new_wpdir/wormpep$release");


#system ("/nfs/disk100/wormpub/bin.ALPHA/seqpress -a $new_wpdir/wormpep_unwrap$release > $new_wpdir/wormpep$release");
#chmod (0444 , "$new_wpdir/wormpep$release") || print LOG "cannot chmod $new_wpdir/wormpep$release\n";

if( $final ){
  #chmod (0444 , "$new_wpdir/wormpep.table$release") || print LOG "cannot chmod $new_wpdir/wormpep.table$release\n";
  ###########################################################################
  # create a blast'able database (indexing) using setdb for Wublast (not formatdb, which is  for blastall)
  
  system ("cp $new_wpdir/wormpep$release $new_wpdir/wormpep_current");
  system ("/usr/local/pubseq/bin/setdb $new_wpdir/wormpep_current > $new_wpdir/wormpep_current.log");
  
  ###########################################################################
}

else {
  `rm -f $new_wpdir/wormpep.table$release`;
}

print LOG &runtime, ": Finished building wormpep & wormpep.table files\n\n";

# count the CDS (with and without alternate splice forms) based on the wormpep file

my @end_dotnames;
open (FILE , "$new_wpdir/wormpep$release") || print LOG  "cannot open the wormpep$release file\n";
while (<FILE>) {
    if (/^\>(\S+)\s+\S+/) {
       my $cds = $1;
       push (@end_dotnames , $cds);
    }
}
close FILE;
my $number_cds = 0;
my $number_total = 0;
my $number_alternate = 0;
my %new_name2x;
foreach (@end_dotnames) {
    $number_total++;
    /^([A-Z0-9_]+)\.(.*)$/i;
    my $cosmid = $1;
    my $cds = $2;
    if ($cds =~ /^[1-9][0-9]?$/) {
        $number_cds++;
    }
    elsif ($cds =~ /(^[1-9][0-9]?)([a-z])/) {
	my $number = $1;
	my $letter = $2;
        my $new_name = $cosmid."_".$number;
        unless (exists ($new_name2x{$new_name})) {
	    $number_cds++;
        } 
	else {
	    $number_alternate++;
        }
        $new_name2x{$new_name} = "x";
    } 
    else {
        print LOG "$_ has a non\-acceptable name in wormpep$release \(has not been counted\)\n";
        next;
    }
}
print LOG "\n\nthere are $number_cds CDS in autoace, $number_total when counting \($number_alternate\) alternate splice_forms\n";

###########################################################################
# create the new wormpep.accession,
# wp.table contains tab-separated:  wpid, list of associated dotnames if active, 
# empty if inactive, and link to active wpid in case of duplicate sequences in wp.fasta
if( $final ) {
  open (WPTABLE , ">$new_wpdir/wormpep.accession$release") || die "cannot create wormpep.accession$release\n";
  
  for (my $i = 1 ; $i <= $wpmax ; $i++) {
    if (defined ($number2accession[$i])) {
      print WPTABLE "$number2accession[$i]\n";
    }
  }
  
  close WPTABLE;
  
  #chmod (0444 , "$new_wpdir/wormpep.accession$release") || print LOG "cannot chmod $new_wpdir/wormpep.accession$release\n";
}
############################################################################
# read in the current wormpep.history file, update it, and read it back out,
# wormpep.history contains tab-separated:  dotname, wpid, start (release), end (release)

open (OLDHISTORY, "$wpdir/wormpep.history$old_release")  || die "cannot open $wpdir/wormpep.history$old_release\n";
open (HISTORY,    ">$new_wpdir/wormpep.history$release") || die "cannot create $wpdir/wormpep.history$release\n";
open (DIFF,       ">$new_wpdir/wormpep.diff$release")    || die "cannot create $wpdir/wormpep.diff$release\n";

my %line;
while (my $line = <OLDHISTORY>) {
    chomp $line;
    my $cdsname = ""; my $wpid = ""; my $start = ""; my $end = "";
    ($cdsname , $wpid , $start , $end) = split (/\t/ , $line);
    my $c_dot = $cdsname;
    $wpid =~ /CE0*([1-9]\d*)/ ; my $num = $1;
    $line{$c_dot} = $line;
    if ((!exists ($cds2num{$c_dot}) && ($end eq ""))) {
       print HISTORY "$c_dot\t$wpid\t$start\t$release\n";
       print DIFF "lost:\t$c_dot\t$wpid\n";
   } elsif (($cds2num{$c_dot} ne $num) && ($end eq "")) {
        print HISTORY "$c_dot\t$wpid\t$start\t$release\n";
        my $new_num = $cds2num{$c_dot};
        my $new_pad = sprintf "%05d" , $new_num;
        print HISTORY "$c_dot\tCE$new_pad\t$release\t$end\n";
        print DIFF "changed:\t$c_dot\t$wpid --> CE$new_pad\n";
     } else {
         print HISTORY "$c_dot\t$wpid\t$start\t$end\n";
       }
}

foreach (keys (%cds2num)) {
    my $empty = "";
    my $c_dot = $_;
    if (!exists ($line{$c_dot})) {
        my $num = $cds2num{$c_dot};
        my $pad = sprintf "%05d" , $num;
        print HISTORY "$c_dot\tCE$pad\t$release\t$empty\n";
        print DIFF "new:\t$c_dot\tCE$pad\n";
        next;
    }
    my $cdsname = ""; my $wpid = ""; my $start = ""; my $end = "";
    ($cdsname , $wpid , $start , $end) = split (/\t/ , $line{$c_dot});
    if ($end ne "") {   
        my $new_num = $cds2num{$c_dot};
        my $new_pad = sprintf "%05d" , $new_num;
        print HISTORY "$c_dot\tCE$new_pad\t$release\t$empty\n";
        print DIFF "reappeared:\t$c_dot\tCE$new_pad\n";
    }
}

close OLDHISTORY;
close HISTORY;
close DIFF;
#chmod (0444 , "$new_wpdir/wormpep.history$release") || print LOG "cannot chmod $new_wpdir/wormpep.history$release\n";
#chmod (0444 , "$new_wpdir/wormpep.diff$release")    || print LOG "cannot chmod $new_wpdir/wormpep.diff$release\n";
my $wpdiff = $wpmax - $old_wpmax;
print LOG "\n\nnew wpmax of wp.fasta$release equals $wpmax\n$wpdiff new sequences have been added\n";


##############################
# mail $maintainer report    #
##############################

unlink ("$new_wpdir/wp.fasta_unwrap$release") || print LOG "cannot delete $new_wpdir/wp.fasta_unwrap$release\n";
unlink ("$new_wpdir/wormpep_unwrap$release")  || print LOG "cannot delete $new_wpdir/wormpep_unwrap$release\n";


close LOG;

#Email log
&mail_maintainer("WormBase Build Report: make_wormpep",$maintainers,$log);


# write the release letter (also does some checks)
if( $final ) {
  &release_wormpep($number_cds,$number_total,$number_alternate);
  chmod (0444 , "$new_wpdir/*") || print LOG "cannot chmod $new_wpdir/ files\n";
}

 ##############################
 # hasta luego                #
 ##############################

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
  
  if ($verbose) {
    print "Usage : make_wormpep [-options]\n";
    print "=============================================\n";
    print " -r <int>     : release version number\n";
    print " -h           : help pages   \n";
    print " -d           : verbose (debug) mode\n";
    print "=============================================\n";
    print "\n";
  }
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
# process_cds_class
# 
# 1) retrieves list of valid ?CDS objects
# 2) gets dna and peptide sequence for each ?CDS, writes dna.fasta file,
# 3) maps the peptide sequences onto wormpep, deleting Peptides containing an X or a *,
# 4) creates new Wormpep IDs if necessary, creating %dot2num
#
##########################################################################################

sub process_cds_class{

  print LOG &runtime, ": connecting to $dbdir\n";
  print     &runtime, ": connecting to $dbdir\n" if ($verbose);

  # grab list of valid CDS names, ignore any temp genes
  my $db = Ace->connect (-path => $dbdir, -program => $tace) || &error(6); # die "cannot connect to autoace\n";
  my @CDSs = $db->fetch (-query => 'FIND elegans_CDS NOT *temp*');
  @CDSs = sort @CDSs;

  print LOG "=> " . scalar(@CDSs) . " CDSs\n";
  print     "=> " . scalar(@CDSs) . " CDSs\n" if ($verbose);


  # write DNA for each CDS to separate file
  print LOG &runtime, ": creating wormpep.dna file\n\n";
  open (DNA , ">$new_wpdir/wormpep.dna$release") || &error(7); # die "cannot create $new_wpdir/wormpep.dna$release\n";


  foreach my $cds (@CDSs) {
    
    print "Processing: $cds\n" if ($verbose);
    
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
    
    my $c_dot = $cds;

    # check current peptide sequence against peptides loaded in from last releases wp.fasta file,
    # if peptide is new add to hash
    unless (exists ($peptide2number{$peptide_seq})){
      $number2peptide[++$wpmax] = $peptide_seq; # adds new sequence to array, increases $wpmax
      $peptide2number{$peptide_seq} = $wpmax;
      $cds2num{$c_dot} = $wpmax;
      my $pad = sprintf "%05d" , $wpmax;
      $number2accession[$wpmax] = "CE$pad\t$c_dot";
    } 
    else {
      $cds2num{$c_dot} = $peptide2number{$peptide_seq};
      my $number = $peptide2number{$peptide_seq};
      $number2accession[$number] .= "\t$c_dot";
    }             
    push (@dotnames, $cds);
    push (@c_dotnames, $c_dot);
    
    # tidy up
    $cds->DESTROY();
    $peptide->DESTROY();
    
    
  }   
  close DNA;
  print LOG &runtime, ": finished writing wormpep.dna file\n\n";

  # close database connection
  $db->close;
  print LOG &runtime, ": finished connection to database\n\n";
  print     &runtime, ": finished connection to database\n\n" if ($verbose);

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
  $log        = "$basedir/logs/$script_name.initial.$rundate.$$" if ($initial);
  $log        = "$basedir/logs/$script_name.final.$rundate.$$"   if ($final);

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




















