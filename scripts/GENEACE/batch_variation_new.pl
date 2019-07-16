#!/software/bin/perl -w
use strict;
use lib $ENV{'CVS_DIR'};
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";

use NameDB_handler;
use Getopt::Long;
use Log_files;
use Ace;
use Wormbase;
=pod

=head batch_variation_new.pl

=item Options:


  -file	     file containing variation information from the nameserver emails. <Mandatory>

  -override  This option forces the geneace status checks to be ignored for pushing 
             merges into the Nameserver when the changes have already been done in geneace
             via a batch adhoc submission. Cannot be used with the -load option

  -force     This option forces the .ace file to be loaded into geneace when there are Gene
             names involved as it's some what easier to resolve issues post merge instead of 
             overlooking when merges have failed because of this.

    FORMAT:

          EMAIL : daniela@wormbase.org
           NAME : Daniela Raciti
           Name : sy682
         Remark : WBPaper00053325, lov-3 whittaker micropublication
        WBVarID : WBVar02148709

e.g. perl batch_merge.pl -file merger.txt


=cut

my ($USER, $test, $file, $debug, $load, $old, $ns, $PASS,$out,$force,$override);
GetOptions(
	   'file:s'     => \$file,
	   'debug:s'    => \$debug,
	   'load'       => \$load,
	   'out:s'      => \$out,
	   'force'      => \$force,
	   'override'   => \$override,
	  ) or die;

if ($debug) {
    $USER = $debug;
}
else {
    $USER = "pad";
}

my $species = 'elegans';
my $log;
if (defined $USER) {$log = Log_files->make_log("NAMEDB:$file", $USER);}
elsif (defined $debug) {$log = Log_files->make_log("NAMEDB:$file", $debug);}
else {$log = Log_files->make_log("NAMEDB:$file");}
my $DB;
my $db;
my $ecount = 0;
my $countp = 0;

if ($load && $override) {
   $log->log_and_die("OPTIONS ERROR: You cannot use the -load and -override option at the same time!!!!!!!\n");
}

my $wormbase = Wormbase->new("-organism" =>$species);
my $database = "/nfs/wormpub/DATABASES/geneace";
$log->write_to("Working.........\n-----------------------------------\n\n\n1) creating variations in file [${file}]\n\n");


my $ace = Ace->connect('-path', $database) or $log->log_and_die("cant open $database: $!\n");

my $outdir = $database."/NAMEDB_Files/";
my $backupsdir = $outdir."BACKUPS/";
$out ||= "batch_variations_${file}.ace";
my $output = $outdir.$out;

my %gene_versions; # remember the latest version used in all genes altered in case a gene is being merged in to more than once

##############################
# warn/notify on use of -load.
##############################
if (!defined$load) {$log->write_to("2) You have decided not to automatically load the output of this script\n\n");}
elsif (defined$load) { $log->write_to("2) Output has been scheduled for auto-loading.\n\n");}

#open file and read
open (FILE,"<$file") or $log->log_and_die("can't open $file : $!\n");
open (ACE,">$output") or $log->log_and_die("cant write output: $!\n");
my($varid,$varname,$user,$remark,$paper);
my $count=0;
while (<FILE>) {
    chomp;
    print "$_\n" if ($debug);
    unless (/\w/) {
	&new_var;
    }
    else { #gather info
	if   (/\s+WBVarID\s+:\s+(WBVar[0-9]{8})/) {
	    $varid = $1;
	} 
	elsif(/\s+Name\s+:\s+(\S+)/) { 
	    $varname = $1; 
	} 
	elsif(/\s+NAME\s+:\s+(\S+\s+\S+)/) {
	    $user = $1;
	}
	elsif(/\s+Remark\s+:\s+(.*)$/){
	    $remark = $1;
	} # ignore this line
	elsif(/\s+NAME/){} # ignore this line
	elsif(/\s+EMAIL/){}
	else { 
	    $log->error("malformed line : $_\n")
	}
    }
}

&new_var; # remember the last one!
close(ACE);
$log->write_to("3) $count variations in file to be created\n\n");
$log->write_to("4) $countp variations created\n\n");
&load_data if ($load);
$log->write_to("5) Check $output file and load into geneace.\n") unless ($load);
$log->mail();
exit(0);

###############################################################################################

sub new_var {
  if($varid and $varname and $user) {
    print "\nCreating $varname into $varid\n\n";
    my $output = "";
    my $ok = 1; # error status

    # process LIVE gene
    my $varidObj = $ace->fetch('Variation', $varname);
    if (defined $varidObj) {
	print "\nVariation already in the database\n";
	$ok = 1;
    }
    else {
	$countp++;
	print "\nVariation ${varid}(${varname}) is new......processing\n";
	if ($remark =~ /(WBPaper[0-9]{8})/) {
	    $paper = $1;
	}
	$output .= "\nVariation : $varid\nPublic_name $varname\nSeqStatus Pending_curation\nAllele\nSpecies \"Caenorhabditis $species\"\nStatus Live\nMethod Allele\n";
	if ($paper) {
	    $output .= "Evidence Paper_evidence $paper\nReference $paper\n";
	}
	
	if ($user) {
	    $output .= "Remark \"Created in the nameserver by $user\"\n";
	}
	if ($remark) {
	    $output .= "Remark \"$remark\"\n";
	}
	$ok = 0;
    }
    undef $varid; undef $varname ;undef $user;undef $paper;
    # we did this one successfully
    if ($ok ne "1") {
      print ACE $output;
      $count++;
    }
  } elsif (!defined($varid && $varname && $user)) {
    $log->write_to("Warning: additional blank line in input file has been ignored\n");
  } else {
    $log->error("missing info on $varid : $varname : $user\n");
  }
  undef $varid; undef $varname ;undef $user;
}

sub load_data {
# load information to $database if -load is specified
$wormbase->load_to_database("$database", "$output", 'batch_merge.pl', $log, undef, 1);
$log->write_to("5) Loaded $output into $database\n\n");
$wormbase->run_command("mv $output $backupsdir".$out. $wormbase->rundate. "\n"); #append date to filename when moving.
$log->write_to("6) Output file has been cleaned away like a good little fellow\n\n");
print "Finished!!!!\n";
}
