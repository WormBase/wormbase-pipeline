#!/software/bin/perl -w
use strict;
use lib $ENV{'CVS_DIR'};
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
#use lib '/nfs/users/nfs_g/gw3/Nameserver-API';

use NameDB_handler;
use Getopt::Long;
use Log_files;
use Ace;
use Wormbase;

=pod

=head batch_split.pl

=item Options:

  -test      use the test nameserver
  -file	     file containing genes to split <Mandatory>

    FORMAT:

splitgene.pl -old WBGene00238231 -new OVOC13433 -who 2062 -id WBGene00255445 -load -species ovolvulus
splitgene.pl -old WBGene00238817 -new OVOC13435 -who 2062 -id WBGene00255447 -load -species ovolvulus
splitgene.pl -old WBGene00239023 -new OVOC13436 -who 2062 -id WBGene00255448 -load -species ovolvulus
splitgene.pl -old WBGene00239126 -new OVOC13416 -who 2062 -id WBGene00255425 -load -species ovolvulus
splitgene.pl -old WBGene00239216 -new OVOC13438 -who 2062 -id WBGene00255450 -load -species ovolvulus
splitgene.pl -old WBGene00239809 -new OVOC13489 -who 2062 -id WBGene00255515 -load -species ovolvulus
splitgene.pl -old WBGene00239863 -new OVOC13439 -who 2062 -id WBGene00255451 -load -species ovolvulus
splitgene.pl -old WBGene00240141 -new OVOC13437 -who 2062 -id WBGene00255449 -load -species ovolvulus


  -debug     limits to specified user <Optional>
  -load      loads the resulting .ace file into geneace.
  -ns        update the nameserver (the content of the secnd WBGeneID column will be ignored)

e.g. perl batch_split.pl -file genesplits.txt


=cut

my ($test, $file, $load, $ns);
my ($help, $debug, $verbose, $store, $wormbase, $species);

GetOptions(
	   "test"       => \$test,
	   "store:s"    => \$store,
	   "species:s"  => \$species,
	   'file:s'     => \$file,
	   'debug:s'    => \$debug,
	   'load'       => \$load,
           'ns'         => \$ns,
	  ) or die;

my $log;

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -organism => $species,
			     -test => $test,
                             );
}


if (defined $debug) {$log = Log_files->make_log("NAMEDB:$file", $debug);}
else {$log = Log_files->make_log("NAMEDB:$file");}


my $db;
if ($ns) {
  $log->write_to("Contacting NameServer .....\n");
  $db = NameDB_handler->new($wormbase, $test);
}

my $database = "/nfs/wormpub/DATABASES/geneace";
$log->write_to("Working.........\n-----------------------------------\n\n\n1) splitting genes in file [${file}]\n\n");

my $ace = Ace->connect('-path', $database) or $log->log_and_die("cant open $database: $!\n");

my $outdir = $database.'/NAMEDB_Files/';
my $backupsdir = $outdir.'BACKUPS/';
my $outname = 'batch_split.ace';
my $outputfile = "$outdir"."$outname";
my $output;
my %gene_versions; # remember the latest version used in all genes altered in case a gene is being split more than once

##############################
# warn/notify on use of -load.
##############################
if (!defined $load) {$log->write_to("You have decided not to automatically load the output of this script\n\n");}
elsif (defined $load) { $log->write_to("Output has been scheduled for auto-loading.\n\n");}

# open file and read
open (FILE,"<$file") or $log->log_and_die("can't open $file : $!\n");
open (ACE,">$outputfile") or $log->log_and_die("cant write output: $!\n");

my($oldgene,$newgene,$newname, $user);

my $malformedcount=0;
my $actualcount=0;
my $count=0;

while (<FILE>) {
  chomp;
  
  # splitgene.pl -old WBGene00243151 -new OVOC13459 -who 2062 -id WBGene00255474 -load -species ovolvulus
  if (/splitgene.pl\s+-old\s+(WBGene\d{8})\s+-new\s+(\S+)\s+-who\s+(\d+)\s+-id\s+(WBGene\d{8})\s+-load\s+-species\s+(\w+)/) { #gather info
    # Captured string ($1) - WBGene00243151
    # Captured string ($2) - OVOC13459
    # Captured string ($3) - 2062
    # Captured string ($4) - WBGene00255474 - this is a dummy field, it is not used if '-ns' is given on the command line
    # Captured string ($5) - ovolvulus
    
    $oldgene = $1;
    $newname = $2;
    $user = $3;
    $newgene = $4;
    $species = $5;
    $actualcount++;
    
    if ($user !~ /^WBPerson/) {
      $user = "WBPerson" . $user;
    }

    &split_gene;
  }
  elsif (/\w+/) {
    $log->error("ERROR: $_ is a malformed line which appears to not include all of the information required\n");
    $malformedcount++;
  }
  else {
    next;
  }
}

close(ACE);

$log->write_to("$actualcount gene pairs in file to be split\n\n");
$log->write_to("$count gene pairs split ($malformedcount malformed_lines)\n\n");
&load_data if $load;
$log->write_to("Check $outputfile file and load into geneace.\n") unless $load;
$log->mail();

exit(0);

###############################################################################################

sub split_gene {
  my $ok;
  if($oldgene and $newgene and $user and $newname) {
    my $species_full_name;
    $output = "";
    $ok = 1; # error status

    if ($ns) {
      # we assume that the new gene biotype is 'CDS'
      $newgene = $db->split_gene($oldgene, $newname, 'CDS', $species);      
    }

    # does the split_into gene already exist?
    my $newgeneObj = $ace->fetch('Gene', $newgene);
    if ($newgeneObj) {
      $log->error("ERROR: $newgene already exists\n");
      $ok = 0;
    }

    # process LIVE gene
    my $oldgeneObj = $ace->fetch('Gene', $oldgene);
    if ($oldgeneObj) {
      # is this a Live gene with no splits to the new gene in the last Split_into?
      my $status = $oldgeneObj->Status->name;
      $species_full_name = $oldgeneObj->Species->name;
      if ($status ne 'Live') {
	$log->error("ERROR: $oldgene is not a Live gene\n");
	$ok = 0;
      }
      # get the last acquires_merge
      my $split_into;
      foreach my $split_into_obj ($oldgeneObj->at('Identity.History.Split_into')) {
	if (defined $split_into_obj) {
	  ($split_into) = $split_into_obj->row;
	}
      }
      if (defined $split_into && $split_into eq $newgene) {
	$log->error("Warning: $oldgene has a tag saying it has already had $newgene split from it - this split will not be done again\n");
	$ok = 0;
      }

      # get the version
      my $ver;
      if (exists $gene_versions{$oldgene}) {
	$ver = $gene_versions{$oldgene};
      } else {
	$ver = $oldgeneObj->Version->name;
      }
      $ver++;
      $gene_versions{$oldgene} = $ver;

      $output .= "\nGene : $oldgene\nVersion $ver\nHistory Version_change $ver now $user Event Split_into $newgene\nSplit_into $newgene\n";
      
    } else {
      $log->error("ERROR: no such gene $oldgene\n");
      $ok = 0;
    }


    # process NEW gene
    my $ver = "1";
    $output .= "\nGene : $newgene\nVersion $ver\nSequence_name $newname\nPublic_name $newname\nSpecies \"$species_full_name\"\nHistory Version_change $ver now $user Event Split_from $oldgene\nSplit_from $oldgene\nLive\nMethod Gene\n\n";
    
    
    
  } else {
    $log->error("ERROR: Missing information to create $newgene\n");
    $ok = 0;
  }
  
  # we did this one successfully
  if ($ok) {
    print ACE $output;
    $count++;
  }
  else {
    $log->error("ERROR: Too many issues with the $oldgene :: $newgene split, not processing\n");
  }
  
  undef $oldgene; undef $newgene ;undef $user; undef $newname;
}

sub load_data {
  # load information to $database if -load is specified
  $wormbase->load_to_database("$database", "$outputfile", 'batch_split.pl', $log, undef, 1);
  $log->write_to("Loaded $outputfile into $database\n\n");
  $wormbase->run_command("mv $outputfile $backupsdir"."$outname". $wormbase->rundate. "\n"); #append date to filename when moving.
  $log->write_to("Output file has been cleaned away like a good little fellow\n\n");
  print "Finished!\n";
}
