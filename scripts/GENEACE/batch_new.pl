#!/software/bin/perl -w
use strict;
use lib '../blib/lib';
use lib '/nfs/WWWdev/SANGER_docs/lib/Projects/C_elegans';
use lib $ENV{'CVS_DIR'};
use NameDB_handler;
use Getopt::Long;
use Log_files;
use Ace;
use Wormbase;
=pod

=head batch_new.pl

=item Options:

  -file	     file containing genes to create <Mandatory>

    FORMAT:

newgene.pl -seq Bm16920 -who 4055 -load -id WBGene00255463 -species brugia
newgene.pl -seq Bm16921 -who 4055 -load -id WBGene00255483 -species brugia
newgene.pl -seq Bm16922 -who 4055 -load -id WBGene00255484 -species brugia
newgene.pl -seq Bm16923 -who 4055 -load -id WBGene00255501 -species brugia

  -debug     limits to specified user <Optional>
  -load      loads the resulting .ace file into geneace.

e.g. perl batch_new.pl -file newgenes.txt


=cut

my ($USER, $test, $file, $debug, $load,$species);
GetOptions(
	   'user:s'     => \$USER,
	   'test'       => \$test,
	   'file:s'     => \$file,
	   'debug:s'    => \$debug,
	   'load'       => \$load,
           'species:s'  => \$species,
	  ) or die;


my $log;
if (defined $USER) {$log = Log_files->make_log("NAMEDB:$file", $USER);}
elsif (defined $debug) {$log = Log_files->make_log("NAMEDB:$file", $debug);}
else {$log = Log_files->make_log("NAMEDB:$file");}
my $DB;
my $db;
my $wormbase = Wormbase->new("-organism" =>$species);
my $database = "/nfs/wormpub/DATABASES/geneace";
$log->write_to("Working.........\n-----------------------------------\n\n\n1) creating genes in file [${file}]\n\n");
$log->write_to("TEST mode is ON!\n\n") if $test;

my $ace = Ace->connect('-path', $database) or $log->log_and_die("cant open $database: $!\n");


my $outdir = $database."/NAMEDB_Files/";
my $backupsdir = $outdir."BACKUPS/";
my $outname = "batch_new.ace";
my $outputfile = "$outdir"."$outname";
my $output;



#generate a list of corespecies short::ful names
my %full_name_data;

my %accessors = ($wormbase->species_accessors);
  foreach my $wb (values %accessors) {
    my $gspecies = $wb->full_name;
    $full_name_data{$wb->species} = $wb->full_name;
  }

#store the IDs being created in the file to see if they have already been used in this batch

my @genes;
my @genenames;


##############################
# warn/notify on use of -load.
##############################
if (!defined$load) {$log->write_to("2) You have decided not to automatically load the output of this script\n\n");}
elsif (defined$load) { $log->write_to("2) Output has been scheduled for auto-loading.\n\n");}

#open file and read
open (FILE,"<$file") or $log->log_and_die("can't open $file : $!\n");
open (ACE,">$outputfile") or $log->log_and_die("cant write output: $!\n");
my($oldgene,$newgene,$newname,$user);
my $count=0;
my $createdcount=0;
while (<FILE>) {
  chomp;
  
  #newgene.pl -seq Bm16920 -who 4055 -load -id WBGene00255463 -species brugia
  if (/newgene.pl\s+-seq\s+(\w+)\s+-who\s+(\d+)\s+\S+\s+-id\s+(WBGene\d{8})\s+-species\s+(\w+)/) { #gather info
    #Captured string ($1) - Bm16922
    #Captured string ($2) - 4055
    #Captured string ($3) - WBGene00255484
    #Captured string ($4) - brugia
    $createdcount++;
    $newname = $1;
    $user = "WBPerson$2";
    $newgene = $3;
    $species = $4;
    &create_gene;
  }
  elsif (/\w+/) {
    $log->error("ERROR: $_ is a malformed line which appears to not include all of the information required\n");
  }
  else {
    next;
  }
}

close(ACE);
$log->write_to("3) $createdcount genes in file to be created\n\n");
$log->write_to("4) $count genes created\n\n");
&load_data if ($load);
$log->write_to("5) Check $outputfile file and load into geneace.\n") unless ($load);
$log->mail();
exit(0);

###############################################################################################

sub create_gene {
  my $ok;
  if($newgene and $user and $newname) {
    $output = "";
    $ok = 1; # error status

    #Does the new gene already exist?
    my $newgeneObj = $ace->fetch('Gene', $newgene);
    if ($newgeneObj) {
      $log->error("ERROR: $newgene already exists\n");
      $ok = 0;
      }
    my $newseqObj = $ace->fetch('Gene_name', $newname);
    if ($newseqObj) {
      my $error_gene = $newseqObj->Sequence_name_for->name; #check this
      $log->error("ERROR: $newname already exists as sequence name for $error_gene\n");
      $ok = 0;
    }

    if ( grep( /$newgene$/, @genes) ) {
      $log->error("ERROR: $newgene has already been seen in this batch\n");
      $ok = 0;
    }
    if ( grep( /$newname$/, @genenames) ) {
      $log->error("ERROR: $newname has already been seen in this batch\n");
      $ok = 0;
    }
    
    # process NEW gene
    my $ver = "1";
    $output .= "\nGene : $newgene\nVersion $ver\nSequence_name $newname\nPublic_name $newname\nSpecies \"$full_name_data{$species}\"\nHistory Version_change $ver now $user Event Created\nLive\nMethod Gene\n\n";
    
  } else {
    $log->error("ERROR: Missing information to create $newgene\n");
    $ok = 0;
  }
  
  # we did this one successfully
  if ($ok) {
    push(@genes,"$newgene");
    push(@genenames, "$newname");
    print ACE $output;
    $count++;
  }
  else {
    $log->error("ERROR: Too many isses with the ($newgene/$newname), not processing\n\n");
  }
  
  undef $newgene ;undef $user; undef $newname;
}

sub load_data {
# load information to $database if -load is specified
$wormbase->load_to_database("$database", "$outputfile", 'batch_new.pl', $log, undef, 1);
$log->write_to("5) Loaded $outputfile into $database\n\n");
$wormbase->run_command("mv $outputfile $backupsdir"."$outname". $wormbase->rundate. "\n"); #append date to filename when moving.
$log->write_to("6) Output file has been cleaned away like a good little fellow\n\n");
print "Finished!!!!\n";
}
