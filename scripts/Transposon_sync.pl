#!/software/bin/perl -w
#
# Transposon_sync.pl                           
# 
# by pad                         
#
# Script to check all Transposon gene connections are still live and all genes have been treated correctly.
#

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
#use Ace;
#use Sequence_extract;
#use Coords_converter;§

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $build, $species);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "build"      => \$build,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -species => $species,
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

# Check options and usage
$log->log_and_die("-species is mandetory...sorry\n") unless ($species);
if ($debug) {print "running in debug mode, only $debug will receive the log report\n";}
if ($verbose) {print "running with verbose option, consider pipint the output somewhere next time\n";}
# in test mode?
if ($test) {print "In test mode which is the same as normal mode\n";}
if ($build) {print "Running in build mode, all queries will be against autoace\n"}
# Display help if required
&usage("Help") if ($help);


#################################
# Set up some useful paths      #
#################################

# Set up top level base directories (these are different if in test mode)
my $basedir         = $wormbase->basedir;     # BASE DIR
my $ace_dir         = $wormbase->autoace;     # AUTOACE DATABASE DIR
my $wormpep_dir     = $wormbase->wormpep;     # CURRENT WORMPEP
my $wormrna_dir     = $wormbase->wormrna;     # CURRENT WORMRNA
my $common_data_dir = $wormbase->common_data; # AUTOACE COMMON_DATA
my $chromosomes_dir = $wormbase->chromosomes; # AUTOACE CHROMSOMES
my $reports_dir     = $wormbase->reports;     # AUTOACE REPORTS
my $gff_dir         = $wormbase->gff;         # AUTOACE GFF
my $gff_splits_dir  = $wormbase->gff_splits;  # AUTOACE GFF SPLIT
my $logs_dir        = $wormbase->logs;        # AUTOACE LOGS

# some database paths
my ($seqdb, $geneace);
if (defined $build) {
  $seqdb = $wormbase->autoace;
  $geneace = $wormbase->autoace;
  $log->log_and_die("autoace is not at the expected location, has the build finished?/n") unless (-e $wormbase->autoace."/database/block1.wrm");
}
else {
  if ($species eq 'elegans') {
    $seqdb = $wormbase->database('camace');
  }
  else {
    $seqdb = $wormbase->database($species);
  }
  $geneace   = $wormbase->database('geneace');
}

# other paths
my $tace            = $wormbase->tace;        # TACE PATH
my $giface          = $wormbase->giface;      # GIFACE PATH




##########################
# MAIN BODY OF SCRIPT
##########################

# main stuff goes here
my $db = Ace->connect(-path=>$seqdb) or  $log->log_and_die("Couldn't connect to $seqdb\n". Ace->error);
my $gdb =  Ace->connect(-path=>$geneace) or  $log->log_and_die("Couldn't connect to $geneace\n". Ace->error);

$log->write_to("Using : $db for annotation data\nUsing : $gdb for primary gene data.\n\n");

# Query the geneace WBGene objects that have been flagged Transposon_in_origin
my @TPWBGenes;
my $ccount = 0;
my $gcount = 0;
my $ecount1 = 0;
my $ecount2 = 0;

my @SeqGenes = $gdb->fetch (-query => "FIND Gene WHERE Transposon_in_origin");
my @Generef;

foreach my $SeqGene(@SeqGenes) {
  $ccount ++;
  print "Checking ".$SeqGene->name."\n" if ($verbose);
  push (@Generef,"$SeqGene");
  if ($SeqGene->Status->name eq 'Suppressed'){
    $log->write_to("$SeqGene - Correct Status\n") if ($verbose);
    next;
  }
  elsif ($SeqGene->Status->name eq 'Dead'){
    my $tmpTP_ID = $SeqGene->Sequence_name->name; 
    $log->write_to("Warning:$SeqGene was Transposon_in_origin but is Dead, check $tmpTP_ID is also dead.\n");
    if ($SeqGene->Species->name eq $wormbase->full_name($species)) {
      my @tmpTPcheck = $db->fetch (-query => "FIND CDS $tmpTP_ID");
      foreach my $tmpTPcheck (@tmpTPcheck) {
	print "test\n"; 
      }
    }
    else {
      $log->write_to("You need to check this manually\n");
    }
    next;
  } 
  else {
      $log->write_to ("$SeqGene - ERROR Not Suppressed but is flagged as a Transposon gene in geneace\n");
      $ecount1++;
  }
}


# Query the species database for WBTransposon associated genes and corresponding CDS object genes

# Build Tablemaker query
my $query = &seqTPs();
$log->write_to("\nRetrieving transposon gene data from $seqdb\n") if ($verbose);
my $command = "Table-maker -p $query\nquit\n";
my $tacefh;

open($tacefh,  "echo '$command' | $tace $seqdb |");
## Gene Gene data ##
print "Table_maker_finished\n" if ($verbose);

while(<$tacefh>) {
  
  unless (/\"/) {next;}
  chomp; 
  s/\"//g;
  s/\\//g;
  #print "$_\n";
  my $F;
  my @F = split"\t",$_;
  
  if (defined $F[0]) {
    if ( grep( /^$F[0]$/, @TPWBGenes ) ) {
      print "Ignore - $F[0] is already in the array\n" if ($verbose);
      next;
    }
    else {
      push (@TPWBGenes, "$F[0]");
    }
  }
  if (defined $F[1]) {
    if ( grep( /^$F[1]$/, @TPWBGenes ) ) {
      print "Ignore - $F[1] is already in the array\n" if ($verbose);
      next;
    }
    else {
      push (@TPWBGenes, "$F[1]"); 
    }
  }
}

#Test to see if the WBGene ID was seen in Geneace as a valid Transposon_origin gene.
my $testgene;
foreach $testgene (@TPWBGenes){
  $gcount++;
  if ( grep( /^$testgene$/, @Generef ) ) {
    $log->write_to("Known $testgene is a Transposon gene\n") if ($verbose);
    next;
  }
  else {
    $log->write_to("$testgene - ERROR Not been processed correctly. This will need a change of class in the nameserver or if the gene is dead in the nameserver then the action has not been done in geneace\n");
    $ecount2++;
    next;
  }
}

# Close log files and exit
$log->write_to("\n\nStatistics\n");
$log->write_to("----------\n");
$log->write_to("$ccount objects retrieved from $gdb and considered a Transposon gene\n");
$log->write_to("$ecount1 ID error(s)\n");
$log->write_to("$gcount objects retrieved from $db and considered a Transposon gene\n");
$log->write_to("$ecount2 ID error(s)\n");
$log->mail();
print "Finished.\n" if ($verbose);
exit(0);

sub seqTPs {
  my ($full_species) = @_;

  my $tmdef2 = "/tmp/gene_tmqueryvar.$$.def";
  open my $qfhv, ">$tmdef2" or 
      $log->log_and_die("Could not open $tmdef2 for writing\n");  

  my $condition = "";

  my $tablemaker_template2 = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 20 
Optional 
Hidden 
Class 
Class Transposon 
From 1 
Condition WBTransposon*
 
Colonne 2 
Width 20 
Optional 
Visible 
Class 
Class Gene 
From 1 
Tag Gene 
 
Colonne 3 
Width 20 
Optional 
Hidden 
Class 
Class CDS 
From 1 
Tag Corresponding_CDS 
 
Colonne 4 
Width 20 
Optional 
Visible 
Class 
Class Gene 
From 3 
Tag Gene 

EOF

  print $qfhv $tablemaker_template2;
  return $tmdef2;
}


##############################################################
#
# Subroutines
#
##############################################################



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




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - Transposon_sync.pl

=head1 USAGE

=over 4

=item Transposon_sync.pl  [-options]

=back

This script queries bothe the gene reference database (geneace) and the species 
primary sequence database to make sure that curators have processed the data correctly 
in order to reflect the true status of the gene and annotations. Checks that the 
status is correct in geneace for a gene that has been flagged as "transposon_in_origin" 
storing these genes in an array for referenec when querying the sequence database. 
The second database query builds a list of all genes associated with a Transposon via 
TP::Gene and CDS::Gene. This list is then compared against the list of known Transposon 
genes from geneace. Any discrepancy in either query is reported.

Transposon_sync_test.pl MANDATORY arguments:

=over 4

=item None at present.

=back

Transposon_sync_test.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything....doesn't do anything in this script

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back


=head1 REQUIREMENTS

=over 4

=item That you are on the same file system as the soecies sequence databases and geneace

=back

=head1 AUTHOR

=over 4

=item Paul Davis (paul.davis@sanger.ac.uk)

=back

=cut
