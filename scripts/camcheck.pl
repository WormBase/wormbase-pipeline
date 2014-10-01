#!/software/bin/perl -w
#
# dbcheck.pl
#
# Cronjob integrity check controls for camace database.
#
# Usage: camcheck.pl
#
# Last updated by: $Author: pad $
# Last updated on: $Date: 2014-10-01 15:43:11 $
#
# see pod documentation (i.e. 'perldoc camcheck.pl') for more information.
#
#################################################################################


#################################################################################
# variables                                                                     #
#################################################################################

use IO::Handle;
use Getopt::Long;
use Ace;
use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Socket;
use Log_files;
use Storable;

##############################
# command-line options       #
##############################
my ($help,$verbose,$debug,$test,$Weekly,$Montly,$database,$Low,$email,$store,$species);

GetOptions(
	   'help'       => \$help,    #  -help
	   'verbose'    => \$verbose, #  -verbose option.
	   'debug:s'    => \$debug,   #  -debug, debug option
	   'test'       => \$test,    #  -test, TEST_BUILD env used.
	   'w'          => \$Weekly,  #  -w, Weekly checks are active
	   'm'          => \$Montly,  #  -m, Montly checks are active
	   'database:s' => \$database,#  -database select which database to run against
	   'low'        => \$Low,     #  -low level checks - not all the small intron gubbins
	   'email:s'    => \$email,   #  -email, Specifiy a mail recepient so that only the person responsible for a spilt database will be notified
	   'store:s'    => \$store,
	   'species:s'  => \$species,  #  -species 'caenorhabditis elegans' used for limiting the genome sequences retrieved
	  );


my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}

##############################
# Paths etc                  #
##############################
unless (defined$species) {$species = "Caenorhabditis elegans";}
my $log = Log_files->make_build_log($wormbase);
my $tace    = $wormbase->tace;                            # tace executable path
my $dbpath  = $wormbase->database('camace');              # db path
my $maintainers;
my $dbname = $dbpath;
$dbname =~ s/.*camace(.*)/camace$1/;
my $rundate = $wormbase->rundate;
my $runtime = $wormbase->runtime;

if ($database) {
  $dbpath = $database;
  unless (-e "$dbpath/database/ACEDB.wrm") {$log->log_and_die("Failed to connect to $dbpath\n")};
}
print "\n$dbpath is being checked.......\n\n" if $debug;

###########################################################
# only email a specific person responsible for a database #
###########################################################

if($email){
  if(($email eq "gw3") || ($email eq "pad")){ 
    $maintainers = "$email\@sanger.ac.uk";
    print "Email recipient \= $maintainers\n\n" if $verbose;
  }
  else{
    $maintainers = "wormbase\@sanger.ac.uk";
  }
}

# Use debug mode?
if($debug){
  print "// Debug mode selected, recipient is \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

# Verbose mode on?
print "// Verbose mode selected\n\n" if ($verbose);
$log->write_to("# camcheck.pl\n");     
$log->write_to("# run details $dbname  : $rundate $runtime\n");
$log->write_to("\n");


#########################################
# Connect with acedb server             #
#########################################

my $db = Ace->connect(-path=>$dbpath) or $log->log_and_die("Connection failure: ".Ace->error());

$log->write_to("CamCheck run STARTED at $runtime\n\n");
$log->write_to("** Weekly edition **\n\n") if ($Weekly);
$log->write_to("** Monthly edition **\n\n") if ($Montly);

#########################################
# Checks                                #
#########################################

&CloneTests;

&check_worm_transcripts;

&CheckPredictedGenes;

&DBBackupandCompare;

&SingleSequenceMap;

$log->write_to("\nCamCheck run ENDED at $runtime\n\n");

$db->close;
$log->mail;
exit(0);

#################################################################################
################################   Subroutines   ################################
#################################################################################



#####################################################################
# Do the following for every sequence in current.cosmid             #
#####################################################################

sub CloneTests {

  my $seq_ace;
  my @clones;
#  my @clones =  $db->fetch(
  my $clone;
  foreach $clone (@clones) {
    
    
    my $obj = $db->fetch(Sequence=>$clone);
    if (!defined ($obj)) {
      $log->write_to("Could not fetch sequence $clone\n");
      next;
    }
    
    #####################################################################
    # skipping the non-canonical genomic sequences                      #
    #####################################################################
    
    my $canonical=$obj->Properties(1);
    if ($canonical !~ /Genomic_canonical/) {
      $log->write_to("Not Genomic_canonical sequence $clone\n");
      next;
    }
    
    #####################################################################
    # Push the sequence as string in $seq2                              #
    #####################################################################
    
    $seq_ace=$obj->asDNA();
    if (!$seq_ace) {
      $log->write_to("$clone NOT_IN_ACEDB $clone\n");
      next;
    }
    $seq_ace =~ s/\>\w+//mg;
    $seq_ace =~ tr/a-z/A-Z/;
    $seq_ace =~ s/\W+//mg;
    
    
    ######################################################################
    # Iterative checks for each clone                                    #
    ######################################################################
    
    &check_overlap_right($obj);
    
    ######################################################################
    # Check for N's in FINISHED sequences                                #
    ######################################################################
    
    print " | N's" if ($verbose);
    &checkchars($obj,$seq_ace);	
    
    ######################################################################
    # Check Sequence versions with EMBL                                  #
    ######################################################################
    
    if ($Weekly) {
      print " | Sequence versions" if ($verbose);
      &check_sequence_version($obj);
    }
    
    ######################################################################
    # last check complete, tidy up                                       #
    ######################################################################
    
    print "]\n" if ($verbose);
    
    ######################################################################
    # Get rid of this sequence object                                    #
    ######################################################################
    
    $obj->DESTROY();
  }
}


####################################
# Check elegans_RNA_gene class     #
####################################

sub check_worm_transcripts{
  my @Transcripts= $db->fetch(-query=>'find elegans_RNA_genes NOT Transcript');
  if(@Transcripts){
    foreach (@Transcripts){
      $log->write_to("ERROR: $_ has no Transcript tag, this will cause errors in the build\n");
    }
  }
  else {
    $log->write_to("\nTranscripts OK\n");
  }
}
#########################
# Check predicted genes #
#########################

sub CheckPredictedGenes {
  $log->write_to("runnning CheckPredictedGenes\n");
  if ($Low) {
    $wormbase->run_script("check_predicted_genes.pl -database $dbpath -basic", $log);
  }
  else {
    $wormbase->run_script("check_predicted_genes.pl -database $dbpath" , $log);
  }
}


#########################
# Check for lost data   #
#########################

sub DBBackupandCompare {
  $log->write_to("runnning DB_backup_and_compare\n");
  $wormbase->run_script("db_backup_and_compare.pl -db camace" , $log);
}

################################################
# Check that all overlap_right data is in place.
################################################
sub check_overlap_right {
  my $obj = shift;
  if (!($obj->name =~ /cTel/)) {
    $log->write_to("No overlap_right info found for $obj!!! This will cause errors in the build.\n") unless defined($obj->at('Structure.Overlap_right'));
  }
}

############################################################
# Check that each clone only belongs to one (Sequence) map #
############################################################

sub SingleSequenceMap {

    my @multimap_clones= $db->fetch(-query=>'find Clone COUNT Map > 1');
    if (@multimap_clones){
	foreach (@multimap_clones){
	    $log->write_to("Clone error - $_ contains two or more map objects\n");
	}
    }
    else {
      $log->write_to("no clones found with more than one Map\n");
    }
  }


#######################################################################
# Odd chars and N's in  finished sequences                            #
#######################################################################
sub checkchars {
  my $obj = shift;
  my $seq_ace = shift;
  my $finished  = $obj->Finished(1);
  
    if ($seq_ace =~ /[^ACGTUMRWSYKVHDBXN]/img) {
	$log->write_to("ACEDBSEQ for $obj contains bad characters\n");
	$seq_ace =~s/[^ACGTUMRWSYKVHDBXN]//img;
    }
    if (($seq_ace =~ /N/g) && ($finished)) { 
	$log->write_to("ACEDBSEQ FINISHED SEQUENCE for $obj contains N \n");
    }
}

########################################################################
# Check sequence version between camace and EMBL via mfetch libraries. #
########################################################################
sub check_sequence_version {
my ($EMBL_acc,$EMBL_sv,$ACE_ver);
  my $obj = shift;
  
  # get details from camace
  my ($acc) = $obj->at('DB_info.Database.EMBL.NDB_AC');
  my ($ver) = $obj->at('DB_info.Database.EMBL.NDB_SV');
  
if ($ver =~ /\S+.(\d+)/) {
  $ACE_ver = $1;
}
# fetch the ID line for the most rescent acc, not version exclusive.
open (SEQUENCE, "/usr/local/pubseq/scripts/mfetch -d embl -f id  -i \"sv:$acc*\" |");

while (<SEQUENCE>) {
  # deal with header line
  chomp;
  
  #if the mfetch query fails you get "no match"
  if (/no match/){
    $log->write_to("$acc could not be retrieved from EMBL.\n");
    next;
  }
  #ID   Z19152; SV 1; linear; genomic DNA; STD; INV; 40909 BP.
  if (/ID   (\S+); SV (\d+);/) {
    ($EMBL_acc,$EMBL_sv) = ($1,$2);
    next;
  }
}
if (($EMBL_sv = $ACE_ver) && ($verbose)) {print "Sequence version matches for $acc [camace = $ACE_ver | EMBL = $EMBL_sv]\n";}
if ($EMBL_sv != $ACE_ver) {$log->write_to("WARNING: Sequence version discrepent $obj [camace = $ACE_ver | EMBL = $EMBL_sv]\n");}
close SEQUENCE;
}

########################################################
# Output: successful network connection in file handle #
########################################################

sub open_TCP_connection  {
    my ($FS,$dest,$port) = @_;
    my $proto = getprotobyname ('tcp');
    socket ($FS,PF_INET,SOCK_STREAM,$proto);
    my $sin = sockaddr_in($port,inet_aton($dest));
    connect($FS,$sin) || return undef;
    my $old_fh = select($FS);
    $| = 1;
    select($old_fh);
}



#######################################################################
# Write HTML page with maintenance job results
#######################################################################


#######################################################################
# Help and error trap outputs                                         #
#######################################################################

sub run_details {
    print "# camcheck.pl\n";     
    print "# run details    : $rundate $runtime\n";
    print "\n";
}



__END__

=pod

=head1 NAME - camcheck.pl

=back

=head1 USAGE

=over 4

=item camcheck.pl [-options]

camcheck.pl performs a number of integrity/consistency checks against
the camace database. The script is based on an iterative loop across
all Genome_sequences objects.

=back

=head2 camcheck.pl MANDATORY arguments:

=over 4

=item none

=back

=head2 camcheck.pl OPTIONAL arguments:

=over 4

=item -h, Help

=item -d <username>, debug mode, sends email to user only

=item -v, verbose - writes longer output to screen

=item -w, Weekly checks are active

=item -m, Montly checks are active

=item -s, Select a specific database to test eg -s ~wormpub/camace_pad

=item -l, Doesnt do all the small intron checking stuff in check_predicted_genes.pl

=item -e, Only emails person responsible for specific split database

=back

=head1 DOCUMENTATION

=over 4

=back

The following checks have been incorporated into camcheck.pl:

=head2 Status tags

=head3 Genome sequences which are not Finished.

Genome sequences which do not have a Finished tag.

=head3 Genome sequences which are Finished but not Annotated.

Genome sequences which are finished but not annotated.

=head2 File storage on /analysis/cosmids

=head3 Date mismatch between the file system and camace.

Inconsistent Date_directory tag in ACEDB with respect to the file
system (/analyis/cosmids/current.versions).

For details of how the date dirtectory structure works:
 http://intweb.sanger.ac.uk/Projects/C_elegans/MANUAL

=head3 Sequence mismatch between the file system and camace. 

This is based on a GCG checksum calculation for the .seq file in 
the date directory and the sequence extracted from ACEDB.

For details of how the date dirtectory structure works:
 http://intweb.sanger.ac.uk/Projects/C_elegans/MANUAL

=head2 Sequence version synchrony with EMBL (Weekly)

=head3 Sequence version (SV) field in EMBL files should be in
synchrony with the DB_info Database fields.

=head2 Gene Models

=head3 Absence of Sequence tag in ?CDS

All CDS Gene Models MUST have a parent ?Sequence.

    i.e. CDS "ZK637.5"
         Sequence "ZK637"

=head3 Absence of Source_exons in ?Sequence

All CDS Gene Models MUST have Source_exons.

    i.e. CDS "ZK637.5"
         Source_exons      1   434
                         483   741
                         950  1159
                        1288  1413

=head3 Absence of Method in ?CDS

All CDS Gene Models MUST have a Method tag. Method tags 
are used in two ways: Firstly, as a means of tagging Gene Models
into classes (e.g. 'curated' for active CDS prediction, 'RNA' for
RNA gene) and secondly as an internal ACEDB description of how
to display the object in the F-map.

    i.e. CDS "ZK637.5"
         Method "curated"

For details of permissible Method tags see:
 http://intweb.sanger.ac.uk/Projects/C_elegans/MANUAL

=head3 Absence of Species tag in ?CDS

All CDS Gene Models MUST have Species tag.

    i.e. CDS "ZK637.5"
         Species "Caenorhabditis elegans"

=head3 Absence of CDS tag in ?CDS

All CDS Gene Models MUST have a Coding CDS tag.

    i.e. CDS "ZK637.5"
         Coding CDS

=head3 Mismatch between Parent coordinates and CDS span length

The coordinate span of the Gene Model (1 => n) based on the 
Source_exon tags MUST be in agreement with the coordinate span
of the Gene Model as a CDS of the parent ?Sequence.

    i.e. Sequence "ZK637"
         CDS ZK637.5 11124 12536

         CDS "ZK637.5"
         Source_exons      1   434
                         483   741
                         950  1159
                        1288  1413

         Span in CDS (ZK637.5) = ( 1413 -     1) + 1 = 1413
         Parent (ZK637)        = (12536 - 11124) + 1 = 1413
                                                       ----

=head3 CDS coordinates exceeding parent sequence length

The coordinates of each CDS must be less than the length of the
parent sequence.

=head2 Valid predicted gene objects

camcheck.pl calls the check_predicted_genes.pl script which performs a
number of integrity checks on each object within the 'Predicted_gene'
class. See L<check_predicted_genes.pl> for more details of these 
checks.

=head2 Valid clone/map relationships

If clones are placed on a sequence map then they should exist only on one
sequence map and not on two different maps belonging to different chromosomes.
camcheck.pl determines if this is the case by using the simple acedb query:
'find Clone COUNT Map > 1'


=head1 AUTHOR - Daniel Lawson (with some contributions by Keith Bradnam)

Email dl1@sanger.ac.uk, krb@sanger.ac.uk

=cut

