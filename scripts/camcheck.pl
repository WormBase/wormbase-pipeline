#!/usr/local/bin/perl5.8.0 -w
#
# dbcheck.pl
#
# Cronjob integrity check controls for camace database.
#
# Usage: camcheck.pl
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2004-11-08 10:00:48 $
#
# see pod documentation (i.e. 'perldoc camcheck.pl') for more information.
#
#################################################################################


#################################################################################
# variables                                                                     #
#################################################################################

use IO::Handle;
use Getopt::Std;
use Ace;
use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Socket;

##############################
# Script variables (run)     #
##############################

my $maintainers = "All";

##############################
# command-line options       #
##############################

use vars qw/ $opt_d $opt_h $opt_w $opt_m $opt_s $opt_l $opt_e $opt_v/;
getopts ('hd:wms:le:v');

&usage(0) if ($opt_h); # perldoc help message

my $debug   = $opt_d;  # for debug options
my $verbose = $opt_v;  # for specifying more output

#  -h, Help
#  -d, Debug, specify user name to receive email
#  -w, Weekly checks are active
#  -m, Montly checks are active
#  -s, select which database to run against
#  -l, low level checks - not all the small intron gubbins in check_predicted_genes.pl
#  -e, Specifiy a mail recepient so that only the person responsible for a spilt database will be notified

##############################
# Paths etc                  #
##############################

my $tace    = &tace;                            # tace executable path
my $dbpath  = "/wormsrv1/camace";               # Database path

if ($opt_s) {
    $dbpath = $opt_s;
    &usage('bad database path') unless (-e "$dbpath/database/ACEDB.wrm");
}

# only email a specific person responsible for a database
if($opt_e){
  if(($opt_e eq "ar2") || ($opt_e eq "pad") || ($opt_e eq "dl1")){ 
     $maintainers = $opt_e;
  }
  else{
     $maintainers = "wormbase\@sanger.ac.uk";
  }
}

# Use debug mode?
if($debug){
  print "// Debug mode seleted, recipient is \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

# Verbose mode on?
print "// Verbose mode selected\n\n" if ($verbose);


########################################
# Open logfile                         #
########################################
my $dbname = $dbpath;
$dbname =~ s/.*camace(.*)/camace$1/;
my $rundate = &rundate;
my $runtime = &runtime;
my $log="/wormsrv2/logs/camcheck.$dbname.$rundate.$$";


open (LOG,">$log") || die "Couldn't write to log file\n";
LOG->autoflush();

print LOG "# camcheck.pl\n";     
print LOG "# run details $dbname  : $rundate $runtime\n";
print LOG "\n";


#########################################
# Connect with acedb server             #
#########################################

my $db = Ace->connect(-path=>$dbpath,
		      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

print LOG "CamCheck run STARTED at $runtime\n\n";
print LOG "** Weekly edition **\n\n"  if ($opt_w);
print LOG "** Monthly edition **\n\n" if ($opt_m);

#########################################
# Checks                                #
#########################################

&CloneTests;

#&check_worm_genes;

&check_worm_transcripts;

&CheckPredictedGenes;

&SingleSequenceMap;

&LinkObjects;

$runtime = &runtime;
print LOG "\nCamCheck run ENDED at $runtime\n\n";

close LOG;

##############################
# mail $maintainer report    #
##############################

&mail_maintainer("camcheck.pl Report:",$maintainers,$log);

##############################
# Write log to wormpub intweb
##############################

&writehtml;

##############################
# hasta luego                #
##############################

$db->close;
exit(0);

#################################################################################
################################   Subroutines   ################################
#################################################################################



#####################################################################
# Do the following for every sequence in current.cosmid             #
#####################################################################

sub CloneTests {

    my $clonepath = "/nfs/disk100/wormpub/analysis/cosmids";
    my $clonefile = "$clonepath"."/current.versions";
    
    open (CLONEFILE,"<$clonefile") || die "Couldn't open $clonefile for reading\n";
    
    my $line;
    while(defined($line=<CLONEFILE>)){
	my $seq_file;
	my $seq_ace;
	chomp ($line);
	$line =~ m/(\w+)\/(\w+)/;
	my $clone = $1;
	my $dir_date = $2;
	
	print "[$clone|$dir_date]    \t" if ($verbose);

	######################################################################
	# Retrieve the first sequence and date FROM DIRECTORY                #
	######################################################################
	
	my $seqpath = "$clonepath"."/"."$line"."/"."$clone.seq";
	
	open SEQPATH,"<$seqpath" || do {print LOG "$clone  \t:NOT_IN_DIRECTORY $seqpath\n"; next; };
	my $line1;
	while (defined($line1=<SEQPATH>)) {  
	    chomp($line1);
	    $seq_file.="$line1";
	}
	close SEQPATH;
	$seq_file =~ tr/a-z/A-Z/;
	$seq_file =~ s/\>\w+//;
	$seq_file =~ s/\W+//mg;
	if ($seq_file =~ /[^ACGTUMRWSYKVHDBXN]/img) {
	    print LOG  "DIRSEQ for $clone contains bad characters\n";
	    $seq_file =~ s/[^ACGTUMRWSYKVHDBXN]//img;
	}
	
	######################################################################
	# Retrieve the second sequence and date FROM ACEDB                   #
	######################################################################
	
	my $obj = $db->fetch(Sequence=>$clone);
	if (!defined ($obj)) {
	    print LOG "Could not fetch sequence $clone\n";
	    next;
	}
	
	#####################################################################
	# skipping the non-canonical genomic sequences                      #
	#####################################################################
	
	my $canonical=$obj->Properties(1);
	if ($canonical !~ /Genomic_canonical/) {
	    print LOG "Not Genomic_canonical sequence $clone\n";
	    next;
	}
  
	#####################################################################
	# Push the sequence as string in $seq2                              #
	#####################################################################
	
	$seq_ace=$obj->asDNA();
	if (!$seq_ace) {
	    print LOG "$clone NOT_IN_ACEDB $clone\n" ;
	    next;
	}
	$seq_ace =~ s/\>\w+//mg;
	$seq_ace =~ tr/a-z/A-Z/;
	$seq_ace =~ s/\W+//mg;
	
	
	######################################################################
	# Iterative checks for each clone                                    #
	######################################################################
	
	######################################################################
	# ?Sequence is Finished but not annotated                            #
	######################################################################
	
	print " [FINISHED/ANNOTATED" if ($verbose);
	&CloneStatus($obj);
	
	######################################################################
	# Check for N's in FINISHED sequences                                #
	######################################################################
	
	print " | N's" if ($verbose);
	&checkchars($obj,$seq_ace);
	
	######################################################################
	# Compare date and checksum                                          #
	######################################################################
	
	print " | DATE" if ($verbose);
	&dateseq($obj,$dir_date);
	
	print " | CHKSUM" if ($verbose);
	&chksum($seq_file,$seq_ace,$clone);
	
	######################################################################
	# Check correctness of gene structure                                #
	######################################################################
	
	print " | CDS_coords" if ($verbose);
	&check_CDSs($obj);
	
	######################################################################
	# Check Sequence versions with EMBL                                  #
	######################################################################
	
	if ($opt_w) {
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
    
    close(CLONEFILE);
}



####################################
# Check worm_genes composite class #
####################################

sub check_worm_genes{

  # look for objects with no Gene tag
  my @genes= $db->fetch(-query=>'find worm_genes NOT Gene');
    if(@genes){
      foreach (@genes){
	print LOG "ERROR: $_ has no Gene tag, please add valid Gene ID from geneace\n";
      }
    }
}

####################################
# Check elegans_RNA_gene class     #
####################################

sub check_worm_transcripts{
  my @Transcripts= $db->fetch(-query=>'find elegans_RNA_genes NOT Transcript');
  if(@Transcripts){
    foreach (@Transcripts){
      print LOG "ERROR: $_ has no Transcript tag, this will cause errors in the build\n";
    }
  }
}
#########################
# Check predicted genes #
#########################

sub CheckPredictedGenes {
    
    # need to close and reopen existing log filehandle either side of system call

    close(LOG);

    my $cpg_call = "/wormsrv2/scripts/check_predicted_genes.pl -database $dbpath -log $log";
    $cpg_call .= " -basic" if $opt_l;
    system("$cpg_call");
    warn "check_predicted_genes.pl did not run correctly: $?" if ($?);
    
# now opened in append mode
    open (LOG,">>$log");
    LOG->autoflush();
    
}


############################################################
# Check that each clone only belongs to one (Sequence) map #
############################################################

sub SingleSequenceMap {

    my @multimap_clones= $db->fetch(-query=>'find Clone COUNT Map > 1');
    if (@multimap_clones){
	foreach (@multimap_clones){
	    print LOG "Clone error - $_ contains two or more map objects\n";
	}
    }
}

###################
# LINK objects    #
###################

sub LinkObjects {

    my $i = $db->fetch_many(-query=> 'find Sequence "SUPERLINK*"');  
    while (my $obj = $i->next) {
	my $link = $obj;
	print "[$link]    \t" if ($verbose);
	print " | CDS_coords" if ($verbose);
	&check_CDSs($obj);
	print "]\n" if ($verbose);
    }
}



####################################
# Coherency check between directory and database
####################################

sub dateseq {
  my $obj = shift;
  my $dir_date = shift;
  my $ace_date = $obj->Date_directory(1);
        
  if ($dir_date != $ace_date) {
    print LOG "DATE mismatch in $obj; dir $dir_date acedb $ace_date\n";
  }
} 

   
################################################
# Coherency check between directory and database
################################################

sub chksum {
    my ($seq_file,$seq_ace,$clone) = @_;
    my ($checksum, $index, $char);

    # calculate checksum routines from Chao-Kung

    my $chk1 = &calculate_chksum($seq_file);
    my $chk2 = &calculate_chksum($seq_ace);

    if ($chk1 != $chk2) {
	print LOG "SEQUENCE mismatch in $clone; dir $chk1 acedb $chk2\t";
	print LOG "=> dir: " . length ($seq_file) . " ace: " . length ($seq_ace) . "\n";

    }
}


sub calculate_chksum {
    my $seq = shift;
    my ($checksum, $index, $char,$chk);
    
    $index = 0;
    
    foreach $char ( split(/[\.\-]*/, $seq)) {
        $index++;
        $checksum += ($index * (unpack("c",$char) || 0) );
        if( $index ==  57 ) {
            $index = 0;
        }
    }

    $chk = $checksum % 10000;
    $checksum=(); 

    return ($chk);
}


####################################
# Finished / Annotated
####################################

sub CloneStatus {
    
    my $obj =shift;
    my $finished  = $obj->Finished(1);
    my $annotated = $obj->Annotated(1);

    if (!$finished) {
	print LOG "NOT_FINISHED $obj\n";
    }    
    if (($finished)&&(!$annotated)){
	print LOG "FINISHED_BUT_NOT_ANNOTATED $obj\n";
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
	print LOG "ACEDBSEQ for $obj contains bad characters\n";
	$seq_ace =~s/[^ACGTUMRWSYKVHDBXN]//img;
    }
    if (($seq_ace =~ /N/g) && ($finished)) { 
	print LOG "ACEDBSEQ FINISHED SEQUENCE for $obj contains N \n";
    }
}

#######################################################################
# Check sequence version between camace and EMBL                      #
#######################################################################

sub check_sequence_version {
    local (*F);
    my $obj = shift;
    my $server = "srs.sanger.ac.uk";

    # get details from camace
    my ($acc) = $obj->at('DB_info.Database.EMBL.NDB_AC');
    my ($ver) = $obj->at('DB_info.Database.EMBL.NDB_SV');

    # get details from EMBL
    my ($EM_acc,$EM_seqver);
    my $querycontent="-e+([EMBLNEW-acc:'$acc'])|(([EMBL-acc:'$acc'])!(EMBL<EMBLNEW))";
    my $request = "/srsbin/cgi-bin/wgetz?$querycontent";

    if (!defined(open_TCP_connection(*F,$server,80))) {
	print "Error connecting to srs.sanger.ac.uk server for object $obj\n";
	exit(-1);
    }
    print F "GET $request HTTP/1.0\n\n";
    print F "Accept: */*\n";
    print F "User-Agent: socketsrs/1.0\n\n";
 
    while (my $return_line=<F>) {
#    print "$return_line";
	if ($return_line =~ /^SV\s+(\S+)\.(\d+)/) {
	    ($EM_acc,$EM_seqver) = ($1,$2);
	    last;
	}
    }
    
    if ($EM_seqver != $ver) {
	print LOG "Sequence version discrepent $obj [camace = $ver | EMBL = $EM_seqver]\n";
    }

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
# Gene length as declared is CDS and in exons list            #
#######################################################################

sub check_CDSs {
  my $obj = shift;

  foreach my $cds ($obj->CDS_child) {
    undef my @num;
    undef my ($method);
    undef my ($parent);
	
    my ($seq, $start, $end) = $cds->row();
    
    unless ($seq =~ /\./) {next;}
    
    my $diff = $end - $start;
    if ($diff < 0) {
      $diff = $start - $end;
    }
    
    my $subseq = $db->fetch(CDS => "$cds");
    if (!defined ($subseq)) {
      print LOG "Cannot fetch CDS $cds\n";
      next;
    }
    
    # All CDS objects must have a 'Sequence' tag to connect to parent object    
    $parent = $subseq->Sequence;
    if ((!defined ($parent))) {
      print LOG "The CDS $cds has no Sequence tag\n";
    }
    
    # check to see if CDS coordinates exceed length of parent clone
    # won't work for parents that are links (they have zero length)
    if($parent !~ m/SUPERLINK/){
      my ($parent_length) = $parent->DNA(2);
      if (($start > $parent_length) || ($end > $parent_length)){
	print LOG "The CDS $cds has coordinates that exceed the length of its parent\n";
      }
    }
    
    # Source Exons
    #
    # All CDS objects must have Source_exons
    
    @num = $subseq->Source_exons(2);
    if (!@num) {
      print LOG "The CDS $cds has no Source_exons\n";
    }
    my @sort = sort numerically @num;
    my $index = $#sort;
    my $length = ($sort[$index])-1;
    if ($diff != $length) {
      print LOG "The CDS $cds belonging to $obj has diff=$diff and length=$length\n";
    }
    
    
    # Method tags
    #
    # All CDS objects must have a Method
    
    $method = $subseq->Method(1);
    if ((!defined ($method))) {
      print LOG "The CDS $cds has no method\n";
    }
    
    # NDB_CDS || HALFWISE || gaze - don't process
    
    if ( ($method eq  "NDB_CDS") || ($method eq  "HALFWISE") || ($method eq  "gaze") ) {
      $subseq->DESTROY();
      $cds->DESTROY();
      $diff="";
      $length="";
      next;
    }
    
    # Species
    # 
    # CDS objects must have a Species
    
    my $species = $subseq->Species(1);
    if ((!defined ($species))) {
      print LOG "The CDS $cds has no Species tag\n";
    }
 

# The following section of code won't work now and needs to be revamped to actually look at 
# ?Transcript and ?Pseudogene objects.
   
    # Transcripts
    
#    if ($cds =~ /\S+\.t\d+/) {
#      if (!defined ($subseq->at('Properties.Transcript'))) {
#	print LOG "The CDS $cds has no Transcript tag\n";
#      }
#    }
#    else {
#      
#      if ($method eq "Pseudogene") {
#	if (!defined ($subseq->at('Properties.Pseudogene'))) {
#	  print LOG "The CDS $cds [$method] has no Pseudogene tag\n";
#	} 
#      }
#      if ( ($method eq "curated") || ($method eq "provisional") ) {
#	if (!defined ($subseq->at('Properties.Coding.CDS'))) {
#	  print LOG "The subsequence $cds [$method] has no CDS tag\n";
#	}
#      }
#    }
    
    
    #	print  "$cds [$parent|$diff|$length|$method]\n";
    
    $subseq->DESTROY();
    $cds->DESTROY();
    $diff="";
    $length="";
    
  }
}


sub numerically {
  $a <=> $b;
}


#######################################################################
# Write HTML page with maintenance job results
#######################################################################

sub writehtml {

    my (@finished);
    my (@annotated);
    my (@date);
    my (@sequence);
    my (@subsequence);
    my (@link);
    my (@clone);
    my (@seqversion);

    my $logdir="/nfs/disk100/wormpub/LocalWWW";

my $HTML_START=<<START;
<HTML>
<HEAD>
<TITLE>Camace Automated db Maintenance log</TITLE>
</HEAD>
<BODY BGCOLOR="WHITE">
START

my $HTML_END=<<END;
</BODY>
</HTML>
END

    open (OUTHTML,">$logdir/camchecklog.html") || die "Couldn't open $logdir/camchecklog.html for writing\n";
    print OUTHTML $HTML_START;
    print OUTHTML "Last ran  : <B>$rundate $runtime</B></P>\n";
    print OUTHTML "<TABLE BORDER=1 WIDTH=100%>\n";

    open (READLOG, "<$log") || die "Couldn't open $log for reading\n";
    while (<READLOG>) {
	push (@finished,$_)    if (/^NOT_FINISHED/);
	push (@annotated,$_)   if (/^FINISHED_BUT_NOT_ANNOTATED/);
	push (@date,$_)        if (/^DATE/);
	push (@sequence,$_)    if (/^SEQUENCE/);
	push (@seqversion,$_)  if (/^Sequence version/);
	push (@subsequence,$_) if (/^The CDS/);
	push (@subsequence,$_) if (/^Gene error - /);
	push (@link,$_)        if (/^SUPERLINK/);
	push (@clone,$_)       if (/^Clone error - /);
    }
    close READLOG;

    # Not Finished
    print OUTHTML "<TABLE BORDER=\"0\" WIDTH=\"100%\">\n";
    print OUTHTML "<TH BGCOLOR=\"darkblue\"><TD><FONT COLOR=\"white\">Not Finished</FONT></TD></TH>\n";
    foreach (@finished) {
	print OUTHTML "<TR><TD><FONT SIZE=\"-1\">$_</FONT></TD></TR>\n";
    }
    print OUTHTML "</TABLE></P>\n";

    # Not Annotated
    print OUTHTML "<TABLE BORDER=\"0\" WIDTH=\"100%\">\n";
    print OUTHTML "<TH BGCOLOR=\"darkblue\"><TD><FONT COLOR=\"white\">Not Annotated</FONT></TD></TH>\n";
    foreach (@annotated) {
	print OUTHTML "<TR><TD><FONT SIZE=\"-1\">$_</FONT></TD></TR>\n";
    }
    print OUTHTML "</TABLE></P>\n";

    # Date Mismatch (not uploaded into camace)
    print OUTHTML "<TABLE BORDER=\"0\" WIDTH=\"100%\">\n";
    print OUTHTML "<TH BGCOLOR=\"darkblue\"><TD><FONT COLOR=\"white\">Date mismatches</FONT></TD></TH>\n";
    foreach (@date) {
	print OUTHTML "<TR><TD><FONT SIZE=\"-1\">$_</FONT></TD></TR>\n";
    }
    print OUTHTML "</TABLE></P>\n";

    # Sequence Mismatch (problem)
    print OUTHTML "<TABLE BORDER=\"0\" WIDTH=\"100%\">\n";
    print OUTHTML "<TH BGCOLOR=\"darkblue\"><TD><FONT COLOR=\"white\">Genomic sequences</FONT></TD></TH>\n";
    foreach (@sequence) {
	print OUTHTML "<TR><TD><FONT SIZE=\"-1\">$_</FONT></TD></TR>\n";
    }
    print OUTHTML "<TR>&nbsp;</TR>\n";

    # Sequence version Mismatch (problem)
    foreach (@seqversion) {
	print OUTHTML "<TR><TD><FONT SIZE=\"-1\">$_</FONT></TD></TR>\n";
    }
    print OUTHTML "</TABLE></P>\n";

    # Subsequence problem
    print OUTHTML "<TABLE BORDER=\"0\" WIDTH=\"100%\">\n";
    print OUTHTML "<TH BGCOLOR=\"darkblue\"><TD><FONT COLOR=\"white\">Subsequences</FONT></TD></TH>\n";
    foreach (@subsequence) {
	print OUTHTML "<TR><TD><FONT SIZE=\"-1\">$_</FONT></TD></TR>\n";
    }
    print OUTHTML "</TABLE></P>\n";

    # SUPERLINK problems
    print OUTHTML "<TABLE BORDER=\"0\" WIDTH=\"100%\">\n";
    print OUTHTML "<TH BGCOLOR=\"darkblue\"><TD><FONT COLOR=\"white\">LINK objects</FONT></TD></TH>\n";
    foreach (@link) {
	print OUTHTML "<TR><TD><FONT SIZE=\"-1\" COLOR=\"red\">$_</FONT></TD></TR>\n";
    }
    print OUTHTML "</TABLE></P>\n";

     # clone problems
    print OUTHTML "<TABLE BORDER=\"0\" WIDTH=\"100%\">\n";
    print OUTHTML "<TH BGCOLOR=\"darkblue\"><TD><FONT COLOR=\"white\">Clone objects</FONT></TD></TH>\n";
    foreach (@clone) {
	print OUTHTML "<TR><TD><FONT SIZE=\"-1\" COLOR=\"red\">$_</FONT></TD></TR>\n";
    }
    print OUTHTML "</TABLE></P>\n";


    print OUTHTML $HTML_END;
    close OUTHTML;

    undef (@finished);
    undef (@annotated);
    undef (@date);
    undef (@sequence);
    undef (@seqversion);
    undef (@subsequence);
    undef (@link);
    undef (@clone);
}

#######################################################################
# Help and error trap outputs                                         #
#######################################################################

sub run_details {
    print "# camcheck.pl\n";     
    print "# run details    : $rundate $runtime\n";
    print "\n";
}


sub usage {
    my $error = shift;

    if ($error == 1) {
        # No WormBase release number file
        print "The WormBase release number cannot be parsed\n\n";
        &run_details;
        exit(0);
    }
    elsif ($error == 0) {
        # Normal help menu
        exec ('perldoc',$0);
    }
    elsif ($error eq "bad database path") {
	print "The database path is invalid. No ACEDB.wrm file found at $opt_s\n\n";
        &run_details;
        exit(0);
    }

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
all Genome_sequences and SUPERLINK* objects.

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

=item -s, Select a specific database to test eg -s ~wormpub/camace_ar2

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

