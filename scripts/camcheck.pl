#!/usr/local/bin/perl5.6.1 -w
#
# camcheck.pl
#
# Cronjob integrity check controls for camace database.
#
# Usage: camcheck.pl
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2002-10-15 15:51:38 $
#
# see pod documentation (i.e. 'perldoc camcheck.pl') for more information.
#
##########################################################################################


#################################################################################
# variables                                                                     #
#################################################################################

$|=1;

use IO::Handle;
use Getopt::Std;
use Ace;
use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Socket;

##############################
# Script variables (run)     #
##############################

my $maintainers = "All";
my $rundate = `date +%y%m%d`; chomp $rundate;
my $runtime = `date +%H:%M:%S`; chomp $runtime;

##############################
# command-line options       #
##############################
use vars qw/ $opt_d $opt_h $opt_w $opt_m/;
getopts ("hdwm");
&usage(0) if ($opt_h);
my $debug = $opt_d;

##############################
# Paths etc                  #
##############################

my $clonepath = "/nfs/disk100/wormpub/analysis/cosmids";
my $clonefile = "$clonepath"."/current.versions";
my $tace      = glob("~wormpub/ACEDB/bin.ALPHA_4/tace");   # tace executable path
my $dbpath    = "/wormsrv1/camace";                           # Database path


# only tell Dan if running debug mode
$maintainers = "dl1\@sanger.ac.uk" if ($debug);


########################################
# Open logfile                         #
########################################

my $log="/wormsrv2/logs/camcheck.$rundate.$$";


open (LOG,">$log");
LOG->autoflush();

print LOG "# camcheck.pl\n";     
print LOG "# run details    : $rundate $runtime\n";
print LOG "\n";


########################################
# Connect with acedb server            #
########################################

my $db = Ace->connect(-path=>$dbpath,
		      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

print LOG "CamCheck run STARTED at $runtime\n\n";
print LOG "** Weekly edition **\n\n"  if ($opt_w);
print LOG "** Monthly edition **\n\n" if ($opt_m);


#####################################################################
# Do the following for every sequence in current.cosmid             #
#####################################################################

open (CLONEFILE,"<$clonefile") || die "Couldn't open $clonefile for reading\n";

my $line;
while(defined($line=<CLONEFILE>)){
  my $seq_file;
  my $seq_ace;
  chomp ($line);
  $line =~ m/(\w+)\/(\w+)/;
  my $clone = $1;
  my $dir_date = $2;
  
  print "[$clone|$dir_date]    \t" if ($debug);

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
    
  print " [FINISHED/ANNOTATED" if ($debug);
  &finannot($obj);

  ######################################################################
  # Check for N's in FINISHED sequences                                #
  ######################################################################
	
  print " | N's" if ($debug);
  &checkchars($obj,$seq_ace);

  ######################################################################
  # Compare date and checksum                                          #
  ######################################################################

  print " | DATE" if ($debug);
  &dateseq($obj,$dir_date);

  print " | CHKSUM" if ($debug);
  &chksum($seq_file,$seq_ace,$clone);

  ######################################################################
  # Check correctness of gene structure                                #
  ######################################################################

  print " | CDS_coords" if ($debug);
  &checkgenes($obj);

  ######################################################################
  # Check Sequence versions with EMBL                                  #
  ######################################################################

  if ($opt_w) {
      print " | Sequence versions" if ($debug);
      &check_sequence_version($obj);
  }
  
  ######################################################################
  # last check complete, tidy up                                       #
  ######################################################################

  print "]\n" if ($debug);

  ######################################################################
  # Get rid of this sequence object                                    #
  ######################################################################

  $obj->DESTROY();
}

close(CLONEFILE);


#########################
# Check predicted genes #
#########################

# need to close and reopen existing log filehandle either side of system call

close(LOG);

system("/wormsrv2/scripts/check_predicted_genes.pl -database $dbpath -log $log");
warn "check_predicted_genes.pl did not run correctly: $?" if ($?);

# now opened in append mode
open (LOG,">>$log");
LOG->autoflush();



############################################################
# Check that each clone only belongs to one (Sequence) map #
############################################################

my @multimap_clones= $db->fetch(-query=>'find Clone COUNT Map > 1');
if (@multimap_clones){
  foreach (@multimap_clones){
    print LOG "Clone error - $_ contains two or more map objects\n";
  }
}

###################
# LINK objects    #
###################

my $i = $db->fetch_many(-query=> 'find Sequence "SUPERLINK*"');  
while (my $obj = $i->next) {
  my $link = $obj;
  print "[$link]    \t" if ($debug);
  print " | CDS_coords" if ($debug);
  &checkgenes($obj);
  print "]\n" if ($debug);
}

$runtime = `date +%H:%M:%S`; chomp $runtime;
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

exit(0);

########################################################################################
####################################   Subroutines   ###################################
########################################################################################

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

sub finannot {    
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
    my $server = "srs.ebi.ac.uk";

    my $acc = $obj->Database(3);
    my $ver = $obj->Database(4);
    my ($EM_acc,$EM_seqver);

    my $querycontent="-e+([EMBLNEW-acc:'$acc'])|(([EMBL-acc:'$acc'])!(EMBL<EMBLNEW))";
    my $request = "/srsbin/cgi-bin/wgetz?$querycontent";

    if (!defined(open_TCP_connection(*F,$server,80))) {
	print "Error connecting to server at \n";
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
# Gene length as declared is subsequence and in exons list            #
#######################################################################

sub checkgenes {
  my $obj = shift;

  foreach my $child ($obj->Subsequence) {
    undef my @num;
    undef my ($method);
    undef my ($source);
	
    my ($seq, $start, $end) = $child->row();
    
    unless ($seq =~ /\./) {next;}
    
    my $diff = $end - $start;
    if ($diff < 0) {
      $diff = $start - $end;
    }
    
    my $subseq = $db->fetch(Sequence => "$child");
    if (!defined ($subseq)) {
      print LOG "Cannot fetch subsequence $child\n";
      next;
    }
    
    # Source
    #
    # All Subsequence objects must have a Source
    
    $source = $subseq->Source(1);
    if ((!defined ($source))) {
      print LOG "The subsequence $child has no Source\n";
    }
    
    # check to see if subsequence coordinates exceed length of parent clone
    # won't work for parents that are links (they have zero length)
    if($source !~ m/SUPERLINK/){
      my ($source_length) = $source->DNA(2);
      if (($start > $source_length) || ($end > $source_length)){
	print LOG "The subsequence $child has coordinates that exceed the length of its parent\n";
      }
    }
    
    # Source Exons
    #
    # All Subsequence objects must have Source_Exons
    
    @num = $subseq->Source_Exons(2);
    if (!@num) {
      print LOG "The subsequence $child has no Source_Exons\n";
    }
    my @sort = sort numerically @num;
    my $index = $#sort;
    my $length = ($sort[$index])-1;
    if ($diff != $length) {
      print LOG "The subsequence $child belonging to $obj has diff=$diff and length=$length\n";
    }
    
    
    # Method tags
    #
    # All Subsequence objects must have a Method
    
    $method = $subseq->Method(1);
    if ((!defined ($method))) {
      print LOG "The subsequence $child has no method\n";
    }
    
    # NDB_CDS || HALFWISE || gaze - don't process
    
    if ( ($method eq  "NDB_CDS") || ($method eq  "HALFWISE") || ($method eq  "gaze") ) {
      $subseq->DESTROY();
      $child->DESTROY();
      $diff="";
      $length="";
      next;
    }
    
    # Species
    # 
    # Most Subsequence objects must have a Species
    
    my $species = $subseq->Species(1);
    if ((!defined ($species))) {
      print LOG "The subsequence $child has no Species tag\n";
    }
    
    # Transcripts
    
    if ($child =~ /\S+\.t\d+/) {
      if (!defined ($subseq->at('Properties.Transcript'))) {
	print LOG "The subsequence $child has no Transcript tag\n";
      }
    }
    else {
      
      if ($method eq "Pseudogene") {
	if (!defined ($subseq->at('Properties.Pseudogene'))) {
	  print LOG "The subsequence $child [$method] has no Pseudogene tag\n";
	} 
      }
      if ( ($method eq "curated") || ($method eq "provisional") ) {
	if (!defined ($subseq->at('Properties.Coding.CDS'))) {
	  print LOG "The subsequence $child [$method] has no CDS tag\n";
	}
      }
    }
    
    
    #	print  "$child [$source|$diff|$length|$method]\n";
    
    $subseq->DESTROY();
    $child->DESTROY();
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
	push (@subsequence,$_) if (/^The subsequence/);
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
        print "The WormBase release number cannot be parsed\n";
#        print "Check File: '$Wormbase_release_file'\n\n";
        &run_details;
        exit(0);
    }
    elsif ($error == 0) {
        # Normal help menu
        exec ('perldoc',$0);
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

=item -d, Debug/Verbose mode

=item -w, Weekly checks are active

=item -m, Montly checks are active

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

=head3 Absence of Source in ?Sequence

All Subsequence Gene Models MUST have a parent ?Sequence.

    i.e. Sequence "ZK637.5"
         Source "ZK637"

=head3 Absence of Source_Exons in ?Sequence

All Subsequence Gene Models MUST have Source_Exons.

    i.e. Sequence "ZK637.5"
         Source_Exons      1   434
                         483   741
                         950  1159
                        1288  1413

=head3 Absence of Method in ?Sequence

All Subsequence Gene Models MUST have a Method tag. Method tags 
are used in two ways: Firstly, as a means of tagging Gene Models
into classes (e.g. 'curated' for active CDS prediction, 'RNA' for
RNA gene) and secondly as an internal ACEDB description of how
to display the object in the F-map.

    i.e. Sequence "ZK637.5"
         Method "curated"

For details of permissible Method tags see:
 http://intweb.sanger.ac.uk/Projects/C_elegans/MANUAL

=head3 Absence of Species tag in ?Sequence

All Subsequence Gene Models MUST have Species tag.

    i.e. Sequence "ZK637.5"
         Species "Caenorhabditis elegans"

=head3 Absence of CDS tag in ?Sequence

All Subsequence Gene Models MUST have a Coding CDS tag.

    i.e. Sequence "ZK637.5"
         Coding CDS

=head3 Mismatch between Parent coordinates and CDS span length

The coordinate span of the Gene Model (1 => n) based on the 
Source_Exon tags MUST be in agreement with the coordinate span
of the Gene Model as a Subsequence of the parent ?Sequence.

    i.e. Sequence "ZK637"
         Subsequence ZK637.5 11124 12536

         Sequence "ZK637.5"
         Source_Exons      1   434
                         483   741
                         950  1159
                        1288  1413

         Span in CDS (ZK637.5) = ( 1413 -     1) + 1 = 1413
         Parent (ZK637)        = (12536 - 11124) + 1 = 1413
                                                       ----

=head3 Subsequence coordinates exceeding parent sequence length

The coordinates of each subsequence must be less than the length of the
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

