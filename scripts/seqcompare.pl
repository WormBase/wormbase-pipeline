#!/usr/local/bin/perl

# Seqcompare.pl - Maintenance Crontab Job for CamAce
# Retrieves sequences from ACEDB; compares them with the one
# contained in wormpub directory; checks sequence content and
# date correctness et al; reports errors by email
# Kills himself if hung up for > 5000 sec
# Ver. 1.00 A.Guffanti, August 1999 uses BioPerl, AcePerl

BEGIN {
  unshift (@INC,"/nfs/disk92/PerlSource/Bioperl/Releases/bioperl-0.05");
}
use Bio::Seq;
use Ace;

$|=1;
$SIG{'ALRM'} = \&timed_out;
@time = localtime();
($MINS,$HOURS,$DAY,$MONTH,$YEAR)=(localtime)[1,2,3,4,5];
$REALMONTH=$MONTH+1;
$REALYEAR=$YEAR+1900;
$TODAY = "$DAY $REALMONTH $REALYEAR at $HOURS:$MINS";

$seqerror="Job report from SeqCompare run $$ STARTED $TODAY\n\n";
$acepath="/nfs/disk100/wormpub/acedb/ace4/cam";
$progrpath="/nfs/disk100/acedb/RELEASE.SUPPORTED/bin.ALPHA_4";
$clonepath="/nfs/disk100/wormpub/analysis/cosmids";
$clonefile="$clonepath"."/current.versions";
$outfile="/tmp/seqcompare.$$";
$db = Ace->connect(-path=>$acepath);

open (CLONEFILE,"<$clonefile");
eval {
  alarm(5000);
  while($line=<CLONEFILE>) {
    $seq1=$seq2="";
    chomp ($line);
    $line =~ m/(\w+)\/(\w+)/;
    $seqname=$1;
    # Retrieve the first sequence and date FROM REPOSITORY DIRECTORY and push the sequence as string in $seq1
    $d_date=$2;
    $seqpath="$clonepath"."/"."$line"."/"."$seqname.seq";
    open SEQPATH,"<$seqpath" || do {$seqerror .= "NOT_IN_DIRECTORY $seqpath\n" ;next;};
    while ($line1=<SEQPATH>) {  
      chomp($line1);
      $seq1.="$line1";
    }
   close SEQPATH;
#   if (!$seq1) {
#      $seqerror .= "NOT_IN_DIRECTORY $seqname\n" ;
#    }
    $seq1=~tr/a-z/A-Z/;
    $seq1=~s/\>\w+//;
    $seq1=~s/\W+//mg;
    if ($seq1 =~ /[^ACGTUMRWSYKVHDBXN]/img) {
      $seqerror .= "DIRSEQ for $seqname contains bad characters\n";
      $seq1=~s/[^ACGTUMRWSYKVHDBXN]//img;
    }
    # Retrieve the second sequence and date FROM ACEDB and push the sequence as string in $seq2
    $obj = $db->fetch(Sequence=>$seqname);
    $finished=$obj->Finished(1);
    $annotated=$obj->Annotated(1);
    if (!$finished) {
      $seqerror .= "NOT_FINISHED $seqname\n";
    }    
    if (($finished)&&(!$annotated)){
      $seqerror .= "FINISHED_BUT_NOT_ANNOTATED $seqname\n";
    }
    $seq2=$obj->asDNA();
    if (!$seq2) {
      $seqerror .= "NOT_IN_ACEDB $seqname\n" ;
      next;
    }
    $seq2=~s/\>\w+//mg;
    $seq2=~tr/a-z/A-Z/;
    $seq2=~s/\W+//mg;
# Check for N's in FINISHED sequences
    if (($seq2 =~ /N/g)&&($finished)){ 
            $seqerror .= "ACEDBSEQ FINISHED SEQUENCE for $seqname contains N \n";
    }
    if ($seq2 =~ /[^ACGTUMRWSYKVHDBXN]/img) {
      $seqerror .= "ACEDBSEQ for $seqname contains bad characters\n";
      $seq2=~s/[^ACGTUMRWSYKVHDBXN]//img;
    }
    $a_date=$obj->at('Origin.Date_directory[1]');
    $canonical=$obj->Properties(1);
    if ($canonical !~ /Genomic_canonical/) {
      next;
    }
    # Create two Bio::Seq objects and checksum them.
    $bioseq1 = Bio::Seq->new(-seq=>$seq1,-ffmt=>'Fasta',-type=>'Dna',);
    $bioseq2 = Bio::Seq->new(-seq=>$seq2,-ffmt=>'Fasta',-type=>'Dna',);
    $chk1=$bioseq1->GCG_checksum;
    $chk2=$bioseq2->GCG_checksum;
    # Compare date and checksum
    if ($d_date!=$a_date) {
     $seqerror .= "DATE mismatch in $seqname; dir $d_date acedb $a_date\n";
    }
    if ($chk1 != $chk2) {
       $seqerror .= "SEQUENCE mismatch in $seqname; dir $chk1 acedb $chk2\n";
    }
  }
  alarm(0);
};				# end eval
close(CLONEFILE);

@time = localtime();
($MINS,$HOURS,$DAY,$MONTH,$YEAR)=(localtime)[1,2,3,4,5];
if ($time[1]=~/^\d{1,1}$/) {
  $time[1]="0"."$time[1]";
}
$REALMONTH=$MONTH+1;
$REALYEAR=$YEAR+1900;
$TODAY = "$DAY $REALMONTH $REALYEAR at $HOURS:$MINS";
$seqerror=$seqerror."\n\nSeqCompare.pl run ENDED on $TODAY\n\n";
open (OUTFILE,">$outfile");
print OUTFILE $seqerror;
close OUTFILE;
`/usr/bin/mailx -s seqcompare_report ag3\@sanger.ac.uk < $outfile`;
sleep 5;
unlink $outfile;
die();

# Timeout subroutine for avoiding dandling processes
sub timed_out {
 $seqerror=$seqerror."*** Timeout - the process took more than 5000 seconds to complete\n";
}
