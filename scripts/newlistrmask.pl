#!/usr/local/bin/perl

#  newlistrmask.pl - Worm genomic DNA retrieval;  repeats retrieval and masking; 
#  Blast search through LSF;  postfiltering of search results in "acedb" mode.
#  ver 2.00, Aug 1999, A.Guffanti. Uses AQL

$|=1;
$path="/nfs/disk100/wormpub/acedb/ace4/cam";
$DIR="/nfs/disk65/ag3/FILTERED";
$SOURCEFILE="/nfs/disk100/wormpub/analysis/cosmids/current.versions";
$SOURCEDIR="/nfs/disk100/wormpub/analysis/cosmids/";
$ENV{'ACEDB'}="/nfs/disk100/wormpub/acedb/ace4/cam";
$exec=&tace;
@aqlarray=();
@genseq=();
$outfile="/tmp/newlistrmask.$$";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time);
$date="$mday"."/"."$mon"."/"."$year".":"."$hour".":"."$min";
$message = "newListrMask STARTED at $date\n";

@skip=('B0393','R10E9','F59B2','F43G6','ZC410','C26E1','F53B2','T06E6','F18A11','C41G6','F19B2','D1053','Y87G2A','Y39C12B','Y50E8A');

open (CURRFILE,"$SOURCEFILE");
  while (<CURRFILE>){
  $file = $_;
  chomp $file;
  ($currseq,$date)=split (/\//, $file);
  $currentfile = "$SOURCEDIR"."$file"."/"."$currseq.fasta";
  @currentmatch=();
  @currentmatch = grep {/$currseq/} @skip;
  if (@currentmatch>0) {
    print "Skipping sequence $currseq\n";
    next;
  }
  print "Now examining $currseq\n";

## Get all the names for genomic sequences with DNA
#$query1 = <<END1;
#aql select q from q in class Genome_Sequence where exists_tag q->DNA
#END1
#&GetNames;
#foreach (@genseq) {
#  $currseq=$_;
#  chomp $currseq;
#  if (length $currseq ==0) {
#    next;
#  }

 # Extract DNA sequence
 $dna="";
 open (CURRENTFILE,"<$currentfile");
  while (<CURRENTFILE>) {
    /\>/ && next;
    chomp $_;   
    $dna.=$_;  
  }  
  close (CURRENTFILE);

# Only for retrieval DNA from GENOMIC sequences
#  $query2 =<<END2;
#query find genome_sequence $currseq
#DNA
#END2
#  &GetDna;
  
# Extract inverted and tandem repeats
$query3=<<END3;
aql select m[1], m[2] from s in object ("Sequence","$currseq"), m in s->Feature where m like "tandem" or m like "inverted"
END3
print "INV/REP QUERY: $query3\n";
eval {
 &GetFeatures($query3);
};
print "INV/REP ERROR MESSAGE: $@\n" if $@;

# Extract CeRepeats						   
$query4=<<END4;
aql select m[3], m[4] from s in object("Sequence","$currseq"), m in s->Motif_homol where m like "CeRep*"
END4
print "CeRep QUERY: $query4\n";
eval {
 &GetFeatures($query4);
};
print "CeRep ERROR MESSAGE: $@\n" if $@;

 &RepBuild;

eval {
 &Mask;
};
print "MASKING ERROR MESSAGE: $@\n" if $@;

# Write out masked DNA  from $dna
 &WriteFasta($dna);
  print "MASKED $currseq\n";
 undef ($dna);
}

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time);
$date="$mday"."/"."$mon"."/"."$year".":"."$hour".":"."$min";
$message .= "newListrMask ENDED at $date\n";
open (OUTFILE,">$outfile");
print OUTFILE $message;
close OUTFILE;
`/usr/bin/mailx -s newlistrmask_report ag3\@sanger.ac.uk < $outfile`;
sleep 5;
unlink $outfile;

die();

###################### SUBS ###################

#-----------------------
# Execute tace
#
sub tace {
#-----------------------
   local($prog);
   local($name); 
   $name=`uname -sr`;
    if ($name=~/^SunOS/) {($prog)=</nfs/disk100/acedb/RELEASE.SUPPORTED/bin.SOLARIS_4/tace>;}
    elsif ($name=~/^IRIX/) {($prog)=</nfs/disk100/acedb/RELEASE.SUPPORTED/bin.SGI_4/tace>;}
    elsif ($name=~/^OSF/)  {($prog)=</nfs/disk100/wormpub/ACEDB/bin_ALPHA/tace>;}
    elsif ($name=~/^Linux/) {($prog)=</nfs/disk100/acedb/RELEASE.SUPPORTED/bin.LINUX_4/tace>;}
    else {print STDERR "No known binary for $uname\n";exit;}
    return $prog;
}

#-----------------------
# Get Sequence Names
#
sub GetNames {
#-----------------------
open (QUERY1,"echo '$query1' | $exec |");
while (<QUERY1>) {
   /FATAL ERROR/ && return;
   /0 Active Objects/ && return; 
   /\/\// && next;
   /Please type \? for a list of commands/ && next;
   s/acedb\>\s*//;
   s/\"//mg;
   s/^\s*//;
   s/\s*$//;
   push(@genseq,$_);
}
close QUERY1;
}

#-----------------------
# Retrieval of DNA sequence
#
sub GetDna {
#-----------------------
$dna="";
open (QUERY2,"echo '$query2' | $exec |");
while (<QUERY2>) {
   /FATAL ERROR/ && return;
   /0 Active Objects/ && return; 
   /^\/\// && next;
  />/ && next;
  /Please type \? for a list of commands/ && next;
  s/acedb\>\s*//;
  s/^\s+//mg;
  s/\s+$//mg;  
  chomp;
  $dna.=$_;
}
close QUERY2;
}

#-------------------------
# AQL - get Inverted, Tandem and
# CeRep features
#
sub GetFeatures {
#-------------------------
$query=shift @_;
chomp $query;
open (QUERY,"echo '$query' | $exec |");
while (<QUERY>) {
   /FATAL ERROR/ && return;
   /\/\// && next;
   /Please type \? for a list of commands/ && next;
   s/acedb\>\s*//;
   s/\"//mg;
   s/^\s+//;
   s/\s+$//;
   if (length($_)==0){next;}
   push (@aqlarray,"$_\n");
}
close QUERY;
}

#----------------------
# Construction of associative
# repeat array for masking
#
sub RepBuild {
#----------------------
  foreach (@aqlarray){
   my ($begin,$end) = split (/\s+/,$_);
  $REPS{$begin}=$end;
   }
 @aqlarray=();
}

#----------------------
# Masking of DNA sequence
#
sub Mask {
#----------------------
  my $t, $len, $begin, $end ;
  foreach $begin (sort numerically keys %REPS) {
    $end=$REPS{$begin};
    if ($begin > $end) { $t = $begin ; $begin = $end ; $end = $t ; }
    $len = $end-$begin ;
    if ($len == 0) {
      next;
    }
    substr($dna,$begin-1,$len) = "N" x $len ;
 #print "BEGIN: $begin\n";
 #print "END: $end\n";
 #print "LENGTH: $len\n";
  }
  %REPS={};
}

#----------------------
# Writing of masked DNA 
# in FastA format
#
sub WriteFasta {
#----------------------
  $outdir=$SOURCEDIR."$file"."/blast2";
  if (! -d $outdir) {
    `mkdir $outdir`;
  }
    `chmod -R a+rx $outdir`;
  print "$outdir now OK\n";
  $dnaout="$outdir"."/"."$currseq".".mask";
  my $len, $pos ;
  open (OUTFILE2,">$dnaout"); 
  print OUTFILE2 ">$currseq\n" ;
  $len = length($dna) ;
  for ($pos = 0 ; $pos < $len-60 ; $pos += 60) {
    print OUTFILE2 substr($dna, $pos, 60) . "\n" ;
  }
  if ($pos < $len) { print OUTFILE2 substr ($seq,$pos) . "\n" ; }
  close OUTFILE2;
}
  
#--------------------------------
# Numerical sorting of hashes
#--------------------------------
sub numerically {$a <=>$b;}


######################################
#
#AcePerl version - too memory hungry !
#
######################################

#use Ace;
#$db=Ace->connect(-path=>$path);

#$query1=<<END;
#query find genome_sequence DNA
#END
#my @list = map $_->name(), $db->fetch(-query => $query1 );
#foreach $name (@list) {
#  my $seq = $db->fetch( Sequence => $name );
#  $dna = $seq->asDNA();
#  $dnaout="$DIR/$name.mask";
#  @inverted1=$seq->at ('feature.inverted');
#  @inverted2=$seq->at ('feature.inverted[2]');
#  @tandem1=$seq->at ('feature.tandem');
#  @tandem2=$seq->at ('feature.tandem[2]');
  
### Extract CeReps through AQL query
# $aqlquery = <<END2;
#aql select s->motif_homol[4],s->motif_homol[5] from s in object ("Sequence","$name") where s->motif_homol[1] like "CeRep*"
#END2

#  my $aql=$db->raw_query($aqlquery);
#  $aql=$db->fetch(-query=>$aqlquery);
#  $aql=~s/\/\/.+//mg;
#  $aql=~s/^\s+//mg;
#  $aql=~s/\s+$//mg;
#  @aqlarray=split ' ',$aql;
#   undef $aql;
  
#  # Mask inverted repeats in $dna
#  &RepBuild(\@inverted1,\@inverted2);
#  &Mask;
#  # Mask tandem repeats in $dna
#  &RepBuild(\@tandem1,\@tandem2);
#  &Mask;
#  # Mask Ce repeats in $dna
#  &CeRepBuild;
#  &Mask;
#  # Write out masked DNA  from $dna
#  # Undef all objects
#  &WriteFasta($dnaout);
#  undef $dna;
#  undef $obj;
#}

##----------------------
## Construction of associative 
## repeat array for masking
## Ce Repeats
##
#sub CeRepBuild {
##----------------------
# my $i,$begin,$end; 
# for ($i=0;$i<=$#aqlarray;$i=$i+2){
#    ($begin,$end) =  splice (@aqlarray,0,2);
#    $REPS{$begin}=$end;
#  }
#}
