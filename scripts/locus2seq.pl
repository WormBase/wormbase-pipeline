#!/usr/local/bin/perl5.6.0 -w
use strict;
use Wormbase;
use Ace;
use utf8;

my $maintainer = "All";
my $rundate    = `date +%y%m%d`; chomp $rundate;
#my $log = "/wormsrv2/logs/$0.$rundate.$$";

my $log = "log.txt";
open(LOG,">$log")|| die "cant open $log";
open (CAMOUT,">../../CAM_locus_seq.ace") || die "cant open CAMOUT";
open (STLOUT,">../../STL_locus_seq.ace") || die "cant open STLOUT";
open (ALLOUT,">../../ALL_locus_seq.ace") || die "cant open ALLOUT";
#print LOG "$0\n";
#print LOG "\n";
#print LOG "=============================================\n";
#print LOG "\n";



#get locus with confirmed CGC names and the corresponding seq
#this uses a table_maker query exported from xace
#my $table = "/wormsrv1/geneace/wquery/locus_seq.def";
#print system(stat $table);
my $command1=<<EOF;
Table-maker -p "/wormsrv1/geneace/wquery/locus_seq.def"
quit
EOF

my %seq_locus;
my $count;$count = 0;
my @entry;
my $seq;
my $locus;
open (GENEACE, "echo '$command1' | tace /wormsrv1/geneace | ");
while (<GENEACE>)
  { 
    @entry = split(/\s+/,$_);
    $locus = $entry[0];
    $seq = $entry[1];

    #this statement is to take in to account the acedb> prompt that is included in the GENACE data
    if (scalar(@entry) > 2)
      {
	$locus = $entry[1];
	$seq = $entry[2];
      }
    #######################################

    if ($locus =~ m/(\w{3}\-\d+\.*\d*)/) #validate cgc naming convention (a v. few genes have ***-*.* eg hmg-1.2
      {
        $locus = $1;#this strips the "'s 
	#if ($count > 10 ){last;}
	if ($seq =~ m/([QWERTYUIOPLKJHGFDSAZXCVBNM0123657894]{2,}\.\w+)/)  #validate sequence name eg XNXNX.XN
	  {
	    $seq = $1;
	    $seq_locus{$seq} .= "$locus ";
	    $count++;
	  }
	else{print LOG "invalid or no sequence found for $locus (found $seq)\n";
	   }
      }
    elsif(scalar(@entry) == 2)#entry no test is to exclude AceDB startup text
      {
	print LOG "$locus is incorrectly marked as cgc approved in genace\n";
      }
  }
close GENEACE;
my $sequence;
#print LOG "#############################\n outputing seq_locus hash\n\n";
#foreach $sequence(keys %seq_locus)
#  {
#    print LOG "$sequence\t$seq_locus{$sequence}\n";
#  }
#print LOG "#############################\n\n\n";
#now find out who did the sequencing

my $autoace = Ace->connect('/wormsrv2/autoace');
my $retrved_seq;
my @lab;
my $CAMcount = 0;
my $STLcount = 0;
my $ALLcount = 0;
my $PROBcount = 0;
my @loci;
foreach $sequence(keys %seq_locus)
  {
    $retrved_seq = $autoace->fetch(Sequence => "$sequence");
    if (defined($retrved_seq))
      {
	@lab = $retrved_seq->at('Origin.From_Laboratory');
	#extract any cases where a sequence contains two loci
	@loci = split(/\s+/,"$seq_locus{$sequence}");
	print "loci is \t@loci\n";
	
	foreach $locus (@loci)
	  {
	    print ALLOUT "Locus: \"$locus\"\nSequence\t\"$sequence\"\n\n";
	    $ALLcount++;
	    
	    print "$sequence\n";
	    if($lab[0] eq "HX")
	      {
		print CAMOUT "Locus: \"$locus\"\nSequence\t\"$sequence\"\n\n";
		$CAMcount++;
	      }
	    elsif($lab[0] eq "RW")
	      {
		print STLOUT "Locus: \"$locus\"\nSequence\t\"$sequence\"\n\n";
		$STLcount++;
	      }
	    else
	      {
		print LOG "locus $locus\t$sequence has undetermined <From_Laboratory> tag - $lab[0]\t \n";
		$PROBcount++;
	      }
	  }
      }
    else
      {
	print LOG "$sequence not found in autoace - locus is $seq_locus{$sequence}\n";
	$PROBcount++;
      }
  }
my $sum = $CAMcount+$STLcount+$PROBcount;
print LOG "found $CAMcount loci on Hinxton sequences.\n
found $STLcount loci on StLouis sequences.\n
found $ALLcount total.\n
$PROBcount have problems\n
ALL should equal sum of others ie $sum and no put into hash - $count \n"; 


close $autoace;
close CAMOUT;
close STLOUT;
close ALLOUT;
close LOG;
