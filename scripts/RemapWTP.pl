
# This was written assuming the file format that the Miegs sent us early June 2002 
# WTP  clone  startpos endpos
# this particular data was not included as it was found to contain an assortment of 
# naming conventions and WTP positionings [email sent to Miegs and WB_dev].

#[020626 ar2].  Script now filters and includes only those WTP's that are chromosome Letter 1-3digit number eg 1A111

#!/usr/local/bin/perl5.6.0 -w
use strict;
use Wormbase;
use Ace;

`rm -f ../../cgc_genes_in_WTP.txt`;

my $maintainer = "All";
my $rundate    = `date +%y%m%d`; chomp $rundate;
my $log = "/wormsrv2/logs/$0.$rundate.$$";

#open(LOG,">$log")|| die "cant open $log";
#print LOG "$0\n";
#print LOG "\n";
#print LOG "=============================================\n";
#print LOG "\n";

#build hash of clones to genomic location
my $query = "CHROMOSOME*";
my $command=<<EOF;
find Sequence $query
follow Subsequence
show -a
quit
EOF

my $db = "/wormsrv2/current_DB";
open (DB, "echo '$command' | tace $db | ");

my @clone_info;
my %clone_pos;
my $superlink;
my @splitline;
my %clone_superlink;
while (<DB>)
  {
    if($_ =~ m/Sequence/)
      {
	@splitline = split(/\s+/,$_);
	$superlink = $splitline[2];
	print "$superlink\n";
      }
    if ($_ =~ m/Subsequence/)
      {	
	@clone_info = split(/\s+/,$_);

	chop $clone_info[1];
	$clone_info[1] = substr($clone_info[1],1);

	$clone_pos{ $clone_info[1] } = $clone_info[2];
	$clone_superlink{ $clone_info[1] } = $superlink;
	#print "$clone_superlink{ $clone_info[1] }\t $clone_info[1]\t $clone_info[2]\n";
      }
  }

close DB;

#my $clonekey;
#sub byseq { $clone_pos{$a} <=> $clone_pos{$b} }
#foreach $clonekey (sort byseq keys %clone_pos)
#{  print "$clonekey - $clone_pos{$clonekey}\n";}

#
my $input = "../../wtpcoords.txt";
my $output = "../../wtpace.ace";
my $crap_data = "../../crapdata.txt";

open (WTPFILE,"<$input") || die "cant open source file";
open (OUT, ">$output") || die "cant write output";
open (CRAP,">$crap_data") || die "cant write crap_data"; 

my $good_count = 0;
my $bad_count = 0;
my @wtp;
my $abs_start;
my $abs_end;
my($wtp, $clone, $rel_start, $rel_end);
my $good;
my $error;
my $cgc_count;
while(<WTPFILE>)
  {
    $good = 0;
    @wtp = split(/\s/,$_);
    if ( (scalar(@wtp) != 4) )
      {
	print CRAP "$_\n";
      }
    
    ($wtp, $clone, $rel_start, $rel_end) = @wtp;
    #print "$wtp, $clone, $rel_start, $rel_end\n";

    chop($wtp);
    chop($clone);
    $wtp = substr($wtp,1);
    $clone = substr($clone,1);

    #extact any CGC names in the file
    #$cgc_count += GetCGCnames($wtp);

    my $teststring .= $wtp.$clone.$rel_start.$rel_end;
    unless( $teststring =~ m/_/)
      {
	if($wtp =~ m/(^\d{1}|X){1}[A-Z]\d+$/)# || ($wtp =~m/^\"[[:lower:]]{3}\-\d+/) )
	  {
	    if($clone =~ m/[A-Z]\d+[A-Z]/)
	      {
		if(( $rel_start =~ m /\d+/) && ( $rel_end =~ m /\d+/))
		  {
		    $abs_start = $rel_start + $clone_pos{$clone};
		    $abs_end = $rel_end + $clone_pos{$clone};
		    
		    print OUT "Sequence : $clone_superlink{$clone}\n";
		    print OUT "WTP \"$wtp\" \"$abs_start\" \"$abs_end\"\n\n";
		    
		    print OUT "WTP : \"$wtp\"\n";
		    print OUT "Method \"WTP\"\n";
		    print OUT "Contact_information \"For more info on transcripts or clones for this gene, please email mieg\@ncbi.nlm.nih.gov.  To ask for clones, please email ykohara\@lab.nig.ac.jp .\"\n\n\n";

		    $good_count++;
		    $good = 1;
		   # if ($good_count > 20){last;}
		  }
		else
		  {
		    print CRAP "position $rel_start or $rel_end no good\n";
		  }
	      }
	    else
	      {
		print CRAP "clone $clone no good\n";
	      }
	  }
	else
	  {
	    print CRAP "wtp $wtp no good\n";
	  }
      }
    else
      {
	print CRAP "got an underscore @wtp\n";
      }
    if( $good == 0) {#print "@wtp $error\n"; 
      $bad_count++;}
  }
    
print "$good_count good ones \t:$bad_count crap ones\n";
close WTPFILE;
close CRAP;
close OUT;

print "\n\nthere are $cgc_count cgc style named loci\n";


sub GetCGCnames
  {
    open (CGC,">>../../cgc_genes_in_WTP.txt");
    my $cgcfile = *CGC;
    my $wtp = shift;
    my $count = 0;
   # if ($wtp =~m/([[:lower:]]{3}\-\d+)/g)
    while ($wtp =~m/([[:lower:]]{3}\-\d+)/g)
      {
	print CGC "\"Locus\" \"$1\"\n";
	$count++;
      }
    close CGC;
    return $count;
  }
