# This was written assuming the file format that the Miegs sent us early June 2002 
# WTP  clone  startpos endpos
# this particular data was not included as it was found to contain an assortment of 
# naming conventions and WTP positionings [email sent to Miegs and WB_dev].


#!/usr/local/bin/perl5.6.0 -w
use strict;
use Wormbase;
use Ace;

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

my $db = "/wormsrv2/autoace";
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

open (WTPFILE,"<$input") || die "cant open source file";
open (OUT, ">$output") || die "cant write output";


my $count;
my @wtp;
my $abs_start;
my $abs_end;
my($wtp, $clone, $rel_start, $rel_end);
while(<WTPFILE>)
  {
    ($wtp, $clone, $rel_start, $rel_end) =  split(/\s/,$_);
    # if( undef($wtp) || undef($clone) || undef($rel_start) || undef($rel_end) || undef($clone_pos{$clone}) )
    #if(  defined($clone_pos{$clone}) )
    #   {die "\n\ndied coz of clone $clone\n\n";}
    
    $abs_start = $rel_start + $clone_pos{$clone};
    $abs_end = $rel_end + $clone_pos{$clone};
    
    print OUT "Sequnce: $clone_superlink{$clone}\n";
    print OUT "WTP $wtp\tclone $clone\tPos: $abs_start\t$abs_end\n\n"; 
    $count++;
   # if ($count > 20){last;}
  }
print "$count\n";
close WTPFILE;
close OUT;
