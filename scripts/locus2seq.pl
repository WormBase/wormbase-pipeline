#!/usr/local/bin/perl5.6.0 -w
#
# locus2seq.pl
#
# written by Anthony Rogers (ar2@sanger.ac.uk)
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2002-07-03 14:04:45 $


use strict;
use Wormbase;
use Ace;
use Getopt::Std;
 
##############################
# command-line options       #
##############################
our $opt_d = "";      # Help/Usage page
our $opt_h = "";      # wormpep release number
getopts ('dh');
&usage if ($opt_h);

# choose geneace database
my $geneace_dir;

if ($opt_d){
  $geneace_dir = $opt_d;
}
else{
  $geneace_dir = "/wormsrv2/geneace";
}
print "\nUsing $geneace_dir as target geneace database\n";

my $maintainer = "All";
my $rundate    = `date +%y%m%d`; chomp $rundate;

my $log = "/wormsrv2/logs/locus2seq.log.$rundate";
open(LOG,">$log")|| die "cant open $log";
print LOG "$0\n";
print LOG "Date: $rundate\n";
print LOG "Databases compared: $geneace_dir /wormsrv2/current_DB\n";
print LOG "=====================================================\n";



#get locus with confirmed CGC names and the corresponding seq
#this uses a table_maker query exported from xace
my $command1=<<EOF;
Table-maker -p "/wormsrv2/geneace/wquery/locus_seq.def"
quit
EOF

my %seq_locus;
my $count;$count = 0;
my @entry;
my $seq;
my $locus;


###############################################
# Grab locus->sequence connections from geneace
################################################

open (GENEACE, "echo '$command1' | tace $geneace_dir | ") || die "Couldn't open pipe to $geneace_dir\n";
while (<GENEACE>)
  { 
    @entry = split(/\s+/,$_);
    $locus = $entry[0];
    $seq = $entry[1];
    print "$entry[0]\t$entry[1]\n";

    #this statement is to take in to account the acedb> prompt that is included in the GENEACE data
    if (scalar(@entry) > 2){
      $locus = $entry[1];
      $seq = $entry[2];
    }
    if ((defined($locus))&&($locus =~ m/(\w{3}\-\d+\.*\d*)/)) #validate cgc naming convention (a v. few genes have ***-*.* eg hmg-1.2
      {
        $locus = $1;#this strips the "'s 
	#if ($count > 1 ){last;}
	if ($seq =~ m/([QWERTYUIOPLKJHGFDSAZXCVBNM0123657894]{2,}\.\w+)/)  #validate sequence name eg XNXNX.XN
	  {
	    $seq = $1;
	    $seq_locus{$seq} .= "$locus ";
	    $count++;
	  }
	else{print LOG  "$locus - invalid or no sequence found  (found $seq)\n";
	   }
      }
    elsif(scalar(@entry) == 2)#entry no test is to exclude AceDB startup text
      {
	print LOG "$locus is incorrectly marked as cgc approved in geenace\n";
      }
  }
close(GENEACE);



##########################################
# Compare data to /wormsrv2/current_DB
##########################################

my $autoace_acefiles_dir = "/wormsrv2/autoace/acefiles";
open (CAMOUT,">$autoace_acefiles_dir/CAM_locus_seq.ace") || die "cant open CAMOUT";
open (STLOUT,">$autoace_acefiles_dir/STL_locus_seq.ace") || die "cant open STLOUT";
open (ALLOUT,">$autoace_acefiles_dir/ALL_locus_seq.ace") || die "cant open ALLOUT";



my $sequence;
my $database = "/wormsrv2/current_DB";
my $autoace = Ace->connect($database) || die "cant open $database\n";
my $retrved_seq;
my @lab;
my $CAMcount = 0;
my $STLcount = 0;
my $ALLcount = 0;
my $PROBcount = 0;
my @loci;
my @searched;
my $remark;
foreach $sequence(keys %seq_locus)
  {
    $retrved_seq = $autoace->fetch(Sequence => "$sequence");
    if (defined($retrved_seq))
      {
	@lab = $retrved_seq->at('Origin.From_Laboratory');
	#extract any cases where a sequence contains two loci
	@loci = split(/\s+/,"$seq_locus{$sequence}");
	
	foreach $locus (@loci)
	  {
	    if(defined($lab[0]))
	      {	    
		if($lab[0] eq "HX")
		  {
		    print CAMOUT "Locus: \"$locus\"\nSequence\t\"$sequence\"\n\n";
		    $CAMcount++;
		    print ALLOUT "Locus: \"$locus\"\nSequence\t\"$sequence\"\n\n";
		    $ALLcount++;
		  }
		elsif($lab[0] eq "RW")
		  {
		    print STLOUT "Locus: \"$locus\"\nSequence\t\"$sequence\"\n\n";
		    $STLcount++;	 
		    print ALLOUT "Locus: \"$locus\"\nSequence\t\"$sequence\"\n\n";
		    $ALLcount++;
		  }
		else
		  {
		    print  LOG "\nlocus $locus  -  $sequence has unkown <From_Laboratory> tag - $lab[0]  - not included in outputs \n";
		  }
	      }	
	    else
	      {
		print LOG "$locus -  $sequence has no lab tag\n";
		$PROBcount++;
		FindSequenceInfo($sequence,$locus);
	      }
	  }
      }
    else
      {
	print LOG "\n$sequence not found in $database  - locus is $seq_locus{$sequence}\n";
	$PROBcount++;
	FindSequenceInfo($sequence,$locus);
      }
  }
my $sum = $CAMcount+$STLcount+$PROBcount;
print LOG "found $CAMcount loci on Hinxton sequences.\n
found $STLcount loci on StLouis sequences.\n
found $ALLcount total.\n
$PROBcount have problems\n
ALL should equal sum of others ie $sum and no put into hash - $count \n
\n\n
wrote output ACE files to $autoace_acefiles_dir"; 



my $interested = "Jean TM, John Spieth, Darin Blasair\n";
print LOG "$interested would like to be informed of this update\n\n";

$autoace->close;
close CAMOUT;
close STLOUT;
close ALLOUT;
close LOG;
#$maintainer = "ar2\@sanger.ac.uk";
&mail_maintainer($0,$maintainer,$log);

#copy the ace files to the FTP site

`gzip -f /$autoace_acefiles_dir/STL_locus_seq.ace` && print LOG "gzip failed on STL";
`gzip -f /$autoace_acefiles_dir/CAM_locus_seq.ace` && print LOG "gzip failed on CAM";
`gzip -f /$autoace_acefiles_dir/ALL_locus_seq.ace` && print LOG "gzip failed on ALL";

`cp /$autoace_acefiles_dir/CAM_locus_seq.ace.gz /nfs/privateftp/ftp-wormbase/pub/data/updated_locus2seq/`;
`mv /$autoace_acefiles_dir/STL_locus_seq.ace.gz /nfs/privateftp/ftp-wormbase/pub/data/updated_locus2seq/`;
`mv /$autoace_acefiles_dir/ALL_locus_seq.ace.gz /nfs/privateftp/ftp-wormbase/pub/data/updated_locus2seq/`;


#inform any interested parties
my $notify = $maintainer;
    open (OUTLOG,  "|/usr/bin/mailx -s \"new genace updates\" $notify ");

        print OUTLOG "updated info linking loci to sequences is available from\n
ftp-wormbase\/pub\/data\/updated_locus2seq\/\n
in the 3 files\n
CAM_locus_seq.ace\t loci in Hinxton sequence.\n
STL_locus_seq.ace\t loci in St Louis sequence.\n
ALL_locus_seq.ace\t all loci.\n
\n
These are loci that have approved cgc names.";

    close OUTLOG;



sub FindSequenceInfo #($sequence - genomic seq and $locus )
  {
    my ($seq,$locus,$database) = @_;

    #allows a database to be passed but will default to current_DB
    unless(defined($database)){
      $database = "/wormsrv2/current_DB";
    }
    my $log = "/wormsrv2/logs/locus2seq.log";
    open(LOG,">>$log")|| die "cant open $log";
    my $test_seq = $seq;
    my $autoace = Ace->connect($database) || die "cant connect to $database\n";
    my $solved = 0;
    my @lab;
    my $foundlab;
    my $seq_used;

    #print LOG "examining $seq\n";
    
    if ($test_seq =~ m/(\w+\.)(\d+)$/ )#if the sequence name ends with ".number" (.1   .123)
      {
	my $pre_catch = $1;
	my $catch = $2;


	#print " . . . . . . ends with digit\n";

	#catch where sequence now has isoforms
	$test_seq = $seq."a";
	my $result = TestSeq($test_seq);
	if(defined($result))
	  {
	    print LOG "$locus - $seq does not exist but has isoforms\n\n";
	    $solved = 1;
	  }



	#catch merged sequence eg F34C23.10 has been amalgamated with F34C23.11 or F34C23.9
	else
	  {
	    $catch++;
	    $test_seq = $pre_catch.$catch;
	    #print LOG "testing with seq++ :$seq -> $test_seq\n";
	    if(defined(TestSeq($test_seq)))
	      {
		print LOG "$locus - $seq may have been merged to $test_seq \n\n";
		$solved = 1;
	      }	  
	    else
	      {
		$catch -= 2;
		$test_seq = $pre_catch.$catch;
		print LOG "testing with seq-- :$seq -> $test_seq\n";
		if(defined(TestSeq($test_seq)))
		  {
		    print LOG "$locus - $seq may have been merged to $test_seq \n\n";
		    $solved = 1;
		  }


	#try with just one digit eg  F24G4.23 try F24G4.2
		else
		  {
		    if( length($catch) > 2 )
		      {
			$test_seq = $pre_catch.(substr($catch,0,1));
			if(defined(TestSeq($test_seq)))
			  {
			    print LOG "$locus - $seq not valid - but found $test_seq \n\n";
			    $solved = 1;
			  }
		      }
		  }
	      }
	  }
      }


    ####################################
    #the sequence name ends with a letter

    else   
      {

	#print "$seq\n";
	#print "ends with letter\n";	
	if ($seq =~ m/(\w+\.\d+)([a-z])/)
	{
	  my $letter = $2;
	  my $main_part = $1;

	  #print "\nletter - $letter\tmain - $main_part\n";

	  #increment the last letter eg  F32F7.1a merged -> F32F7.1b
	  my $letter_inc = $letter++;
	  $test_seq = $main_part.$letter_inc;
	 # print "trying with $test_seq . . \n";
	  if(defined(TestSeq($test_seq)))
	    {
	      print LOG "$locus - $seq merged to $test_seq \n\n";
	      $solved = 1;
	    }

	  #leave off the letter and try
	         #if this leaves a bare . eg F32F7. it doesn't matter - this is handled elsewhere
	  else
	    {
	      if(defined(TestSeq($main_part)))
		{
		  print LOG "$locus - $seq isoform not found try $main_part \n\n";
		  $solved = 1;
		}
	    }
	}
    

    #if all else fails just use the root name  ie lop off anything after "."
    if($solved != 1)
      {
	if ($seq =~ m/(\w+)\.\w+/)
	  {
	    $test_seq = $1;
	    if(defined(TestSeq($test_seq)))
	      {
		print LOG "$locus - $seq is not valid seq but $test_seq is\n\n";
	      }
	    else
	      {	      
		print LOG "$locus - $seq -  CANT FIND THIS AT ALL\n\n\n";
	      }
	  }
	else
	  {	      
	    print LOG "$locus - $seq -  CANT FIND THIS AT ALL\n\n\n";
	  }
      }
  }



#retun the last character in a string
sub LastChar 
  {
    my $myString = $_[0];
    return chop($myString);
  }

sub TestSeq # recieves a sequence | returns LabCode if it exists
  {
    
    my $autoace = Ace->connect('/wormsrv2/current_DB');
    my $seq = $_[0];
    my $seq_obj = $autoace->fetch(Sequence => "$seq");
    my $labtag;
    if (defined($seq_obj))
      {
	#get the FROM LAB TAG
	my @lab = $seq_obj->at('Origin.From_Laboratory');
	$labtag = $lab[0];
	unless (defined($labtag))
	  {
	    undef($seq_obj);
	  }
      }
    $autoace->close;
    return $labtag;
  }
}

sub usage {
    system ('perldoc',$0);
    exit;       
}





__END__

=pod

=head2 NAME - locus2seq.pl

=head1 USAGE

=over 4

=item locus2seq.pl  [-options]

=back

This script makes a list of current locus->sequence connections which are valid and
makes a dump of these to the FTP site, making separate files for just St. Louis and 
Sanger sequences and also a combined file.  This information will help keep stlace and
camace synchronised with changes in geneace.

This script also checks for bogus sequences arising from geneace and reports these via
an email message.  This script usually runs at the end of the build process.

locus2seq.pl MANDATORY arguments:

=over 4

=item none

=back

locus2seq.pl  OPTIONAL arguments:

=over 4

=item -d, database

By default this script will compare /wormsrv2/current_DB to /wormsrv2/geneace.  The
-d flag allows you to compare against another copy of geneace (i.e. /wormsrv1/geneace)

=item -h, Help

=back

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
