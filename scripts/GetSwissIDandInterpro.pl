#!/usr/local/bin/perl5.6.1 -w
#
# GetSwissIdandInterpro.pl
# 
# by Anthony Rogers
#
# This script parses the current wormpep.table file and extracts SwissProt /TrEMBL accn no.s. It then pfetches the
# full FASTA record and extracts any Database IDs (eg XXX_CAEEL ) and Interpro motifs.
# Output is a .ace file as such ;
#
# "Protein : "WP:CE05236"
# Database "SwissProt" "FLO1_CAEEL" "Q17766"
# Motif_homol     "INTERPRO:IPR002666"
#
# pfetch is done in batches of 2000, any greater and nothing comes back!
#
# Last updated by: $Author: ar2 $                      # These lines will get filled in by cvs 
# Last updated on: $Date: 2002-10-21 12:48:52 $        # quickly see when script was last changed and by whom

use strict;
use lib "/wormsrv2/scripts/";                  
use Wormbase;
use Ace;


use Getopt::Std;
#######################################
# command-line options                #
#######################################

use vars qw($opt_d);
# $opt_d debug   -  redirect output

getopts ('d');

my $maintainer = "All";
my $rundate    = `date +%y%m%d`; chomp $rundate;
my $log;
my $temp_acefile;

my $wmpep_ver = &get_wormbase_version();#-1 for testing during builds
my $wormpepdir = "/wormsrv2/WORMPEP/wormpep$wmpep_ver";

if( defined($opt_d) )
  {
    $log = glob("~ar2/testlogs/GetSwissIDandInterpro.WB$wmpep_ver.$rundate");
    $temp_acefile = glob("~ar2/testlogs/SwissprotIDs.ace");
    $wormpepdir = "/wormsrv2/WORMPEP/wormpep$wmpep_ver";
    $maintainer = "ar2\@sanger.ac.uk";
  }
else
  {
    $log = "/wormsrv2/logs/GetSwissIDandInterpro.WB$wmpep_ver.$rundate";
    $temp_acefile = "/wormsrv2/autoace/wormpep_ace/SwissprotIDs.ace";
  }

my $ace_output = *ACE_OUTPUT;

open(LOG,">$log")|| die "cant open $log";
my $errorLog = *LOG;
print LOG "#$0\n";
print LOG "\n";
print LOG "Wormpep version  :$wmpep_ver\n\n";
print LOG "=============================================\n";
print LOG "\n";

#create temp ace file to output to |open and close to ensure new file (as output of data appends)
open (ACE_OUTPUT,">$temp_acefile") || die "cant open $temp_acefile";
close ACE_OUTPUT;
open (ACE_OUTPUT,">>$temp_acefile" || die "cant open $temp_acefile");

my @protein;
my @accession;
my @proteinID;
my %idextract;
my %wormpep_acc;

#open file with linking data
open (INPUT, "$wormpepdir/wormpep.table$wmpep_ver")|| print "$0 cant find file wormpep.table$wmpep_ver";

my $accn_holder;
my %noSWALL;   #CEXXXXX => AAMXXXXX.X
my %noswall_acc;

my $linecount = 0;
my $noCE_count = 0;
my $SWTR_count = 0;
my $AAA_count = 0;
my $useless_count = 0;
my $count = 0;
while (<INPUT>)
  {
    chomp;
    $linecount++;
    if ($_ =~ m /(CE\d+)/)#find wormpep id
      {
	$protein[$count] = $1;
	#if($_ =~ m/:([OPQ]\d\w{4})\s/ ) #TrEMBL accessions
	if($_ =~ m/(SW|TR):(\S+)/ ) #TrEMBL accessions
	  {
	    $accession[$count] = $2;
	    #build hash of wormpep to accession no.
	    $wormpep_acc{$protein[$count]} = $accession[$count];
	    $SWTR_count++;
	    $count++;
	  }
	else
	  {
	    
	    if( $_ =~ m/(((AA\p{IsUpper})|(CA\p{IsUpper}))\d{5})\.\d*/ )
	      {
		#put the AAA accs into hash to be investigated
		$noSWALL{$1} = $protein[$count];
		$AAA_count++;
		$count++;
	      }
	    else{
	      print LOG "cant find anything useful in $_\n";$useless_count++;
	    }
	  }
#	if($count == 2000)
#	  {
#	    #DEBUG########################################
#	    last;#only included for testing on small sample sets
#	  }
      }
    else
      {
	print LOG "no protein in $_\n";$noCE_count++;
      }
  }
#try and get AC for those that have AAM or CAD style protein IDs only
&GetNoAccPeps;


#hash and accession array are built in one go 
#now submit unit slices to outputToAce
my $chunk_size = 2000;
my $last_ind = $#accession;
my $upper = $chunk_size -1;
my $lower = 0;
my @array_chunk;


while ( $lower <= $last_ind ) 
  {
     @array_chunk = @accession[$lower .. $upper];
     outputToAce(\%wormpep_acc, \@array_chunk, \$ace_output, \$errorLog);
     $lower += $chunk_size;
     $upper += $chunk_size;
  }

@array_chunk = @accession[$lower .. $last_ind];
##process the remainders
outputToAce(\%wormpep_acc, \@array_chunk, \$ace_output, \$errorLog);


close ACE_OUTPUT;

print LOG "Files available -\n $temp_acefile\n";
print LOG "\nNOTE! - This script no longer parses the outputs to a database\n\n";
print LOG "SubGetSwissId finished at ",`date`,"\n";
print LOG "linecount = $linecount\n
$noCE_count = noCE_count\n
$SWTR_count = SWTR_count\n
$AAA_count = AAA_count\n
$useless_count =useless_count\n";

close LOG;
#### use Wormbase.pl to mail Log ###########
my $name = "Update Swiss ID's and Interpro motifs";
&mail_maintainer ($name,$maintainer,$log);
#########################################

#===========================================================================
sub outputToAce #(\%wormpep_acc, \@accession \$ace_output, \$errorLog)
  {
    my ($wormpep_acc, $chunk, $ace_ouput, $errorLog) = @_;
    #$ace_output is filehandle of file creatred in calling script

    #construct pfetch command string by concatenating accession no's
    my $submitString = "pfetch -F"." @$chunk";#get full Fasta record (includes Interpro motifs)
    my %idextract;

    #create the fasta record in /wormsrv2/tmp dir (and remove it after use)
    my $fasta_output = "/wormsrv2/tmp/fasta_output";

    open (FASTA,">$fasta_output")||die " cant open fasta_output";# > clears before each use 
    #this submits request and put results in file
    print FASTA `$submitString`;
    close FASTA;
    
    #build hash of accession : Swissprot/Trembl ID
    open (FASTA,"<$fasta_output")|| die "cant open fasta record";#read
    my %interpro;
    my $proteinID;
    my $acc;
    print "opening FASTA file  . . . \n";
    while(<FASTA>)
      {
	chomp;
	my $record .= "$_\n";
	if ($_ =~ /^\/\//)
	  {
	    if (defined $proteinID)
	      {
		print "$proteinID went fine\n";
		undef $proteinID;
		undef $acc;
		$record = "";
	      }
	    else{
	      print "gone thru a record with out picking up protein - \n$record\n";
	      die;
	    }
	  }
	#get the id
	if($_ =~ m/^ID\s+(\w+(\.\d)*)/)
	 { $proteinID = $1;}
	
	#get the accession
	if($_ =~ m/^AC\s+(\w+)/)
	  {
	    $acc = $1;
	    #put ID AC pair in to hash
	    $idextract{$acc} = $proteinID
	  }

	#extract Interpro motifs
	while ($_ =~ m/(IPR\d+)/g){
	  print "found IP $1 for protein -$proteinID-end\n";
	  $interpro{$proteinID} .= "$1 ";}
	
		    
      }
    close FASTA;
    `rm -f $fasta_output`;

    foreach my $peptide(keys %wormpep_acc)
      {
	my $full_accn = $wormpep_acc->{$peptide};
	my @splitInterPros;
	my $databaseID = $idextract{$full_accn};
	my $Database;
	#catch the new SwissProt entries
	if(defined($databaseID))
	  {
	    if( $databaseID =~ /CAEEL/){
	      $Database = "SwissProt";
	    }
	    else {
	      if ($databaseID =~ m/[OPQ]\d\w{4}/ ) {
		$Database = "TrEMBL";
	      }
	      else {
		if( $databaseID =~ m/((AA\p{IsUpper})|(CA\p{IsUpper}))\d{5}/ ) {
		  $Database = "TrEMBLnew";
		}
		else {print $errorLog $peptide,"\t",$full_accn,"\tfrom unknown database\n";}
	      }
	    }
	    
	    
	    print $ace_output "Protein : \"WP:$peptide\"\n";
	    print $ace_output "Database \"$Database\" \"$databaseID\" \"$full_accn\"\n";
	    
	    if (defined($interpro{$databaseID}))
	      {
		@splitInterPros = split(/\s/,"$interpro{$databaseID}");#interpros are stored as concatenated sting value in hash
		my $ip;
		while(defined($ip = pop(@splitInterPros)))
		  {
		    print $ace_output "Motif_homol\t\"INTERPRO:$ip\"\n";#Motif_homol      "INTERPRO:IPR000210"
		  }
	      }
	    print $ace_output "\n";#object separator
	  }
      }
  }

sub GetNoAccPeps
  { 
    #construct pfetch command string by concatenating accession no's
    my @AAAs;
    my $aaaID;
    my $pep;
    my $CE;
    my $fasta_output;
    foreach $aaaID (keys %noSWALL)
      {
	push(@AAAs,$aaaID);
      }
    #create the fasta record in /wormsrv2/tmp dir (and remove it after use)
    $fasta_output = "/wormsrv2/tmp/fasta_output";
    open (FASTA,">$fasta_output")||die " cant open $fasta_output";# > clears before each use 

    #PUT THIS IN 2000 CHUNK LOOP
    my $chunky_size = 2000;
    my $AAAlow = 0; 
    my $AAAup = $chunky_size - 1;
    while ( $AAAs[$AAAlow] )
      {
	#print "@AAAs";
	my @AAAchunk = @AAAs[$AAAlow .. $AAAup];
	my $submitString = "pfetch "." @AAAchunk";
	my %idextract;
	
	#this submits request and put results in file
	print FASTA `$submitString`;
	$AAAlow += $chunky_size;
	$AAAup += $chunky_size;	
      }
	
    close FASTA;
    
    my $fasta_output_clean = "$fasta_output"."_clean";
    `grep -F '>' $fasta_output > $fasta_output_clean`;


    open (FH, "<$fasta_output_clean");
      while(<FH>)
	{
	  # while global match as some proteins have multiple ids
	  while( $_ =~ m/(((AA\p{IsUpper})|(CA\p{IsUpper}))\d{5})\.\d*/g )   #matched AAM23424.2 or CAD23525.1 stylee
	    {
	      $aaaID = $1;
	      if ( defined($noSWALL{$aaaID}) )  #if the $aaaID is one we found earlier
		{
		  # this is an AAMXXXXX.X style id that we know about
		  $CE = $noSWALL{$aaaID};
		  if( $_ =~ m/([OPQ]\d\w{4})\s/ ) {
		    $noswall_acc{$CE} = $1;
		  }
		  else {
		    $noswall_acc{$CE} = $aaaID;
		  }
		  last;
		}
	      else {
		print LOG  "Cant find $aaaID in hash\n";
	      }
	    }
	}
    close FH;
    foreach my $CE (keys %noswall_acc)
      {
	push(@accession, $noswall_acc{$CE});
	$wormpep_acc{$CE} = $noswall_acc{$CE};
      }
  }
