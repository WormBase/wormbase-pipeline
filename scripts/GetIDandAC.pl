#!/usr/local/bin/perl5.6.0 -w
#
# GetSwissIdandInterpro.pl
# 
# by Anthony Rogers
#
# This script parses a list of genes and Protein ids, pfetchs them  and extracts SwissProt /TrEMBL accn no.s. 
# Output is a .ace file as such ;
#
# Sequence : "F22H6.2"
# Database "SwissProt" "FLO1_CAEEL" "Q17766"
# Protein_id "AAM245522" "2"         id and version
#
# pfetch is done in batches of 2000, any greater and nothing comes back!
#
# Last updated by: $Author: ar2 $                  
# Last updated on: $Date: 2002-09-06 12:25:35 $    

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;

my $maintainer   = "All";
my $rundate      = `date +%y%m%d`; chomp $rundate;
my $wmpep_ver    = &get_wormbase_version();# for testing during builds
my $log          = "/wormsrv2/logs/GetIDandAC.$rundate";#error log (email to All)
my $wormpepdir   = "/wormsrv2/WORMPEP/wormpep$wmpep_ver";
my $acefile = "/wormsrv2/autoace/wormpep_ace/WormpepACandIDs.ace";

open(LOG,">$log")|| die "cant open $log";
my $errorLog = *LOG;
print LOG "#$0\n";
print LOG "\n";
print LOG "=============================================\n";
print LOG "\n";

#create temp ace file to output to |open and close to ensure new file (as output of data appends)
open (ACE_OUTPUT,">$acefile") || die "cant open $acefile";

#open file with linking data
open (INPUT, "/nfs/disk100/wormpub/analysis/SWALL/output_autoace")|| print "$0 cant find file";

my %protid_data;    #$protid_data{$protein_id}= @[gene, ver_no,  AC, ID, parent_clone ]
my @AAAprotein_ids;
my $protein_id;
my $gene;
my $parent_clone;

while (<INPUT>)
  {
    chomp;
    ($parent_clone) = (/^EMBL \[(\S+)\|\S+\]/);                # modified to unprocessed format
    ($gene,$protein_id) = (/gene:(\S+)\s+protein_id:(\S+)/);   # (no need for intermediary script)
    $protein_id =~ m/(\w+)\.(\d)/;
    $protein_id = $1;
    my $version = $2;
    $protid_data{$protein_id}[0] =$gene;
    $protid_data{$protein_id}[1] = $version;
    $protid_data{$protein_id}[4] = $parent_clone;
    push(@AAAprotein_ids, $protein_id);
  }

#hash and accession array are built in one go 
#now submit unit slices to outputToAce
my $chunk_size = 2000;
my $last_ind = $#AAAprotein_ids;
my $upper = $chunk_size -1;
my $lower = 0;
my @array_chunk;
my $count;

while ( $lower <= $last_ind ) 
  {
     @array_chunk = @AAAprotein_ids[$lower .. $upper];
     
     my $submitString = "pfetch "." @array_chunk";

     #create the fasta record in /wormsrv2/tmp dir (and remove it after use)
     my $fasta_output = "/wormsrv2/tmp/fasta_output";
     open (FASTA,">$fasta_output")||die " cant open fasta_output";# > clears before each use 
     #this submits request and put results in file
     print FASTA `$submitString`;
     close FASTA;
     
     my $fasta_output_clean = "$fasta_output"."_clean";
    `grep -F '>' $fasta_output >> $fasta_output_clean`;
     
     my $aaaID;
     open (FH, "<$fasta_output_clean");
     while(<FH>)
       {
	 # while global match as some proteins have multiple ids
	  while( $_ =~ m/(((AA\p{IsUpper})|(CA\p{IsUpper}))\d{5})\.\d*/g )   #matched AAM23424.2 or CAD23525.1 stylee
	    {
	      $aaaID = $1;
	      if ( defined($protid_data{$aaaID}) )  #if the $aaaID is one we found earlier
		{
		  $count++;
		  # this is an AAMXXXXX.X style id that we know about
		  if( $_ =~ m/\w+_CAEEL/) {
		    $protid_data{$aaaID}[3] = $&;
		    if( $_ =~ m/([OPQ]\d\w{4})\s/ ) {
		      $protid_data{$aaaID}[2] = $1;
		    }
		    else {
		      print LOG "$aaaID has swissprot but no AC\n";
		    }
		  }
		  elsif( $_ =~ m/([OPQ]\d\w{4})\s/ ) {
		    $protid_data{$aaaID}[2] = $1;
		    $protid_data{$aaaID}[3] = $1;
		  }
		  else { 
		    $protid_data{$aaaID}[2] = $aaaID;
		    $protid_data{$aaaID}[3] = $aaaID;  
		  }
		  last;
		}
	    }
	  unless( defined( $protid_data{$aaaID}[3]) )
	    {
	      print "badone";
	    }
	}
     close FH;

     $lower += $chunk_size;
     $upper += $chunk_size;
  #   if( $lower >= $last_ind ) {
#       last;
#     }
     if( $upper > $last_ind ) {
       $upper = $last_ind;
     }
  }

#write the ace file
my $printcount;
foreach my $key (keys %protid_data)
  {
    my $database;
    if( defined($protid_data{$key}[3]) )
      {
	$printcount++;
	if( $protid_data{$key}[3] =~ m/_CAEEL/ ) {
	  $database = "SwissProt";
	}
	elsif( $protid_data{$key}[3] =~ m/[OPQ]\d\w{4}/ ) {
	  $database = "TREMBL";
	}
	else {
	  $database = "TREMBLNEW";
	}
	print ACE_OUTPUT "Sequence : \"$protid_data{$key}[0]\"\n";
	print ACE_OUTPUT "Database \"$database\" \"$protid_data{$key}[3]\" \"$protid_data{$key}[2]\"\n";
	print ACE_OUTPUT "Protein_id \"$protid_data{$key}[4]\" \"$key\" \"$protid_data{$key}[1]\"\n";
	print ACE_OUTPUT "\n";
      }
  }
close ACE_OUTPUT;
close LOG;

print "processed $count defined IDs\n";
print "printed $printcount objects\n\n";
#### use Wormbase.pl to mail Log ###########
my $name = "Update Swiss ID's and Interpro motifs";
#$maintainer = "ar2\@sanger.ac.uk";
&mail_maintainer ($name,$maintainer,$log);
#########################################
