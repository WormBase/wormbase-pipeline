#!/usr/local/bin/perl5.6.0 -w
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
# Last updated on: $Date: 2002-08-08 14:18:29 $        # quickly see when script was last changed and by whom

use strict;
use Wormbase;
use Ace;

my $maintainer = "All";
my $rundate    = `date +%y%m%d`; chomp $rundate;
my $wmpep_ver = &get_wormbase_version() -1;# for testing during builds
my $wormpepdir = "/wormsrv2/WORMPEP/wormpep$wmpep_ver";
my $log = "/wormsrv2/logs/GetSwissIDandInterpro.WB$wmpep_ver.$rundate";#error log (email to All)
my $temp_acefile = "$wormpepdir/SwissprotIDs.ace";

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

#create temp ace file to output to open and close to ensure new file (as output of data appends)
#open (OLD_IP_OUTPUT,">$temp_clearOldInterPro") || die "cant open $temp_clearOldInterPro";
#close OLD_IP_OUTPUT;
#open (OLD_IP_OUTPUT,">>$temp_clearOldInterPro" || die "cant open $temp_clearOldInterPro");


my @protein;
my @accession;
my @proteinID;
my %idextract;
my %wormpep_acc;
my $count;
$count = 0;

#open file with linking data
open (INPUT, "$wormpepdir/wormpep.table$wmpep_ver")|| print "$0 cant find file wormpep.table$wmpep_ver";

my $accn_holder;
my %noSWALL;   #CEXXXXX => AAMXXXXX.X
my %noswall_acc;

while (<INPUT>)
  {
    chomp;
    if ($_ =~ m /(CE\d+)/)#find wormpep id
      {
	$protein[$count] = $1;
	if($_ =~ m/:([OPQ]\d\w{4})\s/ ) #TrEMBL accessions
	  {
	    $accession[$count] = $1;
	    #build hash of wormpep to accession no.
	    $wormpep_acc{$protein[$count]} = $accession[$count];
	    $count++;
	  }
	else
	  {
	    
	    if( $_ =~ m/(((AA\p{IsUpper})|(CAD))\d{5})\.\d*/ )
	      {
		#put the AAA accs into hash to be investigated
		$noSWALL{$1} = $protein[$count];
		$count++;
	      }
	    else{
	      print LOG "cant find anything useful in $_\n";
	    }
	  }
#	if($count == 20)
#	  {
#	    #DEBUG########################################
#	    last;#only included for testing on small sample sets
#	  }
      }
    else
      {
	print LOG "no protein in $_\n";
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

close LOG;
#### use Wormbase.pl to mail Log ###########
my $name = "Update Swiss ID's and Interpro motifs";
$maintainer = "ar2\@sanger.ac.uk";
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
	if( $databaseID =~ /CAEEL/){
	  $Database = "SwissProt";
	}
	else {
	  if ($databaseID =~ m/[OPQ]\d\w{4}/ ) {
	    $Database = "TrEMBL";
	  }
	  else {
	    if( $databaseID =~ m/((AA\p{IsUpper})|(CAD))\d{5}/ ) {
	      $Database = "TrEMBLnew";
	    }
	    else {print $errorLog $peptide,"\t",$full_accn,"\tfrom unknown database\n";}
	  }
	}
	
	if(defined($databaseID))
	  {
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
    foreach $aaaID (keys %noSWALL)
      {
	push(@AAAs,$aaaID);
      }
    print "@AAAs";
    my $submitString = "pfetch "." @AAAs";
    my %idextract;

    #create the fasta record in /wormsrv2/tmp dir (and remove it after use)
    my $fasta_output = "/wormsrv2/tmp/fasta_output";
    open (FASTA,">$fasta_output")||die " cant open fasta_output";# > clears before each use 
    #this submits request and put results in file
    print FASTA `$submitString`;
    close FASTA;
    
    my $fasta_output_clean = "$fasta_output"."_clean";
    `grep -F '>' $fasta_output > $fasta_output_clean`;


    open (FH, "<$fasta_output_clean");
      while(<FH>)
	{
	  # while global match as some proteins have multiple ids
	  while( $_ =~ m/(((AA\p{IsUpper})|(CAD))\d{5})\.\d*/g )   #matched AAM23424.2 or CAD23525.1 stylee
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
