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
# Last updated on: $Date: 2002-08-05 10:44:30 $        # quickly see when script was last changed and by whom

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
my $old_ip_output =*OLD_IP_OUTPUT;

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
while (<INPUT>)
  {
    chomp;
    if ($_ =~ m /(CE\d+)/)#find wormpep id
      {
	$protein[$count] = $1;
	if($_ =~ m/(SW:\w+|TR:\w+)/)#any "words" prefixed with SW: or TR:
	  {
	    $accn_holder = $1;
	    #print "protein $count $protein[$count] accn is $accn_holder\n";

	    $accession[$count] = $accn_holder;
	    #build hash of wormpep to accession no.
	    $wormpep_acc{$protein[$count]} = $accession[$count];
	    $count++;

	  }
	else
	  {
	    print LOG "no TR: SW: accn in $_\n";
	    if( $_ =~ m/((AA\p{IsUpper})|(CAD))\d{5}\.\d*/ )
	      {
		$noSWALL{$&} = $protein[$count];
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

#hash and accession array are built in one go 
#now submit unit slices to outputToAce
my $chunk_size = 2000;
my $last_ind = $#accession;
my $upper = $chunk_size -1;
my $lower = 0;
my @array_chunk;

#try and get AC for those that have AAM or CAD style protein IDs only
&GetNoAccPeps;

#while ( $lower <= $last_ind ) 
#  {
#     @array_chunk = @accession[$lower .. $upper];
#    outputToAce(\%wormpep_acc, \@array_chunk, \$ace_output, \$old_ip_output, \$errorLog);
#    $lower += $chunk_size;
#    $upper += $chunk_size;
#  }

@array_chunk = @accession[$lower .. $last_ind];
##process the remainders
outputToAce(\%wormpep_acc, \@array_chunk, \$ace_output, \$old_ip_output, \$errorLog);







close OLD_IP_OUTPUT;
close ACE_OUTPUT;

print LOG "Files available -\n $temp_acefile\n
$temp_clearOldInterPro\n";
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
    my ($wormpep_acc, $chunk, $ace_ouput, $old_ip_output, $errorLog) = @_;
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
    print FASTA "This is a temp file created by $0 and can be removed -unless script is running!";
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
	#get the SWISSPROT/TrEMBL id
	if($_ =~ m/^ID\s+(\w+)/)
	 { $proteinID = $1;}
	
	#get the SWISSPROT accession
	if($_ =~ m/^AC\s+(\w+)/)
	  {
	    $acc = $1;

	    #put ID AC pair in to hash
	    $idextract{$acc} = $proteinID;
	  }

	#extract Interpro motifs
	while ($_ =~ m/(IPR\d+)/g){
	  print "found IP $1 for protein -$proteinID-end\n";
	  $interpro{$proteinID} .= "$1 ";}
	
	#print "$interpro{$proteinID}\n";
		    
      }
    close FASTA;
    `rm -f $fasta_output`;
    # ! ! ! ! resetting record separator ! ! ! ! #
   # $/ = '\n'; #\n

    #DEBUG#############################
    #print "idextract contents\n";
    #foreach my $key (keys %idextract){
    #  print $key,"\t",$idextract{$key},"\n";}
    ##################################

    #DEBUG########################################
    #print "wormpep_acc hash after building\n";
    #foreach my $worm (keys %wormpep_acc){
    #  print "$worm","\t",$wormpep_acc{$worm},"\n";}
    #############################################

    foreach my $peptide(keys %wormpep_acc)
      {
	my $full_accn = $wormpep_acc->{$peptide};
	my $accn = substr($full_accn,3);
	my $Database = substr($full_accn,0,2);
	if($Database eq "SW"){$Database = "SwissProt";}
	elsif($Database eq "TR"){$Database = "TrEMBL";}
	else{print $errorLog $peptide,"\t",$full_accn,"\tfrom unknown database\n";}

	

	my @splitInterPros;
	my $databaseID = $idextract{$accn};
	#catch the new SwissProt entries
	if( $databaseID =~ /CAEEL/){
	  $Database = "SwissProt";
	}
	if(defined($databaseID))
	  {
	    #remove any old InterPro entries for any proteins with SP / Tr Id's
	    #print $old_ip_output "Protein : \"WP:$peptide\"\n";
	    #print $old_ip_output "-D Motif_homol\t\"INTERPRO:\"\n";
	    #print $old_ip_output "\n";#object separator

	    print $ace_output "Protein : \"WP:$peptide\"\n";
	    print $ace_output "Database \"$Database\" \"$databaseID\" \"$accn\"\n";
	    
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
    my $submitString = "pfetch "." @AAAs";#get full Fasta record (includes Interpro motifs)
    my %idextract;

    #create the fasta record in /wormsrv2/tmp dir (and remove it after use)
    my $fasta_output = "/wormsrv2/tmp/fasta_output";
    open (FASTA,">$fasta_output")||die " cant open fasta_output";# > clears before each use 
    #this submits request and put results in file
    print FASTA `$submitString`;
    close FASTA;
    
    my $fasta_output_clean = "$fasta_output"."_clean";
    `grep -F '>' $fasta_output > $fasta_output_clean`;

    my %CE_acc;
    my %newSwiss;
    my %newTrEMBL;
    my %TrEMBLnew;
    my $foundID;
    open (FH, "<$fasta_output_clean");
    open (NEWIDS ,">/wormsrv2/tmp/found_IDs");
      while(<FH>)
	{
	  if( $_ =~ m/((AA\p{IsUpper})|(CAD))\d{5}\.\d*/ )   #matched AAM23424.2 or CAD23525.1 stylee
	    {
	      $aaaID = $&;
	      if( $_ =~ m/^>(\w+)/ )
		{   #get ID 
		  $foundID = $1;
		
		  if ( defined($noSWALL{$aaaID}) )  #if the $aaaID is one we found earlier
		    {
		      # this is an AAMXXXXX.X style id that we know about
		      $CE = $noSWALL{$aaaID};
		  
		      #assign ID to peptide
		      if( $foundID =~ m/_CAEEL/ ) {
			$newSwiss{$CE} = $foundID;
		      }
		      else {
			if( $foundID =~ m/[OPQ]\d\w{4}/ ) {
			  $newTrEMBL{$CE} = $foundID;
			}
			else {
			  $TrEMBLnew{$CE} = $foundID;
			}
		      }
		      #put in to hash to check for IP's and put in to .ace file
		      $wormpep_acc{$CE} = $foundID;
		    }
		  else {
		    print  "Cant find $aaaID in hash\n";
		  }
		}
#	      else {
#		print NEWIDS "Cant find $aaaID in hash\n";
#	      }
#	      else{
#		print NEWIDS "Cant match anything in $_\n";
#	      }
	    }
	}
    close FH;
    
    print NEWIDS "Here are the new SwissProt entries\n-------------------------------------------\n";
    foreach $pep (keys %newSwiss)
      {
	print NEWIDS "$pep $newSwiss{$pep}\n";
      }
    print NEWIDS "Here are the new TrEMBL entries\n-------------------------------------------\n";
    foreach $pep (keys %newTrEMBL)
      {
	print NEWIDS "$pep $newTrEMBL{$pep}\n";
      }
    print NEWIDS "Here are the TrEMBLnew entries\n-------------------------------------------\n";
    foreach $pep (keys %TrEMBLnew)
      {
	print NEWIDS "$pep $TrEMBLnew{$pep}\n";
      }
    close NEWIDS;
  }
