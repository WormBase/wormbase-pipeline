#!/usr/local/bin/perl5.6.0 -w
use strict;
use Wormbase;
use Ace;

my $maintainer = "All";
my $rundate    = `date +%y%m%d`; chomp $rundate;
my $wmpep_ver = &get_wormbase_version();$wmpep_ver = 78;#while wormpep.table broken
my $wormpepdir = "/wormsrv2/WORMPEP/wormpep$wmpep_ver";
my $log = "/wormsrv2/logs/GetSwissIDandInterpro.WB$wmpep_ver.$rundate";#"$wormpepdir/error_report";#error log (email to All)
my $temp_acefile = "$wormpepdir/SwissprotIDs.ace";
my $temp_clearOldInterPro = "$wormpepdir/clearOldInterPros.ace";# file to clean up old InterPro domains from the database
my $pepace = "/wormsrv2/pepace";

my $ace_output = *ACE_OUTPUT;
my $old_ip_output =*OLD_IP_OUTPUT;

open(LOG,">$log")|| die "cant open $log";
my $errorLog = *LOG;
print LOG "#$0\n";
print LOG "\n";
print LOG "Wormpep version  :$wmpep_ver\n\n";
print LOG "Only using 78 until latest build version is correct\n" if ($wmpep_ver == 78);
print LOG "=============================================\n";
print LOG "\n";

#create temp ace file to output to |open and close to ensure new file (as output of data appends)
open (ACE_OUTPUT,">$temp_acefile") || die "cant open $temp_acefile";
close ACE_OUTPUT;
open (ACE_OUTPUT,">>$temp_acefile" || die "cant open $temp_acefile");

#create temp ace file to output to open and close to ensure new file (as output of data appends)
open (OLD_IP_OUTPUT,">$temp_clearOldInterPro") || die "cant open $temp_clearOldInterPro";
close OLD_IP_OUTPUT;
open (OLD_IP_OUTPUT,">>$temp_clearOldInterPro" || die "cant open $temp_clearOldInterPro");


my @protein;
my @accession;
my @proteinID;
my %idextract;
my %wormpep_acc;
my $count;
$count = 0;

#open file with linkin data
open (INPUT, "$wormpepdir/wormpep.table$wmpep_ver")|| print "SubGetSwiss.pl cant find file wormpep.table$wmpep_ver";

my $db = Ace->connect(-path  =>  $pepace) || print "Couldn't connect to $pepace\n", Ace->error;
my $accn_holder;
while (<INPUT>)
  {
    chomp;
    if ($_ =~ m /(CE\d+)/)#find wormpep id
      {
	$protein[$count] = $1;

	#get the accession no (this is done before checking pepace so that $_ is not the grep output#
	#print "$_\n";
	if($_ =~ m/(SW:\w+|TR:\w+)/)#any "words" prefixed with SW: or TR:
	  {
	    $accn_holder = $1;
	    #print "protein $count $protein[$count] accn is $accn_holder\n";

	    #check if the pepace entry has a SWISS-PROT entry
	    my $worm_protein = $db->fetch(Protein => "WP:$protein[$count]");
	    my @ace_database = $worm_protein->at('DB_info.Database');#returns list of entries in Database field eg "WORMPEP SwissProt"
	    
	    #if the protein does not have a SwissProt ID in Pepace put it on the the list of queries to be done in next batch;
	    my $swiss = grep $_  eq "SwissProt", @ace_database;
	    unless($swiss)
	      {
		$accession[$count] = $accn_holder;
		#build hash of wormpep to accession no.
		$wormpep_acc{$protein[$count]} = $accession[$count];
		$count++;
	      }	
	  }
	
	if($count == 2000)#limits the no of requests in pfetch call per loop
	  {
	    #DEBUG########################################
	    #print "wormpep_acc hash after building\n";
	    #foreach my $worm (keys %wormpep_acc){
	    #  print "$worm","\t",$wormpep_acc{$worm},"\n";}
	    #############################################
	    
	    outputToAce(\%wormpep_acc, \@accession, \$ace_output, \$old_ip_output, \$errorLog);

	    #reset batch loop variables
	    $count = 0;
	    %idextract = ();
	    @accession = "";
	    %wormpep_acc = ();
	    @protein = "";
	    @accession = "";
	    @proteinID = "";	    

	    #last;#only included for testing on small sample sets
	  }
      }
  }
#process the remainders
outputToAce(\%wormpep_acc, \@accession, \$ace_output, \$old_ip_output, \$errorLog);

close OLD_IP_OUTPUT;
close ACE_OUTPUT;
$db->close();

#update pepace
my $command=<<END;
pparse $temp_clearOldInterPro 
  save 
  pparse $temp_acefile
  save
  quit
END
  
  open (WRITEDB, "| tace $pepace >> $log") || die "cant do the tace command on $pepace";
print WRITEDB $command;
close WRITEDB;

#modification of database completed

#===========================================================================

print LOG "parsed files:\n $temp_clearOldInterPro \n $temp_acefile \n\nto database $db\n\n";
print LOG "SubGetSwissId finished at ",`date`,"\n";

close LOG;
#### use Wormbase.pl to mail Log ###########
my $name = "Update Swiss ID's and Interpro motifs";
#$maintainer = "ar2\@sanger.ac.uk";
&mail_maintainer ($name,$maintainer,$log);
#########################################

#===========================================================================
sub outputToAce #(\%wormpep_acc, \@accession \$ace_output, \$errorLog)
  {
    my ($wormpep_acc, $accession, $ace_ouput, $old_ip_output, $errorLog) = @_;
    #$ace_output is filehandle of file creatred in calling script

    #construct pfetch command string by concatenating accession no's
    my $submitString = "pfetch -F"." @accession";#get full Fasta record (includes Interpro motifs)
    my %idextract;
    
    my $fasta_output = "../fasta_output";
    open (FASTA,">$fasta_output")||die " cant open fasta_output";# > clears before each use 
    print FASTA `$submitString`;
    close FASTA;
    
    #build hash of accession : Swissprot/Trembl ID

    # ! ! ! ! ! CHANGING RECORD SEPARATOR ! ! ! ! ! ###
    $/ = "\/\/"; #//
    open (FASTA,"<$fasta_output")|| die "cant open fasta record";#read
    #print FASTA "This is a temp file created by SubGetSwissId.pl and can be removed -unless script is running!";
    my %interpro;
    my $proteinID;
    my $acc;
    print "opening FASTA file  . . . \n";
    while(<FASTA>)
      {
	chomp;
	#get the SWISSPROT/TrEMBL id
	if($_ =~ m/ID\s+(\w+)/)
	 { $proteinID = $1;}
	
	#get the SWISSPROT accession
	if($_ =~ m/AC\s+(\w+)/)
	  { $acc = $1;}

	#get all InterPro motifs
	$idextract{$acc} = $proteinID ;

	#extract Interpro motifs
	while ($_ =~ m/(IPR\d+)/g){
	  $interpro{$proteinID} .= "$1 ";}
	
	#print "$interpro{$proteinID}\n";
		    
      }
    close FASTA;
    # ! ! ! ! resetting record separator ! ! ! ! #
    $/ = '\n'; #\n

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
	if(defined($databaseID))
	  {
	   #remove any old InterPro entries for any proteins with SP / Tr Id's
	    print $old_ip_output "Protein : \"WP:$peptide\"\n";
	    print $old_ip_output "-D Motif_homol\t\"INTERPRO:\"\n";
	    print $old_ip_output "\n";#object separator

	    print $ace_output "Protein : \"WP:$peptide\"\n";
	    print $ace_output "Database \"$Database\" \"$accn\" \"$databaseID\"\n";
	    
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
