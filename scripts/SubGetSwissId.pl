#!/usr/local/bin/perl5.6.0 -w
use strict;
use Wormbase;
use Ace;
my $maintainer = "All";
my $wmpep_ver = &get_wormbase_version();
$wmpep_ver = 78;#while wormpep.table broken

my $wormpepdir = "/wormsrv2/WORMPEP/wormpep$wmpep_ver";

#error log (email to All)
my $log = "$wormpepdir/error_report";
open(LOG,">$log")|| die "cant open $wormpepdir/error_report";
print LOG "# SubGetSwissId\n";
print LOG "\n";
print LOG "Wormpep version  :$wmpep_ver\n\n";
print LOG "Only using 78 unitl latest build version is correct" if ($wmpep_ver == 78);
print LOG "=============================================\n";
print LOG "\n";

my $errorLog = *LOG;

my $temp_acefile = "$wormpepdir/SwissprotIDs.ace";
#create temp ace file to output to |open and close to ensure new file (as output of data appends)
open (ACE_OUTPUT,">$temp_acefile") || die "cant open $temp_acefile";
close ACE_OUTPUT;
open (ACE_OUTPUT,">>$temp_acefile" || die "cant open $temp_acefile");
my $ace_output = *ACE_OUTPUT;

my @protein;
my @accession;
my @proteinID;
my %idextract;
my %wormpep_acc;
my $count;
$count = 0;





#open file with linkin data
open (INPUT, "$wormpepdir/wormpep.table$wmpep_ver")|| print "SubGetSwiss.pl cant find file wormpep.table$wmpep_ver";

my $db = Ace->connect(-path  =>  '/wormsrv2/pepace') || print "Couldn't connect to pepace\n", Ace->error;
my $accn_holder;
while (<INPUT>)
  {
    chomp;
    if ($_ =~ m /(CE\d+)/)#find wormpep id
      {
	$protein[$count] = $1;

	#get the accession no (this is done before checking pepace so that $_ is not the grep output#
#	print "$_\n";
	if($_ =~ m/(SW:\w+|TR:\w+)/)#any "words" prefixed with SW: or TR:
	  {
	    $accn_holder = $1;
	   # print "protein $count $protein[$count] accn is $accn_holder\n";
	  }
	else{
	  print LOG "protein $protein[$count] has no Swissprot / TrEMBL acc\n";}

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
	    
	if($count == 200)
	  {
	    #DEBUG########################################
	    #print "wormpep_acc hash after building\n";
	    #foreach my $worm (keys %wormpep_acc){
	    #  print "$worm","\t",$wormpep_acc{$worm},"\n";}
	    #############################################
	    
	    outputToAce(\%wormpep_acc, \@accession, \$ace_output, \$errorLog);

	    #reset batch loop variables
	    $count = 0;
	    %idextract = ();
	    @accession = "";
	    %wormpep_acc = ();
	    @protein = "";
	    @accession = "";
	    @proteinID = "";	    

	    last;
	  }
      }
  }
#process the remainders
outputToAce(\%wormpep_acc, \@accession, \$ace_output, \$errorLog);

#update pepace
#my @updates = $db->parse_file($temp_acefile);#
#print "@updates";
print LOG "SubGetSwissId finished at ",`date`,"\n";

#### use Wormbase.pl to mail Log ###########
my $name = "SubGetSwissId";
#$maintainer = "ar2\@sanger.ac.uk";
&mail_maintainer ($name,$maintainer,$log);
#########################################

sub outputToAce #(\%wormpep_acc, \@accession \$ace_output, \$errorLog)
  {
    my ($wormpep_acc, $accession, $ace_ouput, $errorLog) = @_;
    #$ace_output is filehandle of file creatred in calling script

    #construct pfetch command string by concatenating accession no's
    my $submitString = "pfetch "." @accession";#get full Fasta record (includes Interpro motifs)
    my %idextract;
    
    my $fasta_output = "../fasta_output";
    open (FASTA,">$fasta_output")||die " cant open fasta_output";# > clears before each use 
    print FASTA `$submitString`;
    close FASTA;
    
    #build hash of accession : Swissprot/Trembl ID

    open (FASTA,"<$fasta_output")|| die "cant open fasta record";#read
    #print FASTA "This is a temp file created by SubGetSwissId.pl and can be removed -unless script is running!";
   
    my $proteinID;
    print "opening FASTA file  . . . \n";
    while(<FASTA>)
      {
	chomp;
#	print "$_\n";
	if($_ =~ m/>(\w+)/)
	  {
	    $proteinID = $1;
	    my @pfetch = split(/ |>/,$_);
	    my $acc = $pfetch[2];
	    $idextract{$acc} = $proteinID ;
	  }	   
      }		    

    close FASTA;

    print "\n";
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
	elsif($Database eq "TR"){$Database = "Trembl";}
	else{print $errorLog $peptide,"\t",$full_accn,"\tfrom unknown database\n";}
	
	my $databaseID = $idextract{$accn};
	if(defined($databaseID)){
	  print $ace_output "Protein : \"WP:$peptide\"\n";
	  print $ace_output "Database \"$Database\" \"$accn\" \"$databaseID\"\n\n";}
      }
  }
