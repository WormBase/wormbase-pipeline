#!/usr/local/bin/perl5.6.0 -w
use strict;
use Wormbase;
use Ace;

my $wmpep_ver = &get_wormbase_version();
$wmpep_ver = 78;#while wormpep.table broken

my $wormpepdir = "/wormsrv2/WORMPEP/wormpep$wmpep_ver";

#error log (email to All)
open(LOG,">$wormpepdir/error_report");
print LOG $^T,"\n\n";#time


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
	    print "protein $count $protein[$count] accn is $accn_holder\n";
	  }
#	else{
	#  print "protein $count $protein[$count] has no Swissprot acc\n";}

	#check if the pepace entry has a SWISS-PROT entry
	my $worm_protein = $db->fetch(Protein => "WP:$protein[$count]");
	my @ace_database = $worm_protein->at('DB_info.Database');#returns list of entries in Database field eg "WORMPEP SwissProt"

	#if the protein does not have a SwissProt ID in Pepace put it on the the list of queries to be done in next batch;
	my $swiss = grep $_  eq "SwissProt", @ace_database;
	print "$swiss\n";
	unless($swiss)
	  {
      	    $accession[$count] = $accn_holder;
	    #build hash of wormpep to accession no.
	    $wormpep_acc{$protein[$count]} = $accession[$count];
	    $count++;
	  }
	#else{print "$worm_protein already has SwissProt\n";}
	    
	if($count == 2)
	  {
	    #DEBUG########################################
	    #print "wormpep_acc hash after building\n";
	    #foreach my $worm (keys %wormpep_acc){
	    #  print "$worm","\t",$wormpep_acc{$worm},"\n";}
	    #############################################
	    
	    outputToAce(\%wormpep_acc, \@accession, \$ace_output);

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
print $ace_output "$count\n";
outputToAce(\%wormpep_acc, \@accession, \$ace_output);

#update pepace
#my @updates = $db->parse_file($temp_acefile);
#print "@updates";



sub outputToAce #(\%wormpep_acc, \@accession \$ace_output)
  {
    my ($wormpep_acc, $accession, $ace_ouput) = @_;
    #construct pfetch command string by concatenating accession no's
    my $submitString = "pfetch"." @accession";
    my %idextract;
    
    open (FASTA,">fasta_output")||die " cant open fasta_output";# > clears before each use 
    print FASTA `$submitString`;
    close FASTA;
    
    #build hash of accession : Swissprot/Trembl ID
    open (FASTA,"<fasta_output")|| die "cant open fasta record";#read
    print "build idextract hash from A file\n";
    
    while(<FASTA>)
      {
	chomp;
	if($_ =~ m/>(\w+)/)
	  {
	    my $proteinID = $1;
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
     
    print "outputting to ace file $ace_output\n";

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
	else{print LOG $peptide,"\t",$full_accn,"\tfrom unknown database\n";}
	
	my $databaseID = $idextract{$accn};
	if(defined($databaseID)){
	  print $ace_output "Protein : \"WP:$peptide\"\n";
	  print $ace_output "Database \"$Database\" \"$accn\" \"$databaseID\"\n\n";}
	#	else{
	#  print LOG "$peptide : $accn has no databaseID ( $full_accn)\n";}
      }
    print "\n";
  }
