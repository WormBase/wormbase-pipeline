#!/usr/local/bin/perl5.6.0 -w
use strict;
use Wormbase;
use Ace;

my $maintainer = "All";
my $rundate    = `date +%y%m%d`; chomp $rundate;
my $wmpep_ver = &get_wormbase_version();
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
	    print "protein $count $protein[$count] accn is $accn_holder\n";

	    $accession[$count] = $accn_holder;
	    #build hash of wormpep to accession no.
	    $wormpep_acc{$protein[$count]} = $accession[$count];
	    $count++;

	  }
	else
	  {
	    print LOG "no TR: SW: accn in $_\n";
	  }
	
#	if($count == 2000)
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

while ( $lower <= $last_ind ) 
  {
     @array_chunk = @accession[$lower .. $upper];
    outputToAce(\%wormpep_acc, \@array_chunk, \$ace_output, \$old_ip_output, \$errorLog);
    $lower += $chunk_size;
    $upper += $chunk_size;
  }

@array_chunk = @accession[$lower .. $last_ind];
#process the remainders
outputToAce(\%wormpep_acc, \@array_chunk, \$ace_output, \$old_ip_output, \$errorLog);

close OLD_IP_OUTPUT;
close ACE_OUTPUT;
$db->close();

print LOG "Files available -\n $temp_acefile\n
$temp_clearOldInterPro\n";
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

    #create the fasta record in wormpep dir (and remove it after use)
    my $ver = &get_wormbase_version();
    print "\n\n$submitString\n";
    my $fasta_output = "~/fasta_output";#"/wormsrv2/WORMPEP/wormpep$ver/fasta_output";
    open (FASTA,">$fasta_output")||die " cant open fasta_output";# > clears before each use 
    #this submits request and put results in file
    print FASTA `$submitString`;
    close FASTA;
    
    #build hash of accession : Swissprot/Trembl ID
    open (FASTA,"<$fasta_output")|| die "cant open fasta record";#read
    #print FASTA "This is a temp file created by SubGetSwissId.pl and can be removed -unless script is running!";
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
	if( $databaseID =~q /CAEEL/){
	  $Database = "SwissProt";
	}
	if(defined($databaseID))
	  {
	   #remove any old InterPro entries for any proteins with SP / Tr Id's
	    print $old_ip_output "Protein : \"WP:$peptide\"\n";
	    print $old_ip_output "-D Motif_homol\t\"INTERPRO:\"\n";
	    print $old_ip_output "\n";#object separator

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
