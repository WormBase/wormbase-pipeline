#!/usr/local/bin/perl
#
# A perl function library for the Sanger Centre
#
# 
#
#
#
#
# Steven Jones, Sanger Centre 1995


###############################################################################
#Calculate which tace binary to use
###############################################################################

sub tace {
   local($prog);
   local($name); 
   $name=`uname -sr`;
    
    if ($name=~/^SunOS/) {($prog)=<~wormpub/acedb/ace4/bin.SUN_4/tace>;}
    elsif ($name=~/^IRIX/) {($prog)=<~wormpub/acedb/ace4/bin.SGI_4/tace>;}
    elsif ($name=~/^OSF/)  {($prog)=<~wormpub/acedb/ace4/bin.ALPHA_4/tace>;}
    elsif ($name=~/^Linux/)  {($prog)=<~wormpub/acedb/ace4/bin.LINUX/tace>;}

    else {print STDERR "No known binary for $uname\n";exit;}

    return $prog;
}

################################################################################
#A Function which assigns the home directory of the user to $tilda
################################################################################

sub tilda  {
	#A Function which assigns the home directory of the user to $~
	local($name,$passwd,$uid,$gid,$quota,$comment,$gcos,$dir)=getpwnam($_[0]);
	$tilda=$dir;
	return $tilda;
	  }

###############################################################################
# A function which returns an array containing all annotated cosmids
##############################################################################

sub annotcosmids {
	local($exec);
        local($command);
	local(@array);
	undef @array;
	$exec=&tace;
	$ENV{'ACEDB'}="/nfs/disk100/wormpub/acedb/ace4/cam";     
	$command=<<EOF;
	Query Find genome_sequence annotated
	list -f - 
	quit
EOF
open(textace, "echo '$command' | $exec | ");
while (<textace>) {
  if (/^Sequence\s+:\s+\"(\w+)\"/)  {if ($1 ne "?Sequence" &&  $1!~/YAC_*/) {push (@array, $1)}}
	 	 }
	close textace;     
	return @array;
        }

###############################################################################
#A function which retuurns an array of all cosmids to be sequenced (including St. louis) 
###############################################################################

sub canonical_cosmids {

        local($exec);
        local($command);
	local(@array);
	undef @array;
	$exec=&tace;
	$ENV{'ACEDB'}="/nfs/disk100/wormpub/acedb/ace4/cgc";     
	$command=<<EOF;
	find genome_sequence
       	list -f -
	quit
EOF
open(textace, "echo '$command' | $exec  | ");
while (<textace>) {
  if (/^Sequence\s+:\s+\"(\w+)\"/)  {if ($1 ne "?Sequence" &&  $1!~/YAC_*/) {push (@array, $1)}}
	 	 }
	close textace;     
	return @array;
        }


###############################################################################
# A function to return the finisher of a given cosmid 
###############################################################################

sub finisher_cosmid {
       
        local($finisher);
        $ENV{'SQUIRREL'}="/nfs/disk39/squirrel/squirrel";
	open(line,"/usr/local/badger/bin/sq-lookup  worm databases $_[0] |");

	while (<line>) {if (/^\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)$/) { $finisher=$1;}
			}
	close line;
	return $finisher;
}



###############################################################################
#A function to return the chromosomal position of a cosmid clone
###############################################################################



sub celechrom {
	$ENV{'ACEDB'}="/nfs/disk100/wormpub/acedb/ace4/cgc";
	local ($exec);
        local (*textace);
	local ($command);
        local ($chromosome);
	$exec=&tace;
	$chromosome="";
	$command=<<EOF;
	find clone $_[0]
        follow pMap
        follow Map
	list -f -
	quit
EOF

open(textace, "echo '$command' | $exec   | ");
while (<textace>) {if (/^Map\s+:\s+\"(\w+)\"$/)  {$chromosome=$1;}}

#If cannot find cosmid location - it could be buried, funny match to, or exact 
#match to
if ($chromosome eq "") {
$command=<<EOF;
	find clone $_[0]
        follow Approximate_Match_to
        follow pMap
        follow Map
        list -f -
	find clone $_[0]
        follow Funny_Match_to
        follow pMap
        follow Map
        list -f -
        find clone $_[0]
        follow Exact_Match_to
        follow pMap
        follow Map
        list -f - 
	quit
EOF
open(textace, "echo '$command' | $exec  | ");
while (<textace>) {if (/^Map\s+:\s+\"(\w+)\"$/)  {$chromosome=$1;}}
}
close textace;
#return 0 if chromosome is still not found
if ($chromosome ne "") {return $chromosome;} else {return "0"};
}


###############################################################################
#A function to return an array of all the finished cosmids from cam only 
###############################################################################

sub fincosmids {
        local (@array);
        undef @array;
	local($dir);
	$dir=&tilda(wormpub);
	open(XX123,"$dir/analysis/cosmids/current.versions");
	while (<XX123>)
          {if (/^(\S+)\//) {push (@array, $1);}
           }  
        close XX123;
        return  @array;
}

###############################################################################
#A function to return an array of all the finished cosmids from stl and cam
#this is used to updating the embl sequences. 
###############################################################################

sub finallcosmids {

        local($exec);
        local($command);
	local(@array);
	undef @array;
	$exec=&tace;
	$ENV{'ACEDB'}="/nfs/disk100/wormpub/acedb/ace4/cgc";     
	$command=<<EOF;
	find genome_sequence
	where DNA  AND species="Caenorhabditis elegans"
       	list -f -
	quit
EOF
open(textace, "echo '$command' | $exec  | ");
while (<textace>) {
  if (/^Sequence\s+:\s+\"(\w+)\"/)  {if ($1 ne "?Sequence" &&  $1!~/YAC_*/) {push (@array, $1)}}
	 	 }
	close textace;     
	return @array;
        }



###############################################################################
# A function to return an array of the sequence of an  C. elegans cosmid 
###############################################################################

sub celeseq {
        local($cosmid);
	local(@cosmid);
	local(*file);
        local($dir);
	local($_);
	$dir=&tilda(ftp);
	open(file,"$dir/pub/databases/C.elegans_sequences/FINISHED_SEQUENCES/UNSORTED/$_[0].seq");
	$/="";
	while (<file>) {
	             s/>.+\n//g;
                     s/\s+//g;
	             s/\n//g;
	             $cosmid=$_;
		                 }
	$/="\n";
        @cosmid=split("",$cosmid);
	#print "@cosmid\n";
close file;
return @cosmid;


}
###############################################################################
#A function to get accession number for a C. elegans cosmid
###############################################################################

sub celeaccession {
        local($exec);
        $exec=&tace;
        local($command);
	local($accession);
	$ENV{'ACEDB'}="/nfs/disk100/wormpub/acedb/ace4/cam";     
	$command=<<EOF;
	find sequence $_[0]
	show DB_info
 	quit
EOF
open(textace, "echo '$command' | $exec  | ");

while (<textace>) {if (/\s+Database\s+EMBL\s+\S+\s+(\S+)\n/) {$accession=$1;}
              }
	close textace;     
	return $accession;
	
        }

###############################################################################
#A function to return the number of CDS (genes) in the database derived from
#C. elegans genomic sequence
#only gets number from submitted cosmids 
###############################################################################

sub getnumberofgenes {
        local($total);
        local($exec);
        local($command);
	$ENV{'ACEDB'}="/nfs/disk100/wormpub/acedb/ace4/cam";
	$exec=&tace;
        $command=<<EOF;
	find genome_sequence
        where DB_info AND NEXT AND NEXT AND NEXT AND NEXT
        follow subsequence
        where CDS
        list -f -
        find sequence *LINK*
        follow subsequence
        where  CDS
        list -f -
        find sequence *LINK*
        follow subsequence
        follow subsequence
        where CDS
        list -f -
        quit
EOF
      #The sed command is used to handle the multiple entries caused by
      #alternative transcripts so they are counted as just one gene
      #At the moment cannot handles genes with >2 alternates. 

      open(textace, "echo '$command' | $exec | sed 's/b\"/\"/' | sort -u |"); 
      while (<textace>) {   
                            #This line only counts genes which have been 
                            #analysed i.e. .number genes
                            if (/^Sequence\s+\:\s+\"\S+\.\d+/) {$total++;}

                        }
close textace;
return $total;
}

###############################################################################
#A function to determine the number of ESTs hitting the Cambridge dataset
###############################################################################

sub ESThits {
        local($exec);
	local($command);
	local($total);
	$ENV{'ACEDB'}="/nfs/disk100/wormpub/acedb/ace4/cam";
	$exec=&tace;
	$command=<<EOF;
	find genome_sequence
        follow subsequence
        follow matching_cDNA
        where cDNA
        list -f -
        find sequence *LINK*
        follow subsequence
        follow matching_cDNA 
        where cDNA
        list -f -
        find sequence *LINK*
        follow subsequence
        follow subsequence
        follow matching_cDNA
        where cDNA
        list -f -
        quit
EOF
        open(textace, "echo '$command' | $exec  | sort -u |");
	while (<textace>) {if (/^Sequence\s+\:/) {$total++;}
                          }
close textace;
return $total;
}

###############################################################################
#A function to reverse complement a DNA string 
###############################################################################

sub revcomp {
            local(@sequence);
	    local($output);
	    local($intron);
	    $intron=$_[0];
	    $intron=~tr/A-Z/a-z/;
	    $intron=~s/g/C/g;
	    $intron=~s/c/G/g;
	    $intron=~s/a/T/g;
	    $intron=~s/t/A/g;
	    @sequence=split(//,$intron);
	    $output=join('',reverse(@sequence));
	    return $output;
}

###############################################################################
#A function to return the number of coding bases in a cosmid 
###############################################################################


          #note returns -1 if coding percentage could not be determined
          #as currently this relies on the cosmid being in EMBL and being 
          #annotated 

sub codingbases {  
   local(*tempfile);
   local(*sequence);
   local(*textace);
   local($coding=0);
   local($check="");
   local($exec);
   local($command);

   open(tempfile,">/tmp/temp$$");
   open(sequence,"getz -t -d -sf EMBL  \"([emnew-id:CE$_[0]]) \| (([embl-id:CE$_[0]]) \! ([embl-id:CE$_[0]] < emnew))\" |");
   while (<sequence>) {$check=$check.$_;print tempfile;}
    close sequence;
    close tempfile;

   #some of the st. louis entries are in embl as CEU(cosmid name)
   if ($check eq "") {

   open(tempfile,">/tmp/temp$$");
   open(sequence,"getz -t -d -sf EMBL  \"([emnew-id:CEU$_[0]]) \| (([embl-id:CEU$_[0]]) \! ([embl-id:CEU$_[0]] < emnew))\" |");
    while (<sequence>) {$check=$check.$_;print tempfile;}

    close sequence;
    close tempfile;
}
   
    open(tempfile,"cat /tmp/temp$$ |");
    while (<tempfile>) {if (/^FT\s+\/translation=\"(\S+)$/) {$trstart="yes";}
			if ($trstart eq "yes") {s/^FT\s+//;
						s/\/translation=\"//;
                                                #print "sequence: $_";
						$sequence=$sequence.$_;
                                                }
			if (/\"$/) {$trstart="stop";}
		    }

    close tempfile;
    unlink "/tmp/temp$$";
   $sequence=~s/\"//;$sequence=~s/\s+//g;$sequence=~s/\n//g;
   #print "$sequence\n";
   $coding=length($sequence)*3;
    if ($check eq "") {print STDERR "$_[0] not found in embl ignoring\n";return "-1";}
                     else {return $coding;}
}

#############################################################
#return introns and exons 
#############################################################

sub intronexon {
 
#note NEED TO HANDLE ALT SPLICED GENES BETTER HERE 
#will need to change this so it returns arrays instead of printing 
#this only assays all cosmids 

        local(@introns);
        local($exec);
        $exec=&tace;
        local($command);
        local($accession);
        $ENV{'ACEDB'}="/nfs/disk100/wormpub/acedb/ace4/cam";     
        $command=<<EOF;
        find sequence $_[0].*
        where CDS
        show Source_Exons
        quit
EOF
open(textace, "echo '$command' | $exec  | ");

while (<textace>) {if (/^\s+(\d+)\s+(\d+)/) {

    print "$_[0]:exon ",$2-$1,"\n";
    #calc intron
    if ($1=="1") {$newgene="yes";} else {$newgene="no";}
    if ($newgene eq "no") {$intron=$1-$lastexonend; push(@introns,$intron);print "$_[0]:intron $intron\n";
                           $intronsize=$intronsize+($1-$lastexonend);$introns++;}
    $lastexonend=$2;
}
                   
               }
    close textace;


}

###############################################################################
#returns the number of inverted repeats 
###############################################################################

sub invcontent {

       $ENV{'ACEDB'}="/nfs/disk100/wormpub/acedb/ace4/cam";     
       local($command);
       local(*textace);
       local($totinv);
       local($exec);
       $exec=&tace;
       $command=<<EOF;
        find sequence $_[0]
        show -a feature
        quit
EOF
open(textace, "echo '$command' | $exec  | ");

    while(<textace>) {if (/inverted\"\s+(\d+)\s+(\d+)/) {$totinv=$totinv+0.5;}
						       
                     }
    close textace; 
    return $totinv;

}


###############################################################################
#returns the number of tandem repeats
###############################################################################

#Note: there may be some tandem arrays which overlap this function
#should be changed at some point to account for this. 


sub tancontent {

       $ENV{'ACEDB'}="/nfs/disk100/wormpub/acedb/ace4/cam";     
       local($command);
       local(*textace);
       local($tottan);
       local($exec);
       $exec=&tace;
       $command=<<EOF;
        find sequence $_[0]
        show -a feature
        quit
EOF
open(textace, "echo '$command' | $exec  | ");

    while(<textace>) {if (/tandem\"\s+(\d+)\s+(\d+)/) {$tottan++;}}
    close textace; 

    return $tottan;

}

################################################################################
#return percent GC
################################################################################

sub GCcontent {

    local(@seq);
    local($sequence);
    local($percentGC=0);
    local($cosmidsize);
    @seq=&celeseq($_[0]); 
    $cosmidsize=$#seq+1;
    $sequence=join('',@seq);
    $sequence=~s/A|T|N//g;$GC=length($sequence);
    if ($cosmidsize >1) {$percentGC=($GC/$cosmidsize)*100;} 
              else {print STDERR "something wrong with $_[0]\n";}
    return $percentGC;

}

###############################################################################
#return the number of repeat family elements for a cosmid sequence 
###############################################################################


sub repfamcontent {

       $ENV{'ACEDB'}="/nfs/disk100/wormpub/acedb/ace4/cgc";     
       local($command);
       local(*textace);
       local($reptot=0);
       local($exec);
       $exec=&tace;
       $command=<<EOF;
       find sequence $_[0]
        show -a homol
        quit
EOF
open(textace, "echo '$command' | $exec  | ");

    while(<textace>) {if (/^Motif_homol\s+\"(\S+)\"/ && ($1 eq $_[1])) {$reptot++;}
                     }
    close textace; 
    return $reptot;

}

################################################################################
#return chromosome and position of a cosmid sequence 
################################################################################

sub seq_pos {

     $ENV{'ACEDB'}="/nfs/disk100/wormpub/acedb/ace4/cgc";     
       local($command);
       local(*textace);
       local($chrome);
       local($pos);
       local($exec);
       $exec=&tace;
       $command=<<EOF;
       find sequence $_[0]
        show map
        quit
EOF
open(textace, "echo '$command' | $exec  | ");
     while(<textace>) {if (/^\s+Sequence-(\S+)\s+Ends\s+Left\s+(\d+)/) {$chrome=$1;$pos=$2;}
                     }
    close textace; 

     return ($chrome,$pos);
}


###############################################################################
#subroutine to return the briefid of a gene from cgcace  
###############################################################################

sub brief_id {

    $ENV{'ACEDB'}="/nfs/disk100/wormpub/acedb/ace4/cgc";     
       local($command);
       local(*textace);
       local($brief_id);
       local($exec);
       $exec=&tace;
       $command=<<EOF;
       find sequence $_[0]
        show -a brief_identification
        quit
EOF
open(textace, "echo '$command' | $exec  | ");
    while(<textace>) {if (/brief_identification\s+\"(.+)\"$/) {$brief_id=$1;}
                     }
    close textace; 

    return $brief_id;
}

sub get_team {
    local($command);
    local($team);
    local($exec);
    $exec=&tace;
    $ENV{'ACEDB'}="/nfs/disk54/badger/pace-worm";
    $command=<<EOF;
    Query Find Sequence $_[0]
    show -a Finishing_group
    quit
EOF
    open(textace, "echo '$command' | $exec  | ");
    while (<textace>) {if (/^Finishing_group\s+\"(\S+)\"/) {$team=$1;}
    
    }
    close textace;
    return $team;
}

######################################################################################
# take the name of a gene and returns the accession with the highest database match
# will also optionally return the highesthit of a species. 
#  e.g. ($tophit,$topscore,$topspecies)=&tophit("AH6.2","schizosaccharomyces  pombe");
# returns the name of the hit, the blastscore and the species name.
######################################################################################

sub tophit {


use Ace;
#use strict vars;

#make the species name lowercase and remove spaces.
my $opt_species;
$opt_species=$_[1];
$opt_species=~s/\s+//g;
$opt_species=~tr/A-Z/a-z/;

use constant HOST => $ENV{ACEDB_HOST} || 'ics1a.sanger.ac.uk';
use constant PORT => $ENV{ACEDB_PORT} || 100100;
$|=1;

my $db = Ace->connect(-host=>HOST,-port=>PORT) || die "Connection failure: ",Ace->error;
my $sequence= $db->fetch('Sequence',$_[0]);

my $wormpep=$sequence->at('Visible.Corresponding_protein')->pick;

my @hits=$wormpep->at('Homol.Pep_homol')->col(1);
my $hit;
my $highscore=0;
my $hitprotein;
my $matchspecies;
my $species;
my @scores;
my $score;
my $highesthit;

foreach $hit (@hits) {#print $hit->asString;
		      my $hitprotein = $db->fetch('Protein',$hit);	
		      ($species)=$hitprotein->at('Origin.Species');

		      #make a lower case non-whitspace species for easy matching		      
		      $matchspecies=$species;
		      $matchspecies=~s/\s+//g;
                      $matchspecies=~tr/A-Z/a-z/;

		      if ($opt_species ne "" && $opt_species ne $matchspecies) {next;} 

		      $species{$hit}=$species;

		      #print "Species $species\n";
		      @scores=$hit->right->col(1);
		      foreach $score (@scores) {
				if ($score > $highscore) {$highscore=$score;$highesthit=$hit;}
		       }
}

return($highesthit,$highscore,$species{$highesthit});

}

################################################################################
#Return a true value
################################################################################

1;

################################################################################ 












