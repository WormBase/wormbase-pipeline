#!/usr/local/bin/perl
#
# maintain_EMBL_entries 
# v2.1 
# dl
#
# Script to generate unfinished data EMBL entries from db2fasta dumped fasta files
#
#
#
# 980604 dan : PP version.
# 980918 dan : Incorporated SQUIRREL-like database file for keeping better track of projects.
# 980924 dan : Customise for P.falciparum usage.
# 990315 dan : Upgrade to use pace ACEDB database (~/acedb/pace-mal).
# 990315 dan : Improve timestamp handling, no longer need to add EMBL format in command line.
# 990315 dan : Add dbxref:taxon line.
# 990804 dan : Major reworking of tace data and generation of the EMBL file
# 990804 dan : CC line shows date of initial submission

$ENV{'BABEL'} = "/nfs/disk222/malpub";
$ENV{'CURATOR'} = "dl1\@sanger.ac.uk";

require glob('~malpub/analysis/scripts/pf_library.pl'); 
require glob('~sjj/scripts/sanger.pl');
require 'ctime.pl';

#############################################
# command line switch for genomic sequences #
#############################################

while ($ARGV[0] =~ /^-/) {
    $_=shift;
    if (/^-s(.*)/) {
        $single=1;
	$choosen = shift;
    }
    elsif (/^-d(.*)/) {
        $debug = 1;
    }
    else {
        &usage;
    }
}

sub usage {
    print "maintain_EMBL_entries [-options] date\n";
    print " -single x : only parse a single chromosome\n\n";    
    exit;
}


############################
# Debug header information #
############################

if ($debug == 1) {

    print "\nmaintain_EMBL_entries v2.1\n\n";
    if ($single == 1) {
	print "Single splitdb mode\n";
	print "Splitdb = $choosen\n\n";
    }
}

######################
# set some variables #
###################### 

$spacer = ("n" x 800);                              # spacer length between contig sequences
%month_num2txt  = (
	   '01', 'JAN',
	   '02', 'FEB',
	   '03', 'MAR',
	   '04', 'APR',
	   '05', 'MAY',
	   '06', 'JUN',
	   '07', 'JUL',
	   '08', 'AUG',
	   '09', 'SEP',
	   '10', 'OCT', 
	   '11', 'NOV', 
	   '12', 'DEC', 
	   );


########################
# Calculate time stamp #
########################

if ($debug ==1) {
    print "calculate timestamp ...\n";
}
&get_timestamp;

########################
# Calculate date stamp #
########################

if ($debug ==1) {
    print "calculate datestamp ...\n";
}
&get_datestamp;

#########################################
# tace query for active splitdb details #
#########################################

#if ($debug ==1) {
#    print "get splitdb details ...\n";
#}
#&get_active_splits;

###################################
# generate EMBL files for splitdb #
###################################

if ($debug ==1) {
    print "\nwrite EMBL files ...\n";
}
&and_make_embl_files;

exit;

#############################
# get split data from ACEDB #
#############################

sub get_active_splits {
    
    $ENV{'ACEDB'}="/nfs/disk222/malpub/acedb/pace-mal/";
    local ($exec);
    local (*textace);
    local ($command);
    local ($split);
    local ($length);
    $exec=&tace;
    $command=<<EOF;
    find Active_Sequence
    show
    quit
EOF

    open(textace, "echo '$command' | $exec  | ");
    while (<textace>) {
	if (/^Sequence\s+(\S+)/) {
	    $split=$1; 
	    push (@databases, $split); 
	    $active_splits++;
	    if ($debug == 1) {
		print "$split\n";
	    }
	}
	if (/Finisher\s+(\S+.+)/) {
	    $finisher{$split}=$1;
	}
	if (/From_Laboratory\s+(\S+)/) {
	    $lab{$split}=$1;
	}
	if (/DB_info\s+Database\s+EMBL\s+(\S+)\s+(\S+)/) {
	    $EMBL_id{$split}=$1;
	    $EMBL_acc{$split}=$2;
	}
	if (/Chromosome\s+(\d+)/) {
	    $chromosome{$split}=$1;
	}
	if (/Creation_date\s+(\S+)/) {
	    ($year,$month,$day)=split(/-/,$1);
	    $creation_date{$split}=$day . "-" . $month_num2txt{$month} . "-" . $year;
	}
    }
                       
}


sub and_make_embl_files {
    
    if ($debug == 1) {
	print "make EMBL file for $choosen\n";
    }  
    $split = $choosen;
    &make_embl_file;
}

#####################
# make_embl_entries #
#####################

sub make_embl_file {

    $dir = "/nfs/disk100/wormpub/tmp";
    $seq = "";

    open (DNA, "<$dir/$split.shg") || die "Can't open file $split.fasta\n\n";
    while (<DNA>) {
	chomp;
	if (/>/) {
	    $seq .= $spacer;
	    next;
	}
	$seq .= $_;
	
    }
    close (DNA);
    
    # the concatenation routine will add a spacer (800 n's) 5' of the first contig
    # hence you need to strip them out
    
    $seq_format = substr ($seq,800);
    $seq = $seq_format;
    
    # calculate composition statistics for the string
    
    $length = length ($seq);
    $_ = $seq;
    $g = tr /g/g/;
    $a = tr /a/a/;
    $t = tr /t/t/;
    $c = tr /c/c/;
    $n = tr /n/n/;
    
    # make the embl file
 
    if ($debug == 1) {
	print "$split\n";
    }
    &embl_header;
    system ("cp $dir/${split}_${datestamp}.embl /nfs/disk100/wormpub/analysis/TO_SUBMIT/");
    
    # reset variables for next clone
    
    $seq="";
    $g=0; $a=0; $t=0; $c=0; $n=0;

}




###########################
# embl header information #
###########################


sub embl_header {
    #open (EMBL_OUT, ">/nfs/disk100/wormpub/tmp/${split}_${datestamp}.embl");
    print  "ID   CE${split} standard; DNA; INV; $length BP.\n";
    print  "XX\n";
    if ($EMBL_acc{$split} ne "") {print  "AC   $EMBL_acc{$split};\nXX\n";}
    print  "DE   Caenorhabditis elegans DNA *** SEQUENCING IN PROGRESS *** from clone $clone\n"; 
    print  "XX\n";
    print  "KW   HTG; HTGS_PHASE1.\n";
    print  "XX\n";
    print  "OS   Caenorhabditis elegans\n";
    print  "OC   Eukaryota; Metazoa; Nematoda; Secernentea; Rhabditia; Rhabditida;\n";
    print  "OC   Rhabditina; Rhabditoidea; Rhabditidae; Peloderinae; Caenorhabditis.\n";
    print  "XX\n";
    print  "RN   [1]\n";
    print  "RP   1-$length\n";
    print  "RA   Sulston J.;\n";
    print  "RT   ;\n";
    print  "RL   Submitted ($timestamp) to the EMBL/GenBank/DDBJ databases.\n";
    print  "RL   Nematode Sequencing Project, Sanger Centre, Hinxton, Cambridge CB10 1RQ, UK\n";
    print  "RL   and Department of Genetics, Washington University, St. Louis, MO 63110,\n";
    print  "RL   USA. E-mail: jes\@sanger.ac.uk or rw\@nematode.wustl.edu\n";
    print  "XX\n";
    print  "CC   IMPORTANT: This sequence is unfinished and does not necessarily\n";
    print  "CC   represent the correct sequence. Work on the sequence is in progress and\n";
    print  "CC   the release of this data is based on the understanding that the\n";
    print  "CC   sequence may change as work continues. The sequence may be contaminated\n";
    print  "CC   with foreign sequence from E.coli, yeast, vector, phage etc.\n";
    print  "XX\n";
    print  "CC   Order of segments is not known; 800 n's separate segments.\n";
    print  "XX\n";
    print  "FH   Key      Location/Qualifiers\n";
    print  "FH\n";
    print  "FT   source          1..$length\n";
    print  "FT                   /db_xref=\"taxon:6239\"\n";
    print  "FT                   /organism=\"Caenorhabditis elegans\"\n";
    print  "FT                   /chromosome=\"$chromosome{$split}\"\n";
    print  "FT                   /clone=\"$split\"\n";
    print  "XX\n";
    print  "SQ   Sequence $length BP; $a A; $c C; $g G; $t T; $n other;\n";

    
    $no_lines = int ($length / 60) +1; 
    for ($i = 0; $i < $no_lines; $i++) {
	$linestart = $i * 60;
	$newline = "     " . substr($seq,$linestart,10) . " " . substr($seq,($linestart+10),10) . " " .
	    substr($seq,($linestart+20),10) . " " . substr($seq,($linestart+30),10) . " " .
		substr($seq,($linestart+40),10) . " " .  substr($seq,($linestart+50),10);
	print  "$newline\n";
    }

    print  "//\n\n";


}

