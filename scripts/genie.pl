#!/usr/local/bin/perl
#
# Genie v1.0
# dl
# 991102
#
# Script to interactively update an ACEDB using a new dna sequence.
# Requires an open xace running 4.8b or above. A second xace is then 
# opened by this script with the -remote option such that DNA can be 
# extracted from the active zone and a new temp_gene inserted into this
# region.
#
# Caveats : 
# [1] !! THIS ONLY WORKS FOR DNA V DNA !!
# [2] Active zone must contain ALL of the dna query
#
###########################################################################
# 991102 dan : PP version
# 991124 dan : est_genome version
# 991125 dan : release v1.0

########################################
# command-line parsing                 #
########################################

while ($ARGV[0] =~ /^-/) {
    $_=shift;
    if (/^-d(.*)/) {
        $debug=1;
    }
    else {
        &usage;
    }
}
my $dna = shift;
if ($dna eq "") {&usage;}

########################################
# usage subroutine                     #
########################################

sub usage {
    print "Usage: genie.pl [-options] <dna file>/n/n";
    print "Options:\n";
    print "-d   Debug/Verbose mode\n\n";
    exit;
}

#########################################
# set path of xace executable           #
#########################################

$acepath = glob ("~wormpub/acedb/ace4/bin.ALPHA_4");

#########################################
# MainLoop                              #
#########################################

if ($debug == 1) {print "Genie v0.2\n\nMatching DNA sequence $dna\n";}

#########################################
# get sequence from the active region   #
#########################################
# xace -remote command to retrieve the active region

system ("$acepath/xremote.4_8b -remote \"gif seqget \; seqdna\" > active.dna");

if ($debug == 1) {print "Extacted sequence from open xace\n";}

#########################################
# mend sequence output                  #
#########################################

# we need to modify the output from xace to
# [A] rename the fasta header such that it is only the DNA Source
#     !! Also parse the coordinates for later reference !!

open (OUT, ">target.dna") || die "Can't open modified output\n";

open (SEQFILE, "active.dna") || die "Can't open dna output\n";
while (<SEQFILE>) {
    if (/>(\S+)/) {
	($project,$coords) = split (/\//, $1);
	($start,$stop) = split (/-/,$coords);
	print OUT ">$project\n";
    }
    elsif (/\/\//) {
	next;
    }
    else {
	print OUT "$_";
    }
}
close (SEQFILE);
close (OUT);

if ($debug == 1) {print "Modified active sequence ...\nDNA Source : \"$project\"\nCoords : $start -> $stop\n";}

########################################
# tidy up the active dna file          #
########################################

unlink "active.dna";
if ($debug == 1) {print "Removed old dna sequence file\n";}

########################################
# do est_genome job                    #
########################################

system ("est_genome -minscore 100 -genome target.dna -est $dna > e2g.ace");

if ($debug == 1) {print "est2genome prediction complete\n";}

########################################
# mend est_genome output               #
########################################

# we need to modify the est_genome output to 
# [A] reflect the correct subsequence coordinates in the full sequence
# [B] add 'hand_built' method object to the subsequence object

$exon = 1;

open (OUT_GW, ">newgene.ace")  || die "Can't open modified ace output\n";
open (ESTGENOME, "e2g.ace") || die "Can't open ace output\n";
while (<ESTGENOME>) {
   
    if (/Span\s+\d+\s+100.0\s+(\d+)\s+(\d+)\s+\S+\s+\d+\s+\d+\s+\S+/) {
	
	print OUT_GW "Sequence $project\n";
	
	if ($start < $stop) {
            print "forward strand\n";
            # forward strand
            my $sub_start = $start + $1 - 1;
            my $sub_stop  = $start + $2 - 1;
            print OUT_GW "Subsequence temp_GENIE_$project.1 $sub_start $sub_stop\n\n";
        }
        else {
            print "reverse strand\n";
            # reverse strand
            my $sub_start = $start - $1 + 1;
            my $sub_stop  = $start - $2 + 1;
            print OUT_GW "Subsequence temp_GENIE_$project.1 $sub_start $sub_stop\n\n";
        }
    }
    
    if (/Segment\s+\d+\s+100.0\s+(\d+)\s+(\d+)\s+\S+\s+\d+\s+\d+\s+\S+/) {
	if ($exon == 1) {
	    $exon++;
	    $gene_start = $1;
	    $sub_start = $1 - $gene_start + 1;
            $sub_stop  = $2 - $gene_start + 1;
	    print OUT_GW "Sequence : \"temp_GENIE_$project.1\"\nSource_Exons $sub_start $sub_stop\n";
	}
	else {
	    $sub_start = $1 - $gene_start + 1;
            $sub_stop  = $2 - $gene_start + 1;
	    print OUT_GW "Source_Exons $sub_start $sub_stop\n";
	}
    }

}
close (ESTGENOME);

print OUT_GW "\nSequence : \"temp_GENIE_$project.1\"\n";
print OUT_GW "Source \"$project\"\n";
print OUT_GW "Properties Coding CDS\n";
print OUT_GW "Method hand_built\n";
close (OUT_GW);

if ($debug == 1) {print "Modified est_genome output\n";}

#########################################
# tidy up est_genome and dna sequence   #
#########################################

unlink "e2g.ace";
unlink "target.dna";

if ($debug == 1) {print "Removed original est_genome output\nRemoved target dna\n";}

#########################################
# parse ace file into the active xace   #
#########################################

system ("$acepath/xremote.4_8b -remote \"parse newgene.ace \"");
system ("$acepath/xremote.4_8b -remote \"gif seqget $project \; seqdisplay \"");

if ($debug == 1) {print "Parsed ace file into open xace\n";}

#########################################
# tidy up ace file                      #
#########################################

unlink "newgene.ace";

if ($debug == 1) {print "hasta luego\n\n";}
exit;

#########################################
# end of file                           #
#########################################


