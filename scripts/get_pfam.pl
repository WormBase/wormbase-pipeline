#!/gsc/bin/perl 
#
##################
# get_pfam.pl    #
##################
#
# This script interogates an ACEDB database and returns
# all pfam/Interpro/blastx data as appropriate and 
# generates the DB_remark for stlace
#
# NOTE: there is a two-week lack of synchrony with acedb. 
#       new forms in camace not yet appearing in acedb 
#       get their similarities from the previous or "parent"
#       form of the gene (for alternative splicing)
#
#       ie, say B0025.1c is a new form, and not yet in acedb:
#           check B0025
#           if no hits, check B0025.1a
#
#
### DB_remark is generated as follows:  ###
#

###########################
# SEQUENCE OBJECTS        #
###########################
# Real Genes (clone.number)
# --------------------------
# 1. locus --> "C. elegans LOCUS protein"
# 2. pfam motif --> "contains similarity to Pfam domain PFXXXXXX (title)"
#       * take ALL pfam motifs
# 3. interpro motif --> "contains similarity to Interpro domain XXXXXXX (title)"
#       * take ALL interpro motifs, if no Pfam hits
# 4. blastx similarity --> "contains similarity to SPECIES TITLE; ID"
#       * take highest scoring non-worm hit WITH a title, if no Interpro hits
#

#########################
# PSEUDOGENE OBJECTS    #
#########################
# -----------
#  w/ locus: "C. elegans TEXT pseudogene LOCUS"
#       * text is from Pseudogene(1)
#  w/o locus: NOTHING
#

#########################
# TRANSCRIPT OBJECTS    #
#########################

##   fields returned from database:
##       $transcript1 = $obj->Transcript(1); # type of transcript
##       $transcript2 = $obj->Transcript(2); # text
#
# tRNAs
#	w/locus { "C. elegans tRNA $locus"; }
#	w/o locus { "C. elegans predicted tRNA"; }
#   
# misc_RNA genes
#	w/ locus { "C. elegans non-protein coding RNA $locus"; }
#	w/o locus { "C. elegans probable non-coding RNA"; }
#    
# snRNA genes
#	w/locus) { "C. elegans small nuclear RNA $transcript2 $locus"; }
#	w/o locus { "C. elegans small nuclear RNA $transcript2"; }
#    
# snoRNA genes
#	w/locus { "C. elegans small nucleolar RNA $transcript2 $locus"; }
#	w/o locus { "C. elegans small nucleolar RNA $transcript2"; }
#   
# miRNA genes
#	w/locus { "C. elegans microRNA $locus"; }
#	w/o locus { "C. elegans predicted micro RNA"; }
#    
# stRNA genes
#	w/locus { "C. elegans regulatory RNA $locus"; }
#	w/o locus { "C. elegans predicted regulatory RNA"; }
#    
# scRNA genes
#	w/locus { "C. elegans small cytoplasmic RNA $locus"; }
#	w/o locus { "C. elegans predicted small cytoplasmic RNA"; }
    

###########################################################################
# USAGE: get_pfam.pl <sequence (opt.)>
#  if a sequence is given, a DB_remark is generated for only that object;
#  otherwise, all of camace is checked
#
#  - connects to camace for various data
#  - connects to acedb for pfam/interpro/blastx hits
#
###########################################################################


use strict;
use Ace;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{CVS_DIR};
use Wormbase;

# proteins to ignore
#my @ignoreproteins;
#open(IGNORE, "</analysis7/analysis/db/ACEDB/ignore_proteins.ace");
#while(<IGNORE>) {
#    my $prot = $_;
#    my ($igprot) = $prot =~ /\"(\S+)\"/;
#    push(@ignoreproteins, $igprot);
#}

#print "@ignoreproteins\n";



 my $tace = &tace;
# database paths
my %db_list = (
	    'acedb' => glob("~wormpub/DATABASES/current_DB"),
	    'camace'=> glob("~wormpub/camace_ar2")
	);


my $db_name = 'camace';
my $db_path = $db_list{$db_name};
# print "trying to connect to $db_path\n";
my $db = Ace->connect (-path => "$db_path",
		       -program => $tace) || 
    die "cannot connect to camace at $db_path\n";
$db->auto_save(0);


my $db_name1 = 'acedb';
my $db_path1 = $db_list{$db_name1};
# print "trying to connect to $db_path1\n";
my $db1 = Ace->connect (-path => "$db_path1",
		       -program => $tace) || 
    die "cannot connect to acedb at $db_path1\n";
$db1->auto_save(0);


##################################################

# get subsequences

my @sequences;

if($ARGV[0]) {
    push(@sequences, $ARGV[0]);
}
else {
    @sequences = $db->fetch(-query => 'Find CDS *.* Species="Caenorhabditis elegans" & From_Laboratory=HX & Method!=history');
#### this line is for a test that gets only B* clones:
#    @sequences = $db->fetch(-query => 'Find Sequence B*.* Species="Caenorhabditis elegans" & From_Laboratory=RW & Method!=history');
}

#print "checking $#sequences sequences\n";

SUBSEQUENCE: foreach (@sequences) {
    my $gene = $_;
    chomp($gene);
#    print " ### gene: $_    acedb connection : $db1\n";

    my (@motif, @pep, $protein, $locus);

    my $full_string = "";

# get data from camace

    my $obj = $db->fetch(CDS => $gene);
    $locus = $obj->Locus(1);
    my $thismethod = $obj->Method(1);

# get pfam data from acedb

    my $obj1 = $db1->fetch(CDS => $gene);

#    print "testing: $gene $obj1\n";

    if($obj1) {
	$protein = $obj1->Corresponding_protein(1);
	my $protein_obj;
	if($protein) {
	    $protein_obj = $db1->fetch(Protein => $protein);
	}
	if($protein_obj) {
	    @motif = $protein_obj->Motif_homol;
	    unless($motif[0]) {
		@pep = $protein_obj->Pep_homol(1);

#		print "protein hits: @pep\n";

	    }
	}
    }
    else { # if changes haven't carried through to acedb yet,
           # check "parent_clone" (first check with no letter suffix,
           # then ".a" suffix)

	next SUBSEQUENCE unless($gene =~ /\D$/);  # skip if does not end in a letter; ie, skip new genes not alternatively spliced (ie, would skip B0205.15 but not skip B0205.15c)

	my $parent_gene = $gene;
	chop($parent_gene);
	my $obj1 = $db1->fetch(CDS => $parent_gene);
	my $method;
	if($obj1) {  # eg. B0205.15 exists
	    $method = $obj1->Method(1); # test for empty sequence objects
	}
	if(($obj1) && ($method)) {
	    $protein = $obj1->Corresponding_protein(1);
	    my $protein_obj = $db1->fetch(Protein => $protein);
	    if($protein_obj) {
		@motif = $protein_obj->Motif_homol(1);
		unless($motif[0]) {
		    @pep = $protein_obj->Pep_homol(1);
		}
	    }
	}
	else {   # eg. B0205.15 does not exist (or has not method); check if B0205.15a exists
	    my $parent_gene_a = $parent_gene . 'a';
	    my $obj1 = $db1->fetch(CDS => $parent_gene_a);
	    if($obj1) {
		$protein = $obj1->Corresponding_protein(1);
		my $protein_obj = $db1->fetch(Protein => $protein);
		if($protein_obj) {

		    @motif = $protein_obj->Motif_homol(1);
		    unless($motif[0]) {
			@pep = $protein_obj->Pep_homol(1);
		    }
		}
		else { # if all else fails
		    print "$gene not in acedb\n\n";
		}
	    }
	}
    }

    if(($locus) && (!($motif[0]))) {  # locus and no motif
	my $prot = uc($locus);
	$full_string .= "C. elegans $prot protein";
    }
    elsif(($locus) && ($motif[0])) { # locus and motif
	my $prot = uc($locus);
	$full_string .= "C. elegans $prot protein; ";
    }

    if($motif[0]) { # with or without locus

	my %pfamhits;
	my %interprohits;

	foreach(@motif) {
	    if($_ =~ /PFAM/) {
		my $motif_obj = $db1->fetch(Motif => $_);
		my $title = $motif_obj->Title(1);
		my ($pfam_motif) = $_ =~ /\w+\:(\w+)/;
		$pfamhits{$pfam_motif} = $title;
	    }
	    if($_ =~ /INTERPRO/) {
		my $motif_obj = $db1->fetch(Motif => $_);
		my $title = $motif_obj->Title(1);
		my ($interpro_motif) = $_ =~ /\w+\:(\w+)/;
		$interprohits{$interpro_motif} = $title;
	    }
	}

	my @pfamelements = %pfamhits;
	my @interproelements = %interprohits;

	if($#pfamelements == 1) {
	    foreach (keys %pfamhits) {
		$full_string .= "contains similarity to Pfam domain $_ ($pfamhits{$_})";
	    }
	    goto PRINTIT;
	}
	if($#pfamelements > 1) {
	    my $count = 1;
	    $full_string .= "contains similarity to Pfam domains ";
	    foreach (keys %pfamhits) {
		$full_string .= "$_ ($pfamhits{$_})";
		if ($count < $#pfamelements) {
		    $full_string .= ", ";
		}
		$count += 2;
	    }
	    goto PRINTIT;
	}

	if($#interproelements == 1) {
	    foreach (keys %interprohits) {
		$full_string .= "contains similarity to Interpro domain $_ ($interprohits{$_})";
	    }
	    goto PRINTIT;
	}
	if($#interproelements > 1) {
	    my $count = 1;
	    $full_string .= "contains similarity to Interpro domains ";
	    foreach (keys %interprohits) {
		$full_string .= "$_ ($interprohits{$_})";
		if ($count < $#interproelements) {
		    $full_string .= ", ";
		}
		$count += 2;
	    }
	    goto PRINTIT;
	}

    }

#####################################################
# no pfam or interpro hits; getting protein matches
#####################################################

    elsif($pep[0]) {

	my %pepscore;
	my %peptitle;
	my %pepspecies;

PROTEIN:	foreach(@pep) {
    my $thispep = $_;
    next PROTEIN if($thispep =~ /BP\:CBP/);
	    my ($a,$b,$c,$d,$e,$f,$g) = $_->row;
	    next PROTEIN if ($b eq 'wublastp_worm');  # don't want similarities to other worm proteins
#	    next PROTEIN if($_ eq 'BP:CBP20482');
#	    foreach(@ignoreproteins) {
#		if($thispep eq $_) {
#		    print "ignoreing $thispep for $gene\n";
#		    next PROTEIN;
#		}
#	    }
	    $pepscore{$thispep} = $c;
	    my $pep_obj = $db1->fetch(Protein => $thispep);

#	    print "getting title/desc for $thispep ; pep object : $pep_obj ; clone : $gene\n";
    next PROTEIN unless($pep_obj);
#	    my $title = $pep_obj->Title(1);
    my $title = $pep_obj->Description(1);
	    
	    $peptitle{$thispep} = $title;
	    my $species = $pep_obj->Species(1);
	    $pepspecies{$thispep} = $species;

	}

# sort score hash by value
	for (sort { $pepscore{$a} <=> $pepscore{$b} } keys %pepscore) {

	    if($peptitle{$_}) { # takes highest score WITH a title
		if($locus) {    # don't always take a peptide match, so can't add "; " above, must add here
		    $full_string .= "; ";
		}
		$full_string .= "contains similarity to $pepspecies{$_} $peptitle{$_}; $_";
		goto PRINTIT;
	    }
	}
    }

PRINTIT:

    next SUBSEQUENCE if ($full_string eq "");

    print "CDS : $gene\n";
    print "-D DB_remark\n";
    print "DB_remark \"$full_string\"\n\n";


}




# get pseudogene objects

my @pseudogenes = $db->fetch(-query => 'Find Pseudogene');

PSEUDOGENE: foreach (@pseudogenes) {
    my $gene = $_;
    chomp($gene);
#    print "gene: $_\n";

    my ($pseudogene1, $locus, $thismethod);

    my $full_string = "";

# get data from camace

    my $obj = $db->fetch(Pseudogene => $gene);
    $locus = $obj->Locus(1);
    $thismethod = $obj->Method(1);
    $pseudogene1 = $obj->Coding_pseudogene(1); # type of pseudogene

    next PSEUDOGENE if ($thismethod eq 'history');

    if($thismethod eq 'Pseudogene') {
	if(($pseudogene1) || ($locus)) {
	    if (($pseudogene1) && ($locus)) { 
		$full_string .= "C. elegans $pseudogene1 pseudogene $locus"; 
	    }
	    elsif($pseudogene1) { 
		$full_string .= "C. elegans $pseudogene1 pseudogene"; 
	    }
	    elsif($locus) {
		$full_string .= "C. elegans pseudogene $locus"; 
	    }
	}
	else {
	    $full_string .= "C. elegans predicted pseudogene"; 
	}
    }

    next PSEUDOGENE if ($full_string eq "");

    print "Pseudogene : $gene\n";
    print "-D DB_remark\n";
    print "DB_remark \"$full_string\"\n\n";


}



# get transcript objects

#my @transcripts = $db->fetch(-query => 'Find Transcript *.* Species="Caenorhabditis elegans" & From_Laboratory=RW & Method!=history');
my @transcripts = $db->fetch(-query => 'Find Transcript');

TRANSCRIPT: foreach (@transcripts) {
    my $gene = $_;
    chomp($gene);
#    print "gene: $_\n";

    my ($transcript1, $locus, $transcript2, $thismethod);

    my $full_string = "";

# get data from camace

    my $obj = $db->fetch(Transcript => $gene);
    $locus = $obj->Locus(1);
    $thismethod = $obj->Method(1);
    $transcript1 = $obj->Transcript(1); # type of transcript
    $transcript2 = $obj->Transcript(2); # text

    next TRANSCRIPT if ($thismethod eq 'history');


    if($thismethod eq 'tRNAscan-SE-1.23') { # tRNAs
	if ($locus) { $full_string .= "C. elegans tRNA $locus"; }
	else { $full_string .= "C. elegans predicted tRNA"; }
    }
    elsif($transcript1 eq 'misc_RNA') { # RNA genes
	if ($locus) { $full_string .= "C. elegans non-protein coding RNA $locus"; }
	else { $full_string .= "C. elegans probable non-coding RNA"; }
    }
    elsif($transcript1 eq 'snRNA') { # snRNA genes
	if ($locus) { $full_string .= "C. elegans small nuclear RNA $transcript2 $locus"; }
	else { $full_string .= "C. elegans small nuclear RNA $transcript2"; }
    }
    elsif($transcript1 eq 'snoRNA') { # snoRNA genes
	if ($locus) { $full_string .= "C. elegans small nucleolar RNA $transcript2 $locus"; }
	else { $full_string .= "C. elegans small nucleolar RNA $transcript2"; }
    }
    elsif($transcript1 eq 'miRNA') { # miRNA genes
	if ($locus) { $full_string .= "C. elegans microRNA $locus"; }
	else { $full_string .= "C. elegans predicted micro RNA"; }
    }
    elsif($transcript1 eq 'stRNA') { # stRNA genes
	if ($locus) { $full_string .= "C. elegans regulatory RNA $locus"; }
	else { $full_string .= "C. elegans predicted regulatory RNA"; }
    }
    elsif($transcript1 eq 'scRNA') { # scRNA genes
	if ($locus) { $full_string .= "C. elegans small cytoplasmic RNA $locus"; }
	else { $full_string .= "C. elegans predicted small cytoplasmic RNA"; }
    }

    next TRANSCRIPT if ($full_string eq "");

    print "Transcript : $gene\n";
    print "-D DB_remark\n";
    print "DB_remark \"$full_string\"\n\n";


}





