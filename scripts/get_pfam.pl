#!/usr/local/bin/perl5.8.0 -w
#
# get_pfam.pl
#
# by Darin Blasiar
#
# This script interogates an ACEDB database and returns all pfam/Interpro/blastx 
# data as appropriate and generates the DB_remark for stlace
#
# Last updated on: $Date: 2004-06-09 12:34:14 $
# Last updated by: $Author: krb $




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
use Getopt::Long;
use Log_files;

my ($build,$test,$debug);
my ($camace,$stlace);
my $CDS_source_db;
my $remark_target;

GetOptions(
	   "build"    => \$build,
	   "debug=s"  => \$debug,
	   "test"     => \$test,
	   "source:s" => \$CDS_source_db,
	  );


my $tace = &tace;
my $log;
my $basedir;

if( $test ) {
  $basedir = glob("~wormpub/TEST_BUILD");
}
else {
  $basedir = "/wormsrv2";
}

if ( $build ) {
  $CDS_source_db = "$basedir/autoace";
  $remark_target = "$basedir/autoace";
  $file          = "$basedir/wormbase/misc/misc_DB_remark.ace";
  $log           = Log_files->make_build_log() unless $test;
}
else {
  $log = Log_files->make_log("/tmp/$0.$$") unless $test;
}

open (ACE,">$file") or die "cant open output file $file:\t$!\n";

# database connections(s)
# print "trying to connect to $db_path\n";
my $db = Ace->connect (-path => "$CDS_source_db",
		       -program => $tace) || 
  die "cannot connect to database at $CDS_source_db\n";
$db->auto_save(0);

my $db1;
if ( $CDS_source_db eq $remark_target ) {
  $db1 = $db;
} 
else {
  $db1 = Ace->connect (-path => "$remark_target",
		       -program => $tace) || 
			 die "cannot connect to acedb at $remark_target\n";
  $db1->auto_save(0);
}

##################################################

# get subsequences

my @sequences;

if ($ARGV[0]) {
  push(@sequences, $ARGV[0]);
} else {
  @sequences = $db->fetch(-query => 'Find elegans_CDS');
}

#print "checking $#sequences sequences\n";

SUBSEQUENCE: foreach my $cds (@sequences) {
  chomp($cds);

  my (@motif, @pep, $protein, $gene_name, $gene, $cgc_name);

  my $full_string = "";

  # get data from camace


  my $obj = $db->fetch(CDS => $cds);

  # grab Gene ID, and use this to look up Gene object to get CGC_name if present
  $gene_name = $obj->Gene;
  $gene = $db->fetch(Gene => $gene_name);
  if(defined($gene->CGC_name)){
    $cgc_name = $gene->CGC_name;
  }


  my $thismethod = $obj->Method(1);

  # get pfam data from acedb

  my $obj1 = $db1->fetch(CDS => $cds);

  #    print "testing: $cds $obj1\n";

  if ($obj1) {
    $protein = $obj1->Corresponding_protein(1);
    if ($protein) {
      my $protein_obj;
      if ($protein_obj) {
	$protein_obj = $db1->fetch(Protein => $protein);
      }

      (@motif = $protein_obj->Motif_homol) if ($protein_obj->Motif_homol);
      unless($motif[0]) {
	@pep = $protein_obj->Pep_homol(1);
      }
    }
  } else {			# if changes haven't carried through to acedb yet,
				# check "parent_clone" (first check with no letter suffix,
				# then ".a" suffix)

    next SUBSEQUENCE unless($cds =~ /\D$/); # skip if does not end in a letter; ie, skip new genes not alternatively spliced (ie, would skip B0205.15 but not skip B0205.15c)

    my $parent_gene = $cds;
    chop($parent_gene);
    my $obj1 = $db1->fetch(CDS => $parent_gene);
    my $method;
    if ($obj1) {		# eg. B0205.15 exists
      $method = $obj1->Method(1); # test for empty sequence objects
    }
    if (($obj1) && ($method)) {
      $protein = $obj1->Corresponding_protein(1);
      my $protein_obj = $db1->fetch(Protein => $protein);
      if ($protein_obj) {
	@motif = $protein_obj->Motif_homol(1);
	unless($motif[0]) {
	  @pep = $protein_obj->Pep_homol(1);
	}
      }
    } else {			# eg. B0205.15 does not exist (or has not method); check if B0205.15a exists
      my $parent_gene_a = $parent_gene . 'a';
      my $obj1 = $db1->fetch(CDS => $parent_gene_a);
      if ($obj1) {
	$protein = $obj1->Corresponding_protein(1);
	my $protein_obj = $db1->fetch(Protein => $protein);
	if ($protein_obj) {

	  @motif = $protein_obj->Motif_homol(1);
	  unless($motif[0]) {
	    @pep = $protein_obj->Pep_homol(1);
	  }
	} else {		# if all else fails
	  $log->write_to("$cds not in acedb\n\n");
	}
      }
    }
  }

  if (($cgc_name) && (!($motif[0]))) { # locus and no motif
    my $prot = uc($cgc_name);
    $full_string .= "C. elegans $prot protein";
  } elsif (($cgc_name) && ($motif[0])) { # locus and motif
    my $prot = uc($cgc_name);
    $full_string .= "C. elegans $prot protein; ";
  }

  if ($motif[0]) {		# with or without locus

    my %pfamhits;
    my %interprohits;

    my %pfam_count;
    foreach (@motif) {
      if ($_ =~ /PFAM/) {
	my $motif_obj = $db1->fetch(Motif => $_);
	my $title = $motif_obj->Title(1);
	my ($pfam_motif) = $_ =~ /\w+\:(\w+)/;
	$pfamhits{$pfam_motif} = $title;

	my $pointer = $_->right->right;
	$pfam_count{$_->name} = 1;
	while( $pointer->down ) {
	  $pfam_count{$_->name}++;
	  $pointer = $pointer->down;
	}
      }
      if ($_ =~ /INTERPRO/) {
	my $motif_obj = $db1->fetch(Motif => $_);
	my $title = $motif_obj->Title(1);
	my ($interpro_motif) = $_ =~ /\w+\:(\w+)/;
	$interprohits{$interpro_motif} = $title;
      }
    }

    my @pfamelements = %pfamhits;
    my @interproelements = %interprohits;

    if ($#pfamelements == 1) {
      foreach (keys %pfamhits) {
	$full_string .= "contains similarity to Pfam domain $_ ($pfamhits{$_})";
      }
      goto PRINTIT;
    }
    if ($#pfamelements > 1) {
      my $count = 1;
      $full_string .= "contains similarity to Pfam domains ";
      foreach (keys %pfamhits) {
	$full_string .= "$_ ($pfamhits{$_})";
	$full_string .=  "($pfam_count{\"PFAM:$_\"})" if $pfam_count{"PFAM:$_"} > 1;
	if ($count < $#pfamelements) {
	  $full_string .= ", ";
	}
	$count += 2;
      }
      goto PRINTIT;
    }

    if ($#interproelements == 1) {
      foreach (keys %interprohits) {
	$full_string .= "contains similarity to Interpro domain $_ ($interprohits{$_})";
      }
      goto PRINTIT;
    }
    if ($#interproelements > 1) {
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

  elsif ($pep[0]) {

    my %pepscore;
    my %peptitle;
    my %pepspecies;

  PROTEIN:	foreach(@pep) {
      my $thispep = $_;
      next PROTEIN if($thispep =~ /BP\:CBP/);
      my ($a,$b,$c,$d,$e,$f,$g) = $_->row;
      next PROTEIN if (!(defined $b) or $b eq 'wublastp_worm'); # don't want similarities to other worm proteins
      #	    next PROTEIN if($_ eq 'BP:CBP20482');
      #	    foreach(@ignoreproteins) {
      #		if($thispep eq $_) {
      #		    print "ignoreing $thispep for $cds\n";
      #		    next PROTEIN;
      #		}
      #	    }
      $pepscore{$thispep} = $c;
      my $pep_obj = $db1->fetch(Protein => $thispep);

      #	    print "getting title/desc for $thispep ; pep object : $pep_obj ; clone : $cds\n";
      next PROTEIN unless($pep_obj);
      #	    my $title = $pep_obj->Title(1);
      my $title = $pep_obj->Description(1);
	    
      $peptitle{$thispep} = $title;
      my $species = $pep_obj->Species(1);
      $pepspecies{$thispep} = $species;

    }

    # sort score hash by value
    for (sort { $pepscore{$a} <=> $pepscore{$b} } keys %pepscore) {

      if ($peptitle{$_}) {	# takes highest score WITH a title
	if ($cgc_name) {		# don't always take a peptide match, so can't add "; " above, must add here
	  $full_string .= "; ";
	}
	$full_string .= "contains similarity to $pepspecies{$_} $peptitle{$_}; $_";
	goto PRINTIT;
      }
    }
  }

 PRINTIT:

  next SUBSEQUENCE if ($full_string eq "");

  print ACE "CDS : $cds\n";
  print ACE "-D DB_remark\n";
  print ACE "DB_remark \"$full_string\"\n\n";


}




# get pseudogene objects

my @pseudogenes = $db->fetch(-query => 'Find Pseudogene & Method!=history');

PSEUDOGENE: foreach my $pseudogene (@pseudogenes) {
  chomp($pseudogene);

  my ($pseudogene1, $gene_name, $gene, $cgc_name, $thismethod);

  my $full_string = "";

  # get data from camace
  my $obj = $db->fetch(Pseudogene => $pseudogene);

  # grab Gene ID, and use this to look up Gene object to get CGC_name if present
  $gene_name = $obj->Gene;
  $gene = $db->fetch(Gene => $gene_name);
  if(defined($gene->CGC_name)){
    $cgc_name = $gene->CGC_name;
  }

  $thismethod = $obj->Method(1);
  $pseudogene1 = $obj->Coding_pseudogene(1); # type of pseudogene

  next PSEUDOGENE if ($thismethod eq 'history');

  if ($thismethod eq 'Pseudogene') {
    if (($pseudogene1) || ($cgc_name)) {
      if (($pseudogene1) && ($cgc_name)) { 
	$full_string .= "C. elegans $pseudogene1 pseudogene $cgc_name"; 
      } elsif ($pseudogene1) { 
	$full_string .= "C. elegans $pseudogene1 pseudogene"; 
      } elsif ($cgc_name) {
	$full_string .= "C. elegans pseudogene $cgc_name"; 
      }
    } else {
      $full_string .= "C. elegans predicted pseudogene"; 
    }
  }

  next PSEUDOGENE if ($full_string eq "");

  print ACE "Pseudogene : $pseudogene\n";
  print ACE "-D DB_remark\n";
  print ACE "DB_remark \"$full_string\"\n\n";


}



# get transcript objects

my @transcripts = $db->fetch(-query => 'Find Transcript *.* Species="Caenorhabditis elegans" & Method!=history');
#my @transcripts = $db->fetch(-query => 'Find Transcript');

TRANSCRIPT: foreach my $transcript (@transcripts) {
  chomp($transcript);

  my ($transcript1, $cgc_name, $gene_name, $gene, $transcript2, $thismethod);

  my $full_string = "";

  # get data from camace
  my $obj = $db->fetch(Transcript => $transcript);

  # grab Gene ID, and use this to look up Gene object to get CGC_name if present
  $gene_name = $obj->Gene;
  $gene = $db->fetch(Gene => $gene_name);
  if(defined($gene->CGC_name)){
    $cgc_name = $gene->CGC_name;
  }

  $thismethod = $obj->Method(1);
  $transcript1 = $obj->Transcript(1); # type of transcript
  $transcript2 = $obj->Transcript(2); # text

  next TRANSCRIPT if ($thismethod eq 'history' or $thismethod eq "Coding_transcript");


  if ($thismethod eq 'tRNAscan-SE-1.23') { # tRNAs
    if ($cgc_name) {
      $full_string .= "C. elegans tRNA $cgc_name";
    } else {
      $full_string .= "C. elegans predicted tRNA";
    }
  } elsif ($transcript1 eq 'misc_RNA') { # RNA genes
    if ($cgc_name) {
      $full_string .= "C. elegans non-protein coding RNA $cgc_name";
    } else {
      $full_string .= "C. elegans probable non-coding RNA";
    }
  } elsif ($transcript1 eq 'snRNA') { # snRNA genes
    if ($cgc_name) {
      $full_string .= "C. elegans small nuclear RNA $transcript2 $cgc_name";
    } else {
      $full_string .= "C. elegans small nuclear RNA $transcript2";
    }
  } elsif ($transcript1 eq 'snoRNA') { # snoRNA genes
    if ($cgc_name) {
      $full_string .= "C. elegans small nucleolar RNA $transcript2 $cgc_name";
    } else {
      $full_string .= "C. elegans small nucleolar RNA $transcript2";
    }
  } elsif ($transcript1 eq 'miRNA') { # miRNA genes
    if ($cgc_name) {
      $full_string .= "C. elegans microRNA $cgc_name";
    } else {
      $full_string .= "C. elegans predicted micro RNA";
    }
  } elsif ($transcript1 eq 'stRNA') { # stRNA genes
    if ($cgc_name) {
      $full_string .= "C. elegans regulatory RNA $cgc_name";
    } else {
      $full_string .= "C. elegans predicted regulatory RNA";
    }
  } elsif ($transcript1 eq 'scRNA') { # scRNA genes
    if ($cgc_name) {
      $full_string .= "C. elegans small cytoplasmic RNA $cgc_name";
    } else {
      $full_string .= "C. elegans predicted small cytoplasmic RNA";
    }
  }

  next TRANSCRIPT if ($full_string eq "");

  print ACE "Transcript : $transcript\n";
  print ACE "-D DB_remark\n";
  print ACE "DB_remark \"$full_string\"\n\n";


}

close ACE;

$db->close;
$db1->close;


# load the file
&load_to_database($remark_target,$file,"DB_remark");

$log->mail unless $test;

exit(0);






