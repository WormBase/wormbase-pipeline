#!/usr/local/bin/perl5.8.0 -w
#
# get_pfam.pl
#
# by Darin Blasiar
#
# This script interogates an ACEDB database and returns all pfam/Interpro/blastx 
# data as appropriate and generates a suitable DB_remark
#
# Last updated on: $Date: 2014-01-23 13:38:37 $
# Last updated by: $Author: pad $


### DB_remark is generated as follows:  ###
#

###########################
# SEQUENCE OBJECTS        #
###########################
# Real Genes (clone.number)
# --------------------------
# 1. CGC name --> "C. elegans CGC name protein"
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
# scRNA genes
#	w/locus { "C. elegans small cytoplasmic RNA $locus"; }
#	w/o locus { "C. elegans predicted small cytoplasmic RNA"; }
#    
# lincRNA genes
#	w/locus { "C. elegans large intervening non-coding RNA RNA $locus"; }
#	w/o locus { "C. elegans predicted large intervening non-coding RNA"; }
#    
# asRNA genes
#	w/locus { "C. elegans antisense RNA $locus"; }
#	w/o locus { "C. elegans predicted antisense RNA"; }
# piRNA genes
#	w/locus { "C. elegans piwi-associated RNA $locus"; }
#	w/o locus { "C. elegans predicted piwi-associated RNA"; }


    
use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;

#################################
#
#  Command line options stuff
#
#################################

my ($help, $debug, $test, $verbose, $species, $store, $wormbase);

my $database; # which database to query
my $output;   # choose different location of output file
my $gene;     # to test on a single gene

GetOptions("help"   => \$help,
	   "debug=s"    => \$debug,
	   "test"       => \$test,
	   "verbose"    => \$verbose,
	   "species:s" 	=> \$species,
	   "database=s" => \$database,
	   "store:s"    => \$store,
	   "output=s"   => \$output,
	   "gene=s"     => \$gene
	  );



if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}


my $log = Log_files->make_build_log($wormbase);



#################################
#
#  Misc. variables
#
#################################

my $tace            = $wormbase->tace;        # TACE PATH
my $basedir         = $wormbase->basedir;
$database = $wormbase->autoace unless $database;
my $output_file;  # specify output file location
my $runtime;

$species = "elegans" unless $species;

# set up file paths
if ($output) {
  $output_file = $output;
} else {
  $output_file   = $wormbase->acefiles . "/misc_DB_remark.ace";
}
print "Output file is $output_file\n\n";
open (ACE,">$output_file") or die "Can't open output file $output_file:\t$!\n";


# open database connection

my $db = Ace->connect (-path => "$database",
		       -program => $tace) || die "cannot connect to database at $database\n";


##################################################
#
# Main loop through CDS class
#
##################################################

$runtime= $wormbase->runtime;
$log->write_to("$runtime: Processing CDS class\n");
my $object_count = 0;
my $remark_count = 0;

# get CDSs for C. elegans
my $CDSs;
if ( $gene ) {
  $CDSs = $db->fetch_many(-query => "Find CDS $gene");
} else {
  $CDSs = $db->fetch_many(-query => 'Find CDS where Method = "curated" AND Species = "'.$wormbase->full_name.'"');
}

SUBSEQUENCE: while ( my $cds = $CDSs->next ) {
  
  my (@motifs, @peptide_homols, $protein, $gene_name, $gene, $cgc_name, $cgc_protein_name);

  my $full_string = "";

  $gene = $cds->Gene;
  unless( defined $gene) {
    $log->write_to("ERROR :".$cds->name." has no Gene\n");
    next;
  }

  if (defined($gene->CGC_name)) {
    $cgc_name = $gene->CGC_name;
    $cgc_protein_name = uc($cgc_name);
  }
  
  # Find motifs
  $protein = $cds->Corresponding_protein;
  if ($protein) {
    @motifs = $protein->Motif_homol;
    if ($motifs[0]) {
      $full_string .= $wormbase->full_name('-short' => 1)." $cgc_protein_name protein; " if ($cgc_name);
      # process motif information if present
      my %pfamhits;
      my %interprohits;

      my %pfam_count;
      foreach my $motif (@motifs) {
	my $title = $motif->Title;
	if ($motif =~ /PFAM/) {
	  my ($pfam_motif) = $motif =~ /\w+\:(\w+)/;
	  $pfamhits{$pfam_motif} = $title;
	  my $pointer = $motif->right->right;
	  $pfam_count{$motif->name} = 1;
	  while ($pointer->down ) {
	    $pfam_count{$motif->name}++;
	    $pointer = $pointer->down;
	  }
	}
	if ($motif =~ /INTERPRO/) {
	  my ($interpro_motif) = $motif =~ /\w+\:(\w+)/;
	  $interprohits{$interpro_motif} = $title;
	}
	# free up memory
	$motif->DESTROY();
      }

      my @pfamelements = %pfamhits;
      my @interproelements = %interprohits;

      if ($#pfamelements == 1) {
			foreach (keys %pfamhits) {
	  			$full_string .= "contains similarity to Pfam domain $_ ";
	  			$full_string .= "$pfamhits{$_} " if $pfamhits{$_};
			}
      }
      elsif ($#pfamelements > 1) {
	my $count = 1;
	$full_string .= "contains similarity to Pfam domains ";
	foreach (keys %pfamhits) {
	  if ($pfamhits{$_}) {
	    $full_string .= "$_ ($pfamhits{$_}) ";
	    $full_string .=  "($pfam_count{\"PFAM:$_\"})" if $pfam_count{"PFAM:$_"} > 1;
	    if ($count < $#pfamelements) {
	      $full_string .= ", ";
	    }
	    $count += 2;
	  }
	}
      }

      if ($#interproelements == 1) {
	foreach (keys %interprohits) {
	  $full_string .= "contains similarity to Interpro domain $_ ($interprohits{$_}) ";
	}
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
      }

    }
    
    # get peptide homologies if no motif data
    else {
      @peptide_homols = $protein->Pep_homol;
      $full_string .= $wormbase->full_name('-short' => 1)." $cgc_protein_name protein " if ($cgc_name); 
      #####################################################
      # no pfam or interpro hits; getting protein matches
      #####################################################

      if ($peptide_homols[0]) {

	# stored details of match with highest score
	my $max_score = 0;
	my $best_match ;
	my $best_description ;
	my $best_species ;

      PROTEIN: foreach my $protein (@peptide_homols) {

	  # ignore other worm matches
	  next PROTEIN if ($protein->Corresponding_CDS);

	  my ($a,$b,$score,$d,$e,$f,$g) = $protein->row;

	  my $title = $protein->Description;
	  my $protein_species = $protein->Species;
    
	  # replace details if you find better score
	  if (($score > $max_score) && $title && $protein_species && $protein) {
	    $max_score = $score;
	    $best_match = $protein;
	    $best_description = $title;
	    $best_species = $protein_species;
	  }

	  $protein->DESTROY();
	}
    
	if ($cgc_name) {	# don't always take a peptide match, so can't add "; " above, must add here
	  $full_string .= "; ";
	}
	if ($best_species && $best_description && $best_match) {
	  $full_string .= "contains similarity to $best_species $best_description; $best_match ";
	}
      }
    }

    $protein->DESTROY();
  }
  else {
    $log->write_to("ERROR: ".$cds->name." has no Corresponding_protein\n");
  }

  $object_count++;
  next SUBSEQUENCE if ($full_string eq "");
  $remark_count++;

  print ACE "CDS : $cds\n";
  print ACE "-D DB_remark\n";
  print ACE "DB_remark \"$full_string\"\n\n";

  print "$cds\t$full_string\n" if $verbose;
}
$log->write_to("Found $object_count CDS, added remarks to $remark_count\n\n");





###########################################
#
# process pseudogene class
#
###########################################

$runtime= $wormbase->runtime;
$log->write_to("$runtime: Processing pseudogene class\n");
$object_count = 0;
$remark_count = 0;

# get Pseudogenes
my @pseudogenes;
  if ($gene) {
    @pseudogenes = $db->fetch(-query => 'Find Pseudogene $gene');
  }
else {
  @pseudogenes = $db->fetch(-query => 'Find Pseudogene where species = "'.$wormbase->full_name.'" AND (NOT Method = history_pseudogene) AND (NOT Method = Transposon_Pseudogene)');
}

PSEUDOGENE: foreach my $pseudogene (@pseudogenes) {
  
  my ($type, $gene_name, $gene, $cgc_name);

  my $full_string = "";


  # grab Gene ID, and use this to look up Gene object to get CGC_name if present
  if(!defined($pseudogene->Gene)){
    $log->write_to("ERROR: $pseudogene does not have a Gene tag.  This is bad!\n");
    next PSEUDOGENE;
  }
  else{
    $gene_name = $pseudogene->Gene;
    $gene = $db->fetch(Gene => $gene_name);
    if(defined($gene->CGC_name)){
      $cgc_name = $gene->CGC_name;
    }
  }


  # get type of pseudogene
  $type = $pseudogene->Coding_pseudogene;

  if (($type) || ($cgc_name)) {
    if (($type) && ($cgc_name)) { 
      $full_string .= $wormbase->full_name('-short' => 1)." $type pseudogene $cgc_name "; 
    } 
    elsif ($type) { 
      $full_string .= $wormbase->full_name('-short' => 1)." $type pseudogene "; 
    } 
    elsif ($cgc_name) {
      $full_string .= $wormbase->full_name('-short' => 1)." pseudogene $cgc_name "; 
    }
  } 
  else {
    $full_string .= $wormbase->full_name('-short' => 1)." predicted pseudogene "; 
  }

  $object_count++;
  next PSEUDOGENE if ($full_string eq "");
  $remark_count++;

  print ACE "Pseudogene : $pseudogene\n";
  print ACE "-D DB_remark\n";
  print ACE "DB_remark \"$full_string\"\n\n";

  # kill object to free memory
  $pseudogene->DESTROY();
}
$log->write_to("Found $object_count Pseudogenes, added remarks to $remark_count\n\n");



###########################################
#
# process transcript class
#
###########################################

$runtime= $wormbase->runtime;
$log->write_to("$runtime: Processing transcript class\n");
$object_count = 0;
$remark_count = 0;

# get Transcripts
my $transcripts;
if ($gene) {
  $transcripts = $db->fetch_many(-query => 'Find Transcript $gene');
}
else {
  $transcripts = $db->fetch_many(-query => 'Find Transcript where (NOT Method = Coding_transcript) AND (NOT Method = history_transcript) AND (Species = "'.$wormbase->full_name.'")');
}

TRANSCRIPT: while ( my $transcript = $transcripts->next ) {

  my ($description, $cgc_name, $gene_name, $gene, $type, $method);

  my $full_string = "";

  # grab Gene ID, and use this to look up Gene object to get CGC_name if present

  if(!defined($transcript->Gene)){
    $log->write_to("ERROR: $transcript does not have a Gene tag.  This is bad!\n");
    next TRANSCRIPT;
  }
  else{
    $gene_name = $transcript->Gene;
    $gene = $db->fetch(Gene => $gene_name);
    if(defined($gene->CGC_name)){
      $cgc_name = $gene->CGC_name;
    }
  }

  if($transcript->Method(1)){
    $method = $transcript->Method(1);
  }
  else{
    $log->write_to("ERROR: $transcript has no Method set\n");
    print "ERROR: $transcript has no Method set\n";
  }    

  # get type of transcript
  if($transcript->Transcript){
    $type = $transcript->Transcript; 
    ($description = $transcript->Transcript(2)) if ($transcript->Transcript(2)); # text field, not always present
  }
  else{
    $type = "";			# non-coding transcript isoforms have no tag after 'Transcript'
  }

  # set empty text field if $description is empty to prevent -w warnings
  $description = "" if (!defined($description));

  if ($type eq '') { # non-coding transcript isoforms have no tag after 'Transcript'
    if ($cgc_name) {
      $full_string .= $wormbase->full_name('-short' => 1)." non-coding isoform $cgc_name ";
    } 
    else {
      $full_string .= $wormbase->full_name('-short' => 1)." predicted non-coding isoform ";
    } 
  }   
  elsif ($method =~/tRNAscan/) { # tRNAs
    if ($cgc_name) {
      $full_string .= $wormbase->full_name('-short' => 1)." tRNA $cgc_name ";
    } 
    else {
      $full_string .= $wormbase->full_name('-short' => 1)." predicted tRNA ";
    }
  } 
  elsif ($type eq 'ncRNA') { # RNA genes
    if ($cgc_name) {
      $full_string .= $wormbase->full_name('-short' => 1)." non-protein coding RNA $cgc_name ";
    } 
    else {
      $full_string .= $wormbase->full_name('-short' => 1)." probable non-coding RNA ";
    }
  } 
  elsif ($type eq 'snRNA') { # snRNA genes
    if ($cgc_name) {
      $full_string .= $wormbase->full_name('-short' => 1)." small nuclear RNA $description $cgc_name ";
    } 
    else {
      $full_string .= $wormbase->full_name('-short' => 1)." small nuclear RNA $description ";
    }
  } 
  elsif ($type eq 'snoRNA') { # snoRNA genes
    if ($cgc_name) {
      $full_string .= $wormbase->full_name('-short' => 1)." small nucleolar RNA $description $cgc_name ";
    } 
    else {
      $full_string .= $wormbase->full_name('-short' => 1)." small nucleolar RNA $description ";
    }
  } 
  elsif ($type eq 'miRNA') { # miRNA genes
    if ($cgc_name) {
      $full_string .= $wormbase->full_name('-short' => 1)." microRNA $cgc_name ";
    } 
    else {
      $full_string .= $wormbase->full_name('-short' => 1)." predicted micro RNA ";
    }
  } 
  elsif ($type eq 'scRNA') { # scRNA genes
    if ($cgc_name) {
      $full_string .= $wormbase->full_name('-short' => 1)." small cytoplasmic RNA $cgc_name ";
    } else {
      $full_string .= $wormbase->full_name('-short' => 1)." predicted small cytoplasmic RNA ";
    }
  }
  elsif ($type eq 'lincRNA') { # lincRNA genes
    if ($cgc_name) {
      $full_string .= $wormbase->full_name('-short' => 1)." large intervening non-coding RNA $cgc_name ";
    } else {
      $full_string .= $wormbase->full_name('-short' => 1)." predicted large intervening non-coding RNA ";
    }
  }
  elsif ($type eq 'asRNA') { # asRNA genes
    if ($cgc_name) {
      $full_string .= $wormbase->full_name('-short' => 1)." antisense RNA $cgc_name ";
    } else {
      $full_string .= $wormbase->full_name('-short' => 1)." predicted antisense RNA ";
    }
  }
    elsif ($type eq 'piRNA') { # piRNA genes
    if ($cgc_name) {
      $full_string .= $wormbase->full_name('-short' => 1)." piwi-associated RNA $cgc_name ";
    } else {
      $full_string .= $wormbase->full_name('-short' => 1)." predicted piwi-associated RNA ";
    }
  }


  $object_count++;
  next TRANSCRIPT if ($full_string eq "");
  $remark_count++;

  print ACE "Transcript : $transcript\n";
  print ACE "-D DB_remark\n";
  print ACE "DB_remark \"$full_string\"\n\n";

}
$log->write_to("Found $object_count Pseudogenes, added remarks to $remark_count\n\n");

# tidy up
close ACE;
$db->close;

# load the file to autoace
unless ($gene) {
$runtime= $wormbase->runtime;
$log->write_to("$runtime: loading results to $database\n");
$wormbase->load_to_database($database, $output_file, "DB_remark", $log) unless $debug;

##################
# Check the files
##################
$wormbase->check_files($log);
}


$runtime= $wormbase->runtime;
$log->write_to("$runtime: Finished script\n\n");

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);


##############################################################
#
# Subroutines
#
##############################################################



##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################









# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - get_pfam.pl

=head1 USAGE

=over 4

=item get_pfam.pl  [-options]

=back

This script does..

get_pfam.pl MANDATORY arguments:

=over 4

=item None at present.

=back

get_pfam.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Unknown

=back

=cut
