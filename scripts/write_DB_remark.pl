#!/usr/bin/env perl
#
# write_DB_rematk.pl
#
#
# This script interogates an ACEDB database and returns all pfam/Interpro/blastx 
# data as appropriate and generates a suitable DB_remark
#
# Last updated on: $Date: 2015-06-01 11:49:41 $
# Last updated by: $Author: klh $

###########################
# SEQUENCE OBJECTS        #
###########################
# Real Genes (clone.number)
# --------------------------
# 1. CGC name --> "CGC name protein"
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
#
# piRNA genes
#	w/locus { "C. elegans piwi-associated RNA $locus"; }
#	w/o locus { "C. elegans predicted piwi-associated RNA"; }
# circRNA genes
#	w/locus { "C. elegans circular RNA $locus"; }
#	w/o locus { "C. elegans predicted circular RNA"; }


    
use strict;                                      
use Storable;
use Getopt::Long;
use Carp;
use Ace;

use lib $ENV{CVS_DIR};
use Wormbase;
use Log_files;

my ($help, $debug, $test, $verbose, $store, $wormbase);

my ($database, $output_file, $gene, $no_load, $do_cds, $do_transcript, $do_pseudogene, $runtime);

GetOptions("debug=s"        => \$debug,
	   "test"           => \$test,
	   "verbose"        => \$verbose,
	   "database=s"     => \$database,
	   "store:s"        => \$store,
	   "acefile=s"      => \$output_file,
	   "gene=s"         => \$gene,
           "noload"         => \$no_load,
           "cds"            => \$do_cds,
           "transcript"     => \$do_transcript,
           "pseudogene"     => \$do_pseudogene,
	  );



if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

my $log = Log_files->make_build_log($wormbase);


my $tace            = $wormbase->tace;
my $basedir         = $wormbase->basedir;

$database = $wormbase->autoace unless $database;
$output_file = $wormbase->acefiles . "/misc_DB_remark.ace" unless $output_file;

open (my $ace_fh,">$output_file") or die "Can't open output file $output_file:\t$!\n";

my ($object_count, $remark_count) = (0, 0);

if (not $do_cds and not $do_transcript and not $do_pseudogene) {
  $do_cds = $do_transcript = $do_pseudogene = 1;
}

if ($do_cds) {

  $verbose and $log->write_to("Connecting to database...\n");
  my $db = Ace->connect (-path => "$database",
                         -program => $tace) || $log->log_and_die("Cannot connect to database at $database\n");

  $runtime= $wormbase->runtime;
  $log->write_to("$runtime: Processing CDS class\n");

  my $CDSs = ($gene) 
      ? $db->fetch_many(-query => "Find CDS WHERE Gene = \"$gene\"")
      : $db->fetch_many(-query => 'Find CDS WHERE Method = "curated"');
  
  while (my $cds = $CDSs->next ) {
    my (@motifs, @peptide_homols, $protein, $gene_name, $gene, $cgc_name, $cgc_protein_name);
    
    my @full_string;
    
    $gene = $cds->Gene;
    unless( defined $gene) {
      $log->write_to("ERROR: ".$cds->name." has no Gene\n");
      next;
    }
    
    
    if (defined($gene->CGC_name)) {
      $cgc_name = $gene->CGC_name;
      $cgc_protein_name = uc($cgc_name);
    }
    
    push @full_string, "$cgc_protein_name protein" if ($cgc_name); 
    
    # Find motifs
    $protein = $cds->Corresponding_protein;
    if ($protein) {
      @motifs = $protein->Motif_homol;
      if ($motifs[0]) {
        my (%pfamhits, %interprohits, %pfamcount);
        
        foreach my $motif (@motifs) {
          my $title = $motif->Title;
          if ($motif =~ /PFAM/) {
            my ($pfam_motif) = $motif =~ /\w+\:(\w+)/;
            $pfamhits{$pfam_motif} = $title;
            my $pointer = $motif->right->right;
            $pfamcount{$pfam_motif} = 1;
            while ($pointer->down ) {
              $pfamcount{$pfam_motif}++;
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
        
        my (@string_els);
        foreach my $pfam_acc (sort keys %pfamhits) {
          my $pfam_name = $pfamhits{$pfam_acc};
          my $this_string .= "$pfam_acc ($pfam_name"; 
          if ($pfamcount{$pfam_acc} > 1) {
            my $count = $pfamcount{$pfam_acc};
            $this_string .=  " ($count copies)";
          }
          $this_string .= ")";
          push @string_els, $this_string;
        }
        if (@string_els) {
          my $str = (scalar(@string_els) == 1) 
              ? " @string_els" 
              : "s " . join(", ", @string_els);
          push @full_string, "contains similarity to Pfam domain${str}";            
        }
        
        @string_els = ();
        foreach my $interpro_acc (keys %interprohits) {
          my $interpro_name = $interprohits{$interpro_acc};
          push @string_els, "$interpro_acc ($interpro_name)";
        }
        if (@string_els) {
          my $str = (scalar(@string_els) == 1) 
              ? " @string_els" 
              : "s " . join(", ", @string_els);
          push @full_string, "contains similarity to Interpro domain${str}";
        }
      } else {
        #####################################################
        # no pfam or interpro hits; getting protein matches
        #####################################################
        @peptide_homols = $protein->Pep_homol;
        
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
          
          if ($best_species && $best_description && $best_match) {
            push @full_string, "contains similarity to $best_species $best_description; $best_match";
          }
        }
      }
      
      $protein->DESTROY();
    }
    else {
      $log->write_to("ERROR: ".$cds->name." has no Corresponding_protein\n");
    }
    
    $object_count++;
    
    print $ace_fh "\nCDS : $cds\n";
    print $ace_fh "-D DB_remark\n";
    
    next if not @full_string;
    my $full_string = join("; ", @full_string); 

    $remark_count++;
    
    print $ace_fh "\nCDS : $cds\n";
    print $ace_fh "DB_remark \"$full_string\"\n";
    
    $cds->DESTROY();
  }
  $log->write_to("Found $object_count CDS, added remarks to $remark_count\n\n");
  $db->close();
}


###########################################
#
# process pseudogene class
#
###########################################
if ($do_pseudogene) {
  $verbose and $log->write_to("(Re)connecting to database...\n");
  my $db = Ace->connect (-path => "$database",
                         -program => $tace) || $log->log_and_die("Cannot connect to database at $database\n");
  
  $runtime= $wormbase->runtime;
  $log->write_to("$runtime: Processing pseudogene class\n");
  ($object_count, $remark_count) = (0,0);
  
  my @pseudogenes = ($gene)
      ? $db->fetch(-query => 'Find Pseudogene Gene = \"$gene\"')
      : $db->fetch(-query => 'Find Pseudogene WHERE (NOT Method = history_pseudogene) AND (NOT Method = Transposon_Pseudogene)');
  
  foreach my $pseudogene (@pseudogenes) {
    
    my ($type, $gene_name, $gene, $cgc_name);
    
    my $full_string;
    
    # grab Gene ID, and use this to look up Gene object to get CGC_name if present
    if(!defined($pseudogene->Gene)){
      $log->write_to("ERROR: $pseudogene does not have a Gene tag.  This is bad!\n");
      next;
    } 
    
    $gene_name = $pseudogene->Gene;
    $gene = $db->fetch(Gene => $gene_name);
    if(defined($gene->CGC_name)){
      $cgc_name = $gene->CGC_name;
    }
    ($type) = $pseudogene->Type;
    
    if ($type eq 'Coding_pseudogene') {
      $full_string = "coding pseudogene";
    } elsif ($type eq 'RNA_pseudogene') {
      $full_string = "ncRNA pseudogene";
    } else {
      $full_string = "pseudogene";
    }
    
    if ($cgc_name) {
      $full_string .= " $cgc_name";
    }
    
    $object_count++;
    
    print $ace_fh "\nPseudogene : $pseudogene\n";
    print $ace_fh "-D DB_remark\n";
    
    next if not $full_string;
    $remark_count++;
    
    print $ace_fh "\nPseudogene : $pseudogene\n";
    print $ace_fh "DB_remark \"$full_string\"\n";
    
    # kill object to free memory
    $pseudogene->DESTROY();
  }
  $log->write_to("Found $object_count Pseudogenes, added remarks to $remark_count\n\n");
  $db->close();
}


###########################################
#
# process transcript class
#
###########################################
if ($do_transcript) {
  $verbose and $log->write_to("(Re)connecting to database...\n");
  
  my $db = Ace->connect (-path => "$database",
                         -program => $tace) || $log->log_and_die("Cannot connect to database at $database\n");
  
  $runtime= $wormbase->runtime;
  $log->write_to("$runtime: Processing transcript class\n");
  ($object_count, $remark_count) = (0,0);
  
  my $transcripts = ($gene)
      ? $db->fetch_many(-query => "Find Transcript WHERE Gene = \"$gene\"")
      : $db->fetch_many(-query => "Find Transcript WHERE (Method != Coding_transcript) AND (Method != history_transcript)");
  
  while ( my $transcript = $transcripts->next ) {
    my ($description, $cgc_name, $gene_name, $gene, $type, $method);
    
    my $full_string = "";
    
    # grab Gene ID, and use this to look up Gene object to get CGC_name if present
    if(!defined $transcript->Gene){ # Brugia has a lot of 'TIGR_BEST' transcripts with no Gene tag - ignore these, else throw an error.
      $log->write_to("ERROR: $transcript does not have a Gene tag.  This is bad!\n") unless ($transcript->Method eq 'TIGR_BEST');
      next; 
    } 
    
    $gene_name = $transcript->Gene;
    $gene = $db->fetch(Gene => $gene_name);
    if(defined($gene->CGC_name)){
      $cgc_name = $gene->CGC_name;
    }
    
    if($transcript->Method(1)){
      $method = $transcript->Method(1);
    }
    else{
      $log->write_to("ERROR: $transcript has no Method set\n");
    }    
    
    # get type of transcript
    if($transcript->Transcript){
      $type = $transcript->Transcript; 
      ($description = $transcript->Transcript(2)) if ($transcript->Transcript(2));
    }
    else{
      # non-coding transcript isoforms have no tag after 'Transcript'
      $type = "";			
    }
    
    # set empty text field if $description is empty to prevent -w warnings
    $description = "" if (!defined($description));
    
    if (not $type) { # non-coding transcript isoforms have no tag after 'Transcript'
      if ($cgc_name) {
        $full_string = "non-coding isoform of $cgc_name";
      } 
      else {
        $full_string = "predicted non-coding isoform";
      } 
    } elsif ($method =~/tRNAscan/) { # tRNAs
      if ($cgc_name) {
        $full_string = "tRNA $cgc_name";
      } else {
        $full_string = "predicted tRNA";
      }
    } elsif ($type eq 'ncRNA') { # RNA genes
      if ($cgc_name) {
        $full_string = "non-protein coding RNA $cgc_name";
      } else {
        $full_string = "probable non-coding RNA";
      }
    } elsif ($type eq 'snRNA') { # snRNA genes
      if ($cgc_name) {
        $full_string = "small nuclear RNA $description $cgc_name";
      } else {
        $full_string = "small nuclear RNA $description";
      }
    } elsif ($type eq 'snoRNA') { # snoRNA genes
      if ($cgc_name) {
        $full_string = "small nucleolar RNA $description $cgc_name";
      } else {
        $full_string = "small nucleolar RNA $description";
      }
    } elsif ($type eq 'miRNA') { # miRNA genes
      if ($method eq 'miRNA_primary_transcript') {
        $type = "miRNA primary transcript";
      } else {
        $type = $method;
      }
      if ($cgc_name) {
        $full_string = "$type $cgc_name";
      } else {
        $full_string = "predicted $type";
      }
    } elsif ($type eq 'scRNA') { # scRNA genes
      if ($cgc_name) {
        $full_string = "small cytoplasmic RNA $cgc_name";
      } else {
        $full_string = "predicted small cytoplasmic RNA";
      }
    } elsif ($type eq 'lincRNA') { # lincRNA genes
      if ($cgc_name) {
        $full_string = "long intervening non-coding RNA $cgc_name";
      } else {
        $full_string = "predicted long intervening non-coding RNA";
      }
    } elsif ($type eq 'asRNA') { # asRNA genes
      if ($cgc_name) {
        $full_string = "antisense RNA $cgc_name";
      } else {
        $full_string = "predicted antisense RNA";
      }
    } elsif ($type eq 'piRNA') { # piRNA genes
      if ($cgc_name) {
        $full_string = "piwi-associated RNA $cgc_name";
      } else {
        $full_string = "predicted piwi-associated RNA";
      }
    }elsif ($type eq 'circRNA') { # circRNA genes
      if ($cgc_name) {
        $full_string = "circular RNA $cgc_name";
      } else {
        $full_string = "predicted circular RNA";
      }
    }
    
    
    $object_count++;
    
    print $ace_fh "\nTranscript : $transcript\n";
    print $ace_fh "-D DB_remark\n";
    
    next if not $full_string;
    
    $remark_count++;
    
    print $ace_fh "\nTranscript : $transcript\n";
    print $ace_fh "DB_remark \"$full_string\"\n";
    
  }
  $log->write_to("Found $object_count Transcripts, added remarks to $remark_count\n\n");
  $db->close;
}

close($ace_fh);

# load the file to autoace
if (not $no_load) {
  $runtime= $wormbase->runtime;
  $log->write_to("$runtime: loading results to $database\n");
  $wormbase->load_to_database($database, $output_file, "DB_remark", $log) unless $debug;

  $wormbase->check_files($log);
}

$runtime= $wormbase->runtime;
$log->write_to("$runtime: Finished script\n\n");

$log->mail();
exit(0);

__END__

