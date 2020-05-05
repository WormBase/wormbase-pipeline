#!/usr/bin/env perl
#
# release_letter.pl                            
# 
# Last updated by: $Author: klh $               
# Last updated on: $Date: 2015-06-02 11:06:13 $

# Generates a release letter at the end of build.
#
# Three subroutines are called during the build - 
#  release_wormpep by make_wormpep
#  release_composition by dump_chromosomes.pl
#  release_databases by dbcomp
#
# These write to a file in autoace/RELEASE_LETTER and are incorperated in to the letter at the end. 
# This allows for overwriting during the build when errors are fixed and scripts rerun


use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use JSON;

use Log_files;
use Storable;
use Ace;
use Modules::Remap_Sequence_Change;
use Species;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $store, $wormbase, $database);
my ($release_letter);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "store:s"    => \$store,
            "database=s" => \$database,
            "letter=s"   => \$release_letter,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

my $log = Log_files->make_build_log($wormbase);


##############
# variables  #                                                                   
##############

my $reports_dir = $wormbase->reports;     # AUTOACE REPORTS
my $tace        = $wormbase->tace;
my $ver         = $wormbase->get_wormbase_version;
my $old_ver     = $ver -1;
$database       = $wormbase->autoace if not defined $database;

my $gene_stats_file      = "$reports_dir/genedata";
my $wormpep_stats_file   = "$reports_dir/wormpep";
my $assemblies_file      = "$reports_dir/ASSEMBLIES.json";

$log->log_and_die("Expecting to find $assemblies_file; looks like you never ran make_assembly_manifest.pl") if not -e $assemblies_file;
$log->log_and_die("Expecting to find $gene_stats_file; looks like you never ran genestats.pl") if not -e $gene_stats_file;
$log->log_and_die("Expecting to find $wormpep_stats_file; looks like you never build wormpep") if not -e $wormpep_stats_file;

$log->write_to("Connecting to database...\n");
my $db = Ace->connect(-path    => $database,
                      -program => $tace) || $log->log_and_die("Connection failure: ",Ace->error);

my (%cds_status, 
    %operon_data,
    %gene_counts,
    %cgc_counts,
    %go_stats,
    );

my %accessors = ($wormbase->species_accessors);
$accessors{$wormbase->species} = $wormbase;
 
######################################
# Read in pre-queries assembly data
######################################
my $assembly_stats = &read_json($assemblies_file);

  
######################################
#  Generate Operon stats 
######################################
$log->write_to("Retrieving Operon Data...\n");
my $operon_query = "Find Operon WHERE method = \"Operon\" AND Species = \"*elegans\"";
my $gene_query = "Find Gene where Contained_in_operon AND Species = \"*elegans\"";
$operon_data{operons}          = $db->fetch(-query=> "$operon_query");
$operon_data{genes_in_operons} = $db->fetch(-query=> "$gene_query");


######################################
# Generate wormpep stats for each species
######################################
foreach my $spe (keys %accessors) {
  $log->write_to("Generating $spe wormpep overview stats...\n");
  $cds_status{$spe} = &get_cds_confirmation_status($accessors{$spe}, $db);
}   


######################################
#  Generate gene count stats
######################################
foreach my $spec (keys %accessors) {
  $log->write_to("Retrieving basic coding/non-coding gene stats for $spec...\n");
  my $sp_name = $accessors{$spec}->full_name();
  my $gene_count_query = "Find Gene where Species = \"$sp_name\" AND Live";
  my $gene_count = $db->fetch(-query=> "$gene_count_query");
  # now count coding genes, rather than coding isoforms
  my $Coding_count_query = "Find Gene where Species = \"$sp_name\" AND Corresponding_CDS";
  my $Coding_count = $db->fetch(-query=> "$Coding_count_query");
  $gene_counts{$spec}->{genes} = $gene_count;
  $gene_counts{$spec}->{coding_genes} = $Coding_count;
}
  
######################################
# Generate "useful" stats
######################################
foreach my $spec (keys %accessors) {
  next if $spec eq 'elegans';
  $log->write_to("Retrieving CGC name stats for  $spec...\n");
  my $sp_name = $accessors{$spec}->full_name();
  my $query = "Find Gene WHERE (Corresponding_CDS OR Corresponding_transcript OR Corresponding_pseudogene) AND CGC_name AND Species = \"$sp_name\""; 
  my $cgc_count = $db->fetch(-query => $query);
  $cgc_counts{$spec} = $cgc_count;
}

######################################
# Generate GO stats - elegans only
######################################
$go_stats{elegans} = &get_GO_stats($accessors{elegans});

$db->close;

#######################################
# Finally, write it all out
#######################################

$release_letter = "$reports_dir/letter.WS$ver" if not defined $release_letter;

open (my $rlfh ,">$release_letter");
printf $rlfh "New release of WormBase WS$ver\n\n";
printf $rlfh "WS$ver was built by [INSERT NAME HERE]\n\n";
printf $rlfh "-==============================================================================-\n";
printf $rlfh "-========= FTP site structure =================================================-\n";
printf $rlfh "-==============================================================================-\n";
printf $rlfh "The WS$ver build directory includes:\n";
printf $rlfh "species/ DIR              -  contains a sub dir for each WormBase species (G_SPECIES)\n";
printf $rlfh "species/G_SPECIES DIR     -  contains a sub dir for each NCBI genome sequencing BioProject (BIOPROJECT) for the species, with the following files:\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.genomic.fa.gz                  - Unmasked genomic DNA\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.genomic_masked.fa.gz           - Hard-masked (repeats replaced with Ns) genomic DNA\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.genomic_softmasked.fa.gz       - Soft-masked (repeats lower-cased) genomic DNA\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.protein.fa.gz                  - Current live protein set\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.CDS_transcripts.fa.gz          - Spliced cDNA sequence for the CDS portion of protein-coding transcripts\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.mRNA_transcripts.fa.gz         - Spliced cDNA sequence for the full-length (including UTRs) mRNA for transcripts\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.ncrna_transcripts.fa.gz        - Spliced cDNA sequence for non-coding RNA transcripts\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.pseudogenic_transcripts.fa.gz  - Spliced cDNA sequence for pseudogenic transcripts\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.transposon_transcripts.fa.gz   - Spliced cDNA sequence for mRNAs and pseudogenes located in Transposons\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.transposons.fa.gz              - DNA sequence of curated and predicted Transposons\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.intergenic_sequences.fa.gz     - DNA sequence between pairs of adjacent genes\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.annotations.gff[2|3].gz        - Sequence features in either GFF2 or GFF3 format\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.protein_annotation.gff3.gz     - Sequence features in proteins in GFF3 format\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.canonical_geneset.gtf.gz       - Genes, transcripts and CDSs in GTF (GFF2) format\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.ests.fa.gz                     - ESTs and mRNA sequences extracted from the public databases\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.best_blastp_hits.txt.gz        - Best blastp matches to human, fly, yeast, and non-WormBase Uniprot proteins\n";
printf $rlfh "     - G_SPECIES.BIOPROJECT.WS$ver.*pep_package.tar.gz            - latest version of the [worm|brig|bren|rema|jap|ppa|brug]pep package (if updated since last release)\n";
printf $rlfh "     - annotation/                    - contains additional annotations:\n";
printf $rlfh "        - G_SPECIES.BIOPROJECT.WS$ver.confirmed_genes.txt.gz              - DNA sequences of all genes confirmed by EST &/or cDNA\n";
printf $rlfh "        - G_SPECIES.BIOPROJECT.WS$ver.cDNA2orf.txt.gz                     - Latest set of ORF connections to each cDNA (EST, OST, mRNA)\n";      
printf $rlfh "        - G_SPECIES.BIOPROJECT.WS$ver.geneIDs.txtgz                       - list of all current gene identifiers with CGC & molecular names (when known)\n";
printf $rlfh "        - G_SPECIES.BIOPROJECT.WS$ver.PCR_product2gene.txt.gz             - Mappings between PCR products and overlapping Genes\n";
printf $rlfh "        - G_SPECIES.BIOPROJECT.WS$ver.*oligo_mapping.txt.gz               - Oligo array mapping files\n";
printf $rlfh "        - G_SPECIES.BIOPROJECT.WS$ver.knockout_consortium_alleles.xml.gz  - Table of Knockout Consortium alleles\n";
printf $rlfh "        - G_SPECIES.BIOPROJECT.WS$ver.SRA_gene_expression.tar.gz          - Tables of gene expression values computed from SRA RNASeq data\n";
printf $rlfh "        - G_SPECIES.BIOPROJECT.WS$ver.TSS.wig.tar.gz                      - Wiggle plot files of Transcription Start Sites from the papers WBPaper00042246, WBPaper00042529, WBPaper00042354\n";
printf $rlfh "        - G_SPECIES.BIOPROJECT.WS$ver.repeats.fa..gz                      - Latest version of the repeat library for the genome, suitable for use with RepeatMasker\n";

printf $rlfh "acedb DIR                -  Everything needed to generate a local copy of the The Primary database\n";
printf $rlfh "     - database.WS$ver.*.tar.gz   - compressed acedb database for new release\n";
printf $rlfh "     - models.wrm.WS$ver          - the latest database schema (also in above database files)\n";
printf $rlfh "     - WS$ver-WS$old_ver.dbcomp        - log file reporting difference from last release\n";
printf $rlfh "     - *Non_C_elegans_BLASTX/     - This directory contains the blastx data for non-elegans species\n";
printf $rlfh "                                                    (reduces the size of the main database)\n";
printf $rlfh "MULTI_SPECIES DIR - miscellaneous files with data for multiple species\n";
printf $rlfh "     - wormpep_clw.WS$ver.sql.bz2 - ClustalW protein multiple alignments\n";
printf $rlfh "ONTOLOGY DIR             - gene_associations, obo files for (phenotype GO anatomy) and associated association files\n";
printf $rlfh "\n\n";
printf $rlfh "Release notes on the web:\n";
printf $rlfh "-------------------------\n";
printf $rlfh "http://www.wormbase.org/about/release_schedule\n\n\n";

printf $rlfh "-=====================================================================================-\n";
printf $rlfh "-=========== C. elegans data summary =================================================-\n";
printf $rlfh "-=====================================================================================-\n\n";

#
# Genome version
#
printf $rlfh "Genome version\n";
printf $rlfh "--------------\n\n";

my ($ele_ass) = grep { $_->{is_canonical} } @{ $assembly_stats->{c_elegans}->{assemblies}};
printf $rlfh "The version of the C. elegans reference genome included with this release is:\n\n";
printf $rlfh "Version Name: %s\n", $ele_ass->{assembly_name};
printf $rlfh "INSDC accession: %s\n", $ele_ass->{assembly_accession};
printf $rlfh "UNSC name: %s\n", $ele_ass->{ucsc_name} if exists $ele_ass->{ucsc_name};
print $rlfh "\n";
printf $rlfh "This version has been present in WormBase since release %s\n\n", $ele_ass->{appeared_in};

#
# Gene summary stats
#
&slurp_in_file($rlfh, $gene_stats_file);

#
# Wormpep summary stats
#
&slurp_in_file($rlfh, $wormpep_stats_file);

#
# C. elegans CDS confirmation status
#
&write_cds_confirmation_status($rlfh, $cds_status{elegans}, "C. elegans");

#
# C. elegans operon stats
#
printf $rlfh "\nC. elegans Operons Stats\n";
print  $rlfh "------------------------\n";
printf $rlfh "Live Operons        %d\n", $operon_data{operons};
printf $rlfh "Genes in Operons    %d\n\n", $operon_data{genes_in_operons};

#
# C. elegans GO stats
#
$log->write_to("Writing GO stats...\n");
&write_GO_stats($rlfh, $go_stats{elegans});


#
# Assorted summaries for other species
#
printf $rlfh "\n-=============================================================================-\n";
printf $rlfh "-=========== Other core species data summary =================================-\n";
printf $rlfh "-=============================================================================-\n";

#
# Genome versions
#
print $rlfh "\nAssembly versions\n";
printf $rlfh "----------------\n";
foreach my $acc (sort { $a->full_name cmp $b->full_name } values %accessors) {
  next if $acc->species eq 'elegans';
  my $g_species = $acc->full_name(-g_species => 1);
  my ($ass) = grep { $_->{is_canonical} } @{ $assembly_stats->{$g_species}->{assemblies}};
  printf $rlfh "%-28s %-25s (project %s, current since %s)\n", $acc->full_name, $ass->{assembly_name},  $ass->{bioproject}, $ass->{appeared_in};
}


printf $rlfh "\nApproved gene symbols\n";
printf $rlfh "---------------------\n";
foreach my $spec (sort { $accessors{$a}->full_name cmp $accessors{$b}->full_name } keys %cgc_counts) {
  next if $spec eq 'elegans';
  my $count = $cgc_counts{$spec};
  printf $rlfh "%-28s %5d\n", $accessors{$spec}->full_name, $count;
}

print $rlfh "\nGene counts\n";
print $rlfh "-----------\n";
foreach my $spec (sort { $accessors{$a}->full_name cmp $accessors{$b}->full_name } keys %gene_counts) {
  next if $spec eq 'elegans';
  my $gcount = $gene_counts{$spec}->{genes};
  my $cod_count = $gene_counts{$spec}->{coding_genes};
  printf $rlfh "%-28s %5d (%d coding)\n", $accessors{$spec}->full_name, $gcount, $cod_count;
}

foreach my $spec (sort { $accessors{$a}->full_name cmp $accessors{$b}->full_name } keys %cds_status) {
  next if $spec eq 'elegans';
  print STDERR "Species=$spec\n";
  &write_cds_confirmation_status($rlfh, $cds_status{$spec}, $accessors{$spec}->full_name);
}


printf $rlfh "\n\n-==============================================================================-\n";
printf $rlfh "-=========== News for this release ============================================-\n";
printf $rlfh "-==============================================================================-\n\n";
printf $rlfh "New data sets\n";
printf $rlfh "--------------\n\n\n";
printf $rlfh "New/updated reference genomes\n";
print  $rlfh "------------------------------------\n\n\n";
printf $rlfh "Proposed Changes / Forthcoming Data\n";
print  $rlfh "------------------------------------\n\n\n";
printf $rlfh "Model Changes\n";
printf $rlfh "--------------\n\n";
print  $rlfh "Model changes for this release are documented here:\n\n";
print  $rlfh "http://wiki.wormbase.org/index.php/WS${ver}_Models.wrm\n\n";
printf $rlfh "For more information mail help\@wormbase.org\n\n";

# Installation guide
printf $rlfh "-==============================================================================-\n";
printf $rlfh "-=========== Installation guide ===============================================-\n";
printf $rlfh "-==============================================================================-\n";
printf $rlfh "\n\n";
printf $rlfh "Quick installation guide for UNIX/Linux systems\n";
printf $rlfh "-----------------------------------------------\n\n";
printf $rlfh "1. Create a new directory to contain your copy of WormBase,\n\te.g. /users/yourname/wormbase\n\n";
printf $rlfh "2. Unpack and untar all of the database.*.tar.gz files into\n\tthis directory. You will need approximately 50-60 Gb of disk space.\n\n";
printf $rlfh "3. Obtain and install a suitable acedb binary for your system\n\t(available from www.acedb.org).\n\n";
printf $rlfh "4. Use the acedb 'xace' program to open your database, e.g.\n\ttype 'xace /users/yourname/wormbase' at the command prompt.\n\n";
printf $rlfh "5. See the acedb website for more information about acedb and\n\tusing xace.\n\n";

printf $rlfh "____________  END _____________\n";
close($rlfh);

$log->write_to("DONT FORGET TO FILL IN THE LAST FEW FIELDS IN THE LETTER\n found at $release_letter\n");

$wormbase->check_file($release_letter, $log,
                      minsize => 5500);

$log->mail();
exit(0);


##############################################################
sub slurp_in_file {
  my ($outfh, $file) = @_;

  open(my $read_fh, $file) or $log->log_and_die("Could not open $file for reading\n");
  while(<$read_fh>) {
    print $outfh $_;
  }
  print $outfh "\n";
}

##############################################################
sub read_json {
  my ($json_file) = @_;

  my $json = JSON->new->allow_nonref;
 
  my $json_text = "";
  open(my $json_fh, $json_file) or $log->log_and_die("Could not open $json_file for reading\n");
  while(<$json_fh>) {
    $json_text .= $_;
  }
  my $obj_ref = $json->decode( $json_text );

  return $obj_ref;
}

#########################################
sub get_cds_confirmation_status {
  my ($wb, $db) = @_;
  
  my $full_name = $wb->full_name;
    
  my $query = "Find Gene WHERE (Corresponding_CDS OR Corresponding_transcript OR Corresponding_pseudogene) AND CGC_name AND Species = \"$full_name\"";
  my $gene_seq_count = $db->fetch(-query=> "$query");
  
  # wormpep status overview
  my %results;

  #species specific wormpeps
  my $wormpep_datafile = $wb->wormpep."/".$wb->pepdir_prefix."pep".$wb->version;
  open(my $wpfh, $wormpep_datafile) or $log->log_and_die("Could not open $wormpep_datafile\n");
  while(<$wpfh>) {
    next if not /^\>/;
    /Confirmed/ and $results{confirmed}++;
    /Partially_confirmed/ and $results{supported}++;
    /Predicted/ and $results{predicted}++;
    $results{total}++;
  }

  return \%results;
}

#########################################
sub write_cds_confirmation_status {
  my ($fh, $hash, $species) = @_;

  my $total = $hash->{total};
  my $conf = $hash->{confirmed};
  my $conf_perc = ($conf / $total) * 100;
  my $supp = $hash->{supported};
  my $supp_perc = ($supp / $total) * 100;
  my $pred = $hash->{predicted};
  my $pred_perc = ($pred / $total) * 100;

  printf $fh "\n$species Gene model confirmation status (based on the EST/mRNA/RNASeq evidence)\n";
  print  $fh "------------------------------------------------------------\n";
  printf $fh "Confirmed             %5d (%2.1f%%)	Every base of every exon has transcription evidence (mRNA/EST/RNASeq)\n", $conf, $conf_perc;
  printf $fh "Partially_confirmed   %5d (%2.1f%%)	Some, but not all exon bases are covered by transcript evidence\n", $supp, $supp_perc;
  printf $fh "Predicted             %5d (%2.1f%%)	No coverage by mRNA/EST/RNASeq evidence\n", $pred, $pred_perc;

}

#############################################
sub get_GO_stats {
  my ($wb) = @_;

  my (%go_stats, %gene2code,%go2code);
  my (%assocs_by_code,%assocs_by_source,%assocs_by_paper);
  my ($iea_interpro, $iea_other);

  my $go_annot_file = $wb->ontology . "/gene_association.WS${ver}.wb." . $wb->full_name(-g_species => 1 );
  open(my $gaf, $go_annot_file) or $log->log_and_die("Could not open $go_annot_file for reading");
  while(<$gaf>) {
    next if /^\!/;
    my @l = split(/\t/, $_);
    my ($gene, $term, $ref, $code, $source,) = @l[1,4,5,6,14];

    $term =~ s/GO://;
      
    $gene2code{$gene}->{$code} = 1;
    $go2code{$term}->{$code} = 1;
    $assocs_by_code{$code}++;
    $assocs_by_source{$source}++;
    $assocs_by_paper{$ref}++;

    if ($code eq 'IEA') {
      if ($ref eq 'GO_REF:0000002') {
        $iea_interpro++;
      } else {
        $iea_other++;
      }
    }
  }
  
  #
  # Classify *assocation* as IEA vs. non-IEA
  #
  my $total_assocs = 0;
  my $iea_assocs = 0;
  my $non_iea_assocs = 0;
  foreach my $code (keys %assocs_by_code) {
    $total_assocs += $assocs_by_code{$code};
    if ($code eq 'IEA') {
      $iea_assocs += $assocs_by_code{$code};
    } else {
      $non_iea_assocs += $assocs_by_code{$code};
    }
  }
  
  #
  # Classify *genes* as having IEA-only terms, non-IEA-only terms, or both
  #
  my $IEA_gc = 0;
  my $non_IEA_gc = 0;
  my $IEA_and_non_IEA_gc = 0;

  foreach my $gene (sort keys %gene2code) {
    my @non_iea = grep { $_ ne 'IEA' } keys %{$gene2code{$gene}};
    my @iea = grep { $_ eq 'IEA' } keys %{$gene2code{$gene}};
    
    if (@non_iea and @iea) {
      $IEA_and_non_IEA_gc++;
    } elsif (@non_iea) {
      $non_IEA_gc++;
    } elsif (@iea) {
      $IEA_gc++;
    }
  }
  
  #
  # Classify *GO_terms* as having IEA-only genes, non-IEA-only genes, or both
  #
  my $IEA_go_count = 0;
  my $non_IEA_go_count = 0;
  my $IEA_and_non_IEA_go_count = 0;

  foreach my $go_term (keys %go2code) {
    my @codes = keys %{$go2code{$go_term}};
    my @iea = grep { $_ eq 'IEA' } @codes;
    my @non_iea = grep { $_ ne 'IEA' } @codes;
      
    if (@iea and @non_iea) {
      $IEA_and_non_IEA_go_count++;
    } elsif (@iea) {
      $IEA_go_count++;
    } elsif (@non_iea) {
      $non_IEA_go_count++;
    }
  }

  #
  # encode as JSON and write it out
  # 
  $go_stats{annotations}->{by_provider_counts} = \%assocs_by_source;
  $go_stats{annotations}->{by_code_counts}     = \%assocs_by_code;
  #$go_stats{annotations}->{by_paper_counts}    = \$assocs_by_paper;

  $go_stats{summaries}->{annotation_counts}->{total} = $total_assocs;
  $go_stats{summaries}->{annotation_counts}->{iea} = $iea_assocs;
  $go_stats{summaries}->{annotation_counts}->{non_iea} = $non_iea_assocs;
  $go_stats{summaries}->{annotation_counts}->{iea_interpro} = $iea_interpro;
  $go_stats{summaries}->{annotation_counts}->{iea_non_interpro} = $iea_other;
  
  $go_stats{summaries}->{gene_counts}->{iea_only} = $IEA_gc;
  $go_stats{summaries}->{gene_counts}->{non_iea_only} = $non_IEA_gc;
  $go_stats{summaries}->{gene_counts}->{iea_and_non_iea} = $IEA_and_non_IEA_gc;

  $go_stats{summaries}->{go_term_counts}->{iea_only} = $IEA_go_count;
  $go_stats{summaries}->{go_term_counts}->{non_iea_only} = $non_IEA_go_count;
  $go_stats{summaries}->{go_term_counts}->{iea_and_non_iea} = $IEA_and_non_IEA_go_count;

  return \%go_stats;
}



#########################################
sub write_GO_stats {
  my ($fh, $go_stats ) = @_;
  
  my %summa  = %{$go_stats->{summaries}};
  my %annots = %{$go_stats->{annotations}};

  printf $fh "C. elegans GO annotation status\n";
  printf $fh "-------------------------------\n\n";
  
  printf $fh "GO_codes - used for assigning evidence\n";
  printf $fh "  IBA Inferred by Biological aspect of Ancestor\n";
  printf $fh "  IC  Inferred by Curator\n";
  printf $fh "  IDA Inferred from Direct Assay\n";
  printf $fh "  IEA Inferred from Electronic Annotation\n";
  printf $fh "  IEP Inferred from Expression Pattern\n";
  printf $fh "  IGI Inferred from Genetic Interaction\n";
  printf $fh "  IKR Inferred from Key Residues\n";
  printf $fh "  IMP Inferred from Mutant Phenotype\n";
  printf $fh "  IPI Inferred from Physical Interaction\n";
  printf $fh "  IRD Inferred from Rapid Divergence\n";
  printf $fh "  ISM Inferred from Sequence Model\n";
  printf $fh "  ISO Inferred from Sequence Orthology\n";
  printf $fh "  ISS Inferred from Sequence (or Structural) Similarity\n";
  printf $fh "  NAS Non-traceable Author Statement\n";
  printf $fh "  ND  No Biological Data available\n";
  printf $fh "  RCA Inferred from Reviewed Computational Analysis\n";
  printf $fh "  TAS Traceable Author Statement\n\n";
  
  printf $fh "Number of gene<->GO_term associations    %5d\n", $summa{annotation_counts}->{total};
  printf $fh "  Breakdown by annotation provider:\n";
  printf $fh "    %-20s %5d\n", "WormBase", $annots{by_provider_counts}->{"WB"};
  foreach my $k (sort { $annots{by_provider_counts}->{$b} <=> $annots{by_provider_counts}->{$a} } keys %{$annots{by_provider_counts}}) {
    next if $k eq 'WB';
    printf $fh "    %-20s %5d\n",  $k,  $annots{by_provider_counts}->{$k};
  }
  printf $fh "  Breakdown by evidence code:\n";
  printf $fh "    IEA     %d\n", $summa{annotation_counts}->{iea};
  printf $fh "      Interpro2GO %5d\n", $summa{annotation_counts}->{iea_interpro}; 
  printf $fh "      Other       %5d\n", $summa{annotation_counts}->{iea_non_interpro}; 
  printf $fh "    non-IEA %d\n", $summa{annotation_counts}->{non_iea};
  foreach my $k (sort keys %{$annots{by_code_counts}}) {
    next if $k eq 'IEA';
    printf $fh "      %-3s   %5d\n", $k, $annots{by_code_counts}->{$k};
  }
  
  printf $fh "\nGenes Stats:\n";
  printf $fh "  Genes with GO_term connections  %5d \n", $summa{gene_counts}->{non_iea_only} + $summa{gene_counts}->{iea_only} + $summa{gene_counts}->{iea_and_non_iea};
  printf $fh "    Non-IEA-only annotation            %5d\n", $summa{gene_counts}->{non_iea_only};
  printf $fh "    IEA-only annotation                %5d\n", $summa{gene_counts}->{iea_only};
  printf $fh "    Both IEA and non-IEA annotations   %5d\n", $summa{gene_counts}->{iea_and_non_iea};
  
  printf $fh "\nGO_term Stats:\n";
  printf $fh "  Distinct GO_terms connected to Genes  %5d\n", $summa{go_term_counts}->{non_iea_only} + $summa{go_term_counts}->{iea_only} + $summa{go_term_counts}->{iea_and_non_iea};
  printf $fh "    Associated by non-IEA only              %5d\n", $summa{go_term_counts}->{non_iea_only};
  printf $fh "    Associated by IEA only                  %5d\n", $summa{go_term_counts}->{iea_only};
  printf $fh "    Associated by both IEA and non-IEA      %5d\n", $summa{go_term_counts}->{iea_and_non_iea};
}

##########################################

__END__

