#!/usr/local/bin/perl5.8.0 -w
#
# release_letter.pl                            
# 
# by Anthony Rogers                             
#
# Last updated by: $Author: klh $               
# Last updated on: $Date: 2013-04-02 11:32:38 $

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
use Log_files;
use Storable;
use Ace;
use Modules::Remap_Sequence_Change;
use Species;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($opt_c, $opt_l);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"      => \$store,
	    "c"          => \$opt_c,
	    "l"          => \$opt_l,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# establish log file.
my $log = Log_files->make_build_log($wormbase);

# in test mode?
if ($test) {
  $log->write_to("In test mode\n") if ($verbose);
}


##############
# variables  #                                                                   
##############

my $basedir         = $wormbase->basedir;     # BASE DIR
my $ace_dir         = $wormbase->autoace;     # AUTOACE DATABASE DIR
my $wormpep_dir     = $wormbase->wormpep;     # CURRENT WORMPEP
my $reports_dir     = $wormbase->reports;     # AUTOACE REPORTS
my $tace = $wormbase->tace;
my $db = Ace->connect(-path  => $ace_dir,
		      -program =>$tace) || $log->log_and_die("Connection failure: ",Ace->error);
my $ver     = $wormbase->get_wormbase_version;
my $old_ver = $ver -1;

my $date        = `date`;

$wormbase->release_composition($log) if defined($opt_c);

# make the release letter
if( defined($opt_l)) {
  my $release_letter = "$reports_dir/letter.WS$ver";
  open (RL,">$release_letter");
  printf RL "New release of WormBase WS$ver\n\n";
  printf RL "WS$ver was built by [INSERT NAME HERE]\n";
  printf RL "-===================================================================================-\n";
  printf RL "The WS$ver build directory includes:\n";
  printf RL "species/ DIR              -  contains a sub dir for each WormBase species (G_SPECIES)\n";
  printf RL "species/G_SPECIES DIR     -  contains a sub dir for each NCBI genome sequencing BioProject (BIOPROJECT) for the species, with the following files:\n";
  printf RL "     - G_SPECIES.BIOPROJECT.WS$ver.genomic.fa.gz                - Unmasked genomic DNA\n";
  printf RL "     - G_SPECIES.BIOPROJECT.WS$ver.genomic_masked.fa.gz         - Hard-masked (repeats replaced with Ns) genomic DNA\n";
  printf RL "     - G_SPECIES.BIOPROJECT.WS$ver.genomic_softmasked.fa.gz     - Soft-masked (repeats lower-cased) genomic DNA\n";
  printf RL "     - G_SPECIES.BIOPROJECT.WS$ver.protein.fa.gz                - Current live protein set\n";
  printf RL "     - G_SPECIES.BIOPROJECT.WS$ver.cds_transcripts.fa.gz        - Spliced cDNA sequence for the CDS portion of protein-coding transcripts\n";
  printf RL "     - G_SPECIES.BIOPROJECT.WS$ver.ncrna_transcripts.fa.gz      - Spliced cDNA sequence for non-coding RNA transcripts\n";
  printf RL "     - G_SPECIES.BIOPROJECT.WS$ver.intergenic_sequences.fa.gz   - DNA sequence between pairs of adjacent genes\n";
  printf RL "     - G_SPECIES.BIOPROJECT.WS$ver.annotations.gff[2|3].gz      - Sequence features in either GFF2 or GFF3 format\n";
  printf RL "     - G_SPECIES.BIOPROJECT.WS$ver.ests.fa.gz                   - ESTs and mRNA sequences extracted from the public databases\n";
  printf RL "     - G_SPECIES.BIOPROJECT.WS$ver.best_blastp_hits.txt.gz      - Best blastp matches to human, fly, yeast, and non-WormBase Uniprot proteins\n";
  printf RL "     - G_SPECIES.BIOPROJECT.WS$ver.*pep_package.tar.gz          - latest version of the [worm|brig|bren|rema|jap|ppa|brug]pep package (if updated since last release)\n";
  printf RL "     - annotation/                    - contains additional annotations:\n";
  printf RL "        - G_SPECIES.BIOPROJECT.WS$ver.confirmed_genes.txt.gz              - DNA sequences of all genes confirmed by EST &/or cDNA\n";
  printf RL "        - G_SPECIES.BIOPROJECT.WS$ver.cDNA2orf.txt.gz                     - Latest set of ORF connections to each cDNA (EST, OST, mRNA)\n";      
  printf RL "        - G_SPECIES.BIOPROJECT.WS$ver.geneIDs.txtgz                       - list of all current gene identifiers with CGC & molecular names (when known)\n";
  printf RL "        - G_SPECIES.BIOPROJECT.WS$ver.PCR_product2gene.txt.gz             - Mappings between PCR products and overlapping Genes\n";
  printf RL "        - G_SPECIES.BIOPROJECT.WS$ver.*oligo_mapping.txt.gz               - Oligo array mapping files\n";
  printf RL "        - G_SPECIES.BIOPROJECT.WS$ver.knockout_consortium_alleles.xml.gz  - Table of Knockout Consortium alleles\n";
  printf RL "        - G_SPECIES.BIOPROJECT.WS$ver.SRA_gene_expression.tar.gz          - Tables of gene expression values computed from SRA RNASeq data\n";

  printf RL "acedb DIR                -  Everything needed to generate a local copy of the The Primary database\n";
  printf RL "     - database.WS$ver.*.tar.gz   - compressed acedb database for new release\n";
  printf RL "     - models.wrm.WS$ver          - the latest database schema (also in above database files)\n";
  printf RL "     - WS$ver-WS$old_ver.dbcomp   - log file reporting difference from last release\n";
  printf RL "     - *Non_C_elegans_BLASTX/          - This directory contains the blastx data for non-elegans species\n";
  printf RL "                                                    (reduces the size of the main database)\n";
  printf RL "COMPARATIVE_ANALYSIS DIR - comparative analysis files\n";
  printf RL "     - compara.WS$ver.tar.bz2     - gene-tree and alignment GFF files\n";
  printf RL "     - wormpep_clw.WS$ver.sql.bz2 - ClustalW protein multiple alignments\n";
  printf RL "ONTOLOGY DIR             - gene_associations, obo files for (phenotype GO anatomy) and associated association files\n";
  printf RL "\n\n";
  printf RL "Release notes on the web:\n-------------------------\n";
  printf RL "http://www.wormbase.org/wiki/index.php/Release_Schedule\n\n\n\n";
  
  # make the chromosomal sequence changes file
  $log->write_to("\nGenerating C. elegans chromosomal sequence changes Data\n\n");
  open (CC, "> $reports_dir/chromosome_changes") || die "Can't open file $reports_dir/chromosome_changes\n";
  my $assembly_mapper = Remap_Sequence_Change->new($ver-1, $ver, $wormbase->species, $wormbase->genome_diffs);
  my $text = $assembly_mapper->write_changes( $wormbase );
  printf CC $text;
  close(CC);

  my @release_files = ("$reports_dir/chromosome_changes","$reports_dir/genedata","$reports_dir/wormpep","$reports_dir/composition");
  
  #include all the pre-generated reports
  my $file = shift(@release_files);
  while (defined($file)) {
    
    open (READIN, "<$file") || die "cant open $file\n";
    printf RL "C. elegans ";
    #    if ($file eq  $reports_dir."/composition") {printf RL "C. elegans ";}
    while(<READIN>) {
      print RL "$_";
    }
    close READIN;
    printf RL "\n\n";
    $file = shift(@release_files);
  }
  
  ## For all curated/gemones with a gene set do the following ##
  my %accessors = ($wormbase->species_accessors);
  $accessors{$wormbase->species} = $wormbase;
  my $elegansAccessor= $accessors{elegans};
  delete $accessors{'elegans'}; #this aviods duplicating stats for elegans as we also pull in a preprepared report.
  my @wormpep_species = (keys%accessors);
  my $wormpep_species;
  my ($species,$name);
  foreach $wormpep_species(@wormpep_species) {
    $species = $accessors{$wormpep_species};
    $name = $species->full_name;

    $log->write_to("Getting $name genomic sequence info\n\n");
    #Do a genome composition to follow the elegans composition above.
    printf RL "$name Genome sequence composition:\n----------------------------\n";
    $file = $species->chromosomes."/composition.all";
    
    if (defined($file)) {
      open (READIN, "<$file") || die "cant open $file\n";
      while(<READIN>) {
	printf RL "$_";
      }
      close READIN;
      printf RL "\n\n";
    }
  }

  ######################################
  #  Tier II Species Gene Stats        #
  ######################################
  $log->write_to("Retrieving Tier II Species Gene Stats\n\n");
  my %tierII_species = ($wormbase->species_accessors);
  my @tierII = (keys%tierII_species);
  my $tierII;
  printf RL "\n\nTier II Gene counts\n";
  printf RL "---------------------------------------------\n";
  foreach $tierII(@tierII) {
    my $gene_count_query = "Find Gene where Species = \"*${tierII}*\" AND Live";
    my $gene_count = $db->fetch(-query=> "$gene_count_query");
    #my $Coding_count_query = "Find CDS where Species = \"*${tierII}*\" AND method = \"curated\"";
    # now count coding genes, rather than coding isoforms
    my $Coding_count_query = "Find Gene where Species = \"*${tierII}*\" AND Corresponding_CDS";
    my $Coding_count = $db->fetch(-query=> "$Coding_count_query");
    printf RL "$tierII Gene count $gene_count (Coding ${Coding_count})\n";
  }
  printf RL "---------------------------------------------\n\n\n";


  ######################################
  #  CDS Status stats                  #
  ######################################

  #add back in elegans for this step.
  push(@wormpep_species,"elegans");
  foreach $wormpep_species(@wormpep_species) {
    $species = $wormpep_species eq 'elegans' ? $elegansAccessor : $accessors{$wormpep_species};
    $name = $species->full_name;
    my %wp_status;
    my $wormpep_datafile;
    # Find out Gene->CDS, Transcript, Pseudogene connections
    
    my $query = "Find Gene WHERE (Corresponding_CDS OR Corresponding_transcript OR Corresponding_pseudogene) AND CGC_name AND Species = \"$name\"";
    my $gene_seq_count = $db->fetch(-query=> "$query");

    # wormpep status overview
    $log->write_to("Generating $name wormpep overview stats\n\n");
    #species specific wormpeps
    $wormpep_datafile = $species->wormpep."/".$species->pepdir_prefix."pep".$species->version;
    $wp_status{Confirmed}  = `grep Confirmed    $wormpep_datafile | wc -l`;
    $wp_status{Supported}  = `grep confirmed    $wormpep_datafile | wc -l`;
    $wp_status{Predicted}  = `grep Predicted    $wormpep_datafile | wc -l`;
    $wp_status{Gene}       = $gene_seq_count;
    $wp_status{Uniprot}    = `grep 'UniProt:'        $wormpep_datafile | wc -l`;
    $wp_status{Protein_ID} = `grep 'protein_id' $wormpep_datafile | wc -l`;
    $wp_status{Total}      = $wp_status{Confirmed} + $wp_status{Supported} + $wp_status{Predicted}; 
    
    printf RL "\n\n";
    printf RL "-------------------------------------------------\n";
    printf RL "$name Protein Stats:\n";
    printf RL "-------------------------------------------------\n";
    printf RL "Status of entries: Confidence level of prediction (based on the amount of transcript evidence)\n";
    printf RL "-------------------------------------------------\n";
    printf RL "Confirmed            %6d (%2.1f%%)\tEvery base of every exon has transcription evidence (mRNA, EST etc.)\n", $wp_status{Confirmed}, (($wp_status{Confirmed}/$wp_status{Total}) * 100);
    printf RL "Partially_confirmed  %6d (%2.1f%%)\tSome, but not all exon bases are covered by transcript evidence\n", $wp_status{Supported}, (($wp_status{Supported}/$wp_status{Total}) * 100);
    printf RL "Predicted            %6d (%2.1f%%)\tNo transcriptional evidence at all\n", $wp_status{Predicted}, (($wp_status{Predicted}/$wp_status{Total}) * 100);
    printf RL "\n\n\n";
    printf RL "Status of entries: Protein Accessions\n-------------------------------------\nUniProtKB accessions %6d (%2.1f%%)\n\n", $wp_status{Uniprot}, (($wp_status{Uniprot}/$wp_status{Total}) * 100) if ($wp_status{Uniprot} >0);
    printf RL "Status of entries: Protein_ID's in EMBL\n---------------------------------------\nProtein_id           %6d (%2.1f%%)\n\n", $wp_status{Protein_ID}, (($wp_status{Protein_ID}/$wp_status{Total}) * 100) if ($wp_status{Protein_ID} >0);
    printf RL "Gene <-> CDS,Transcript,Pseudogene connections\n";
    printf RL "----------------------------------------------\n";
    printf RL "$name entries with WormBase-approved Gene name %6d\n", $wp_status{Gene};
    printf RL "\n\n";
  }
  
  ######################################
  #  Get the Operon stats    - Paul    #
  ######################################
      $log->write_to("\nRetrieving Operon Data\n");
  my $operon_query = "Find Operon WHERE method = \"Operon\" AND Species = \"*elegans\"";
  my $gene_query = "Find Gene where Contained_in_operon AND NEXT";
  my $Operon_count = $db->fetch(-query=> "$operon_query");
  my $Operon_genes = $db->fetch(-query=> "$gene_query");
  

  printf RL "C. elegans Operons Stats\n";
  printf RL "---------------------------------------------\n";
  printf RL "Description: These exist as closely spaced gene clusters similar to bacterial operons\n";
  printf RL "---------------------------------------------\n";
  printf RL "| Live Operons        $Operon_count                |\n";
  printf RL "| Genes in Operons    $Operon_genes                |\n";
  printf RL "---------------------------------------------\n\n\n";



  ######################################
  #  Get the GO annotations  - Ranjana #
  ######################################

  #Table maker definition
  #SCRIPT:GeneGO_codes.def
  # WBgene00000001 GO:0000001 IEA
  if (-e "$basedir/autoace/wquery/SCRIPT:GeneGO_codes.def") {
    my $def = "$basedir/autoace/wquery/SCRIPT:GeneGO_codes.def";
    my $command = "Table-maker -p $def\nquit\n";
    $log->write_to("\nRetrieving GO::Gene data, using Table-maker...\n");
    my $count = "0";
    my (%Gene2code,%GO2code,%GO2gene,%gene2GO);

    open (TACE, "echo '$command' | $tace $ace_dir | ") || die "Cannot query acedb. $command  $tace\n";
    while (<TACE>) {
      chomp;
      s/\"//g;
      s/GO://g;
      unless (/(WBGene\d+)\s+(\d+)\s+(\w+)/) {next;}
      if (/(WBGene\d+)\s+(\d+)\s+(\w+)/){
	#query with WBGene and GO_code is returned
	$Gene2code{$1} = $3;
	#query with GO term and GO_code is returned
	$GO2code{$2} = $3;
	#query with GO code and get Genes that it connects to.
	# this stores all of the values in an array which makes it complicated to get a total.
	# so count the number times something is pushed onto the value arrays and you will have the 
	#total number of Gene::GO connections
	$count++;
	push @{$gene2GO{$1}} , $2;
	#query with WBGene gived GO_terms connected.
	push @{$GO2gene{$2}} , $1;
      } 
    }
    close TACE;

    # no. of genes with GO terms
    my $query1 = "Find Gene WHERE GO_term";
    my $gc = $db->fetch(-query=> "$query1");
    
    #RNAi_mapping::Gene GO associations
    my $rnai = `cat $basedir/autoace/acefiles/RNAi_mappings.ace | tr -d "\n" | sed 's/Gene : /@/g' | tr "@" "\n" | grep GO_term | sed 's/GO_term /@/' | tr "@" "\n" | tr -d '"' | grep WBGene | sort -u | wc -l`;
    chomp $rnai;
    #Citace::Gene GO associations
    my $citace = `cat $basedir/autoace/acefiles/primaries/citace/caltech_GO_term.ace | tr -d "\n" | sed 's/Gene : /@/g' | tr "@" "\n" | grep GO_term | sed 's/GO_term /@/' | tr "@" "\n" | tr -d '"' | grep WBGene | sort -u | wc -l`;
    chomp $citace;
    #Inherit_GO_terms:Gene GO associations
    my $inherit = `cat $basedir/autoace/acefiles/inherited_GO_terms.ace | tr -d "\n" | sed 's/Gene : /@/g' | tr "@" "\n" | grep GO_term | sed 's/GO_term /@/' | tr "@" "\n" | tr -d '"' | grep WBGene | sort -u | wc -l`;
    chomp $inherit;
    # of genes with IEA GO terms
    my $IEA_gc  = grep {/IEA/} values %Gene2code;
    
    # need to know how many Gene keys for next $non_IEA_gc calc.
    my $Genekeys = keys %Gene2code;

    # of genes with non_IEA GO terms
    my $non_IEA_gc = $Genekeys - $IEA_gc;

    # How many Gene -> GO_term connections are there?
    my $Genevalues = values %Gene2code;

    #and the same for the annotations
    # of GO annotations
    my $query2  = "Find GO_term";
    my $GO_annotations = $db->fetch(-query=> "$query2");
    # GO2Gene
    my $query3  = "Find GO_term WHERE Gene";
    my $GOwithgene = $db->fetch(-query=> "$query3");
    # of IEA GO annotations
    my $IEA = grep {/IEA/} values %GO2code;
    my $IC = grep {/IC/} values %GO2code;
    my $IDA = grep {/IDA/} values %GO2code;
    my $IEP = grep {/IEP/} values %GO2code;
    my $IGI = grep {/IGI/} values %GO2code;
    my $IMP = grep {/IMP/} values %GO2code;
    my $IPI = grep {/IPI/} values %GO2code;
    my $ISS = grep {/ISS/} values %GO2code;
    my $NAS = grep {/NAS/} values %GO2code;
    my $ND = grep {/ND/} values %GO2code;
    my $RCA = grep {/RCA/} values %GO2code;
    my $TAS = grep {/TAS/} values %GO2code;

    # need to know how many GO_terms keys for next $non_IEAno calc.
    my $GOkeys = keys %GO2code;

    # How many GO_term -> Gene connections are there?
    my $GOvalues = values %GO2code;

    # of non-IEA GO annotations
    my $nonIEAno = $GOkeys - $IEA; 

    printf RL "GO Annotation Stats WS$ver\n--------------------------------------\n\n";

    printf RL "GO_codes - used for assigning evidence\n";
    printf RL "--------------------------------------\n";
    printf RL "IC  Inferred by Curator\n";
    printf RL "IDA Inferred from Direct Assay\n";
    printf RL "IEA Inferred from Electronic Annotation\n";
    printf RL "IEP Inferred from Expression Pattern\n";
    printf RL "IGI Inferred from Genetic Interaction\n";
    printf RL "IMP Inferred from Mutant Phenotype\n";
    printf RL "IPI Inferred from Physical Interaction\n";
    printf RL "ISS Inferred from Sequence (or Structural) Similarity\n";
    printf RL "NAS Non-traceable Author Statement\n";
    printf RL "ND  No Biological Data available\n";
    printf RL "RCA Inferred from Reviewed Computational Analysis\n";
    printf RL "TAS Traceable Author Statement\n";
    printf RL "------------------------------------------------\n\n";

    printf RL "Total number of Gene::GO connections:  $count\n\n"; 

    printf RL "Genes Stats:\n";
    printf RL "----------------\n";
    printf RL "Genes with GO_term connections         $gc  \n";
    printf RL "           IEA GO_code present         $IEA_gc  \n";
    printf RL "       non-IEA GO_code present         $non_IEA_gc  \n\n";
    
    printf RL "Source of the mapping data             \n";
    printf RL "Source: *RNAi (GFF mapping overlaps)   $rnai  \n";
    printf RL "        *citace                        $citace  \n";
    printf RL "        *Inherited (motif & phenotype) $inherit  \n\n";

    printf RL "GO_terms Stats:\n";
    printf RL "---------------\n";
    printf RL "Total No. GO_terms                     $GO_annotations  \n";
    printf RL "GO_terms connected to Genes            $GOwithgene  \n";
    printf RL "GO annotations connected with IEA      $IEA  \n";
    printf RL "GO annotations connected with non-IEA  $nonIEAno  \n";
    printf RL "   Breakdown  IC - $IC   IDA - $IDA   ISS - $ISS \n";
    printf RL "             IEP - $IEP   IGI - $IGI   IMP - $IMP \n";
    printf RL "             IPI - $IPI  NAS - $NAS     ND  - $ND  \n";
    printf RL "             RCA - $RCA   TAS - $TAS   \n\n\n";
  }
  else {
    $log->write_to("\nERROR - GeneGO_codes.def abscent from autoace/wquery\nThese stats will be missing from the release letter\n\n");
  }




  
  printf RL "-===================================================================================-\n";
  $log->write_to("\nUseful Gene Stats\n");
  printf RL "\nUseful Stats:\n---------\n\n";
  
  my $gene_seq_cgc_q = "Find Gene WHERE (Corresponding_CDS OR Corresponding_transcript OR Corresponding_pseudogene) AND CGC_name";
  my $ele_seq_cgc_q = "Find Gene WHERE (Corresponding_CDS OR Corresponding_transcript OR Corresponding_pseudogene) AND CGC_name AND Species = \"Caenorhabditis elegans\"";
  my $briggsae_seq_cgc_q = "Find Gene WHERE (Corresponding_CDS OR Corresponding_transcript OR Corresponding_pseudogene) AND CGC_name AND Species = \"Caenorhabditis briggsae\""; 
  my $remanei_seq_cgc_q = "Find Gene WHERE (Corresponding_CDS OR Corresponding_transcript OR Corresponding_pseudogene) AND CGC_name AND Species = \"Caenorhabditis remanei\""; 
  my $japonica_seq_cgc_q = "Find Gene WHERE (Corresponding_CDS OR Corresponding_transcript OR Corresponding_pseudogene) AND CGC_name AND Species = \"Caenorhabditis japonica\""; 
  my $brenneri_seq_cgc_q = "Find Gene WHERE (Corresponding_CDS OR Corresponding_transcript OR Corresponding_pseudogene) AND CGC_name AND Species = \"Caenorhabditis brenneri\""; 
  my $pristionchus_seq_cgc_q= "Find Gene WHERE (Corresponding_CDS OR Corresponding_transcript OR Corresponding_pseudogene) AND CGC_name AND Species = \"Pristionchus pacificus\""; 
  my $brugia_seq_cgc_q= "Find Gene WHERE (Corresponding_CDS OR Corresponding_transcript OR Corresponding_pseudogene) AND CGC_name AND Species = \"Brugia malayi\""; 
  
  my $gene_seq_cgc = $db->fetch(-query=> "$gene_seq_cgc_q");
  my $ele = $db->fetch(-query=> "$ele_seq_cgc_q");
  my $briggsae_seq_cgc = $db->fetch(-query=> "$briggsae_seq_cgc_q");
  my $remanei_seq_cgc = $db->fetch(-query=> "$remanei_seq_cgc_q");
  my $japonica_seq_cgc = $db->fetch(-query=> "$japonica_seq_cgc_q");
  my $brenneri_seq_cgc = $db->fetch(-query=> "$brenneri_seq_cgc_q");
  my $pristionchus_seq_cgc = $db->fetch(-query=> "$pristionchus_seq_cgc_q");
  my $brugia_seq_cgc = $db->fetch(-query=> "$brugia_seq_cgc_q");
  
  printf RL "Genes with Sequence and WormBase-approved Gene names\n";
  printf RL "WS$ver $gene_seq_cgc ($ele elegans / $briggsae_seq_cgc briggsae / $remanei_seq_cgc remanei / $japonica_seq_cgc japonica / $brenneri_seq_cgc brenneri / $pristionchus_seq_cgc pristionchus / $brugia_seq_cgc brugia)\n\n\n";
  # Close the database connection now we have finished with it
  $db->close;


printf RL "-===================================================================================-\n";
# User filled sections
printf RL "\n\n\n";
  printf RL "New Data:\n---------\n\n\n";
  printf RL "Genome sequence updates:\n-----------------------\n\n\n";
  printf RL "New Fixes:\n----------\n\n\n";
  printf RL "Known Problems:\n---------------\n\n\n";
  printf RL "Other Changes:\n--------------\n\n";
  printf RL "Proposed Changes / Forthcoming Data:\n-------------------------------------\n\n\n";
  printf RL "Model Changes:\n------------------------------------\n\n\n";
  printf RL "For more info mail help\@wormbase.org\n";

  # Installation guide
  printf RL "-===================================================================================-\n";
  printf RL "\n\n";
  printf RL "Quick installation guide for UNIX/Linux systems\n-----------------------------------------------\n\n";
  printf RL "1. Create a new directory to contain your copy of WormBase,\n\te.g. /users/yourname/wormbase\n\n";
  printf RL "2. Unpack and untar all of the database.*.tar.gz files into\n\tthis directory. You will need approximately 2-3 Gb of disk space.\n\n";
  printf RL "3. Obtain and install a suitable acedb binary for your system\n\t(available from www.acedb.org).\n\n";
  printf RL "4. Use the acedb 'xace' program to open your database, e.g.\n\ttype 'xace /users/yourname/wormbase' at the command prompt.\n\n";
  printf RL "5. See the acedb website for more information about acedb and\n\tusing xace.\n\n";
  
  
  printf RL "____________  END _____________\n";
  close(RL);


  $log->write_to("DONT FORGET TO FILL IN THE LAST FEW FIELDS IN THE LETTER\n found at $release_letter\n");
  

##################
# Check the files
##################

  $wormbase->check_file($release_letter, $log,
                  minsize => 5500);
}



# say goodbye
$log->mail();
$log->write_to("Finished.\n") if ($verbose);
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


__END__

=pod

=head2 NAME - release_letter.pl

=head1 USAGE

=over 4

=item release_letter.pl  [-c -d -l]

=back

This script:

will kindly take away the pain of having to write a release letter at the end of the build
If the script is being run at the end of a build in which the 3 sub parts have been generated then use the B<-l> option.
Otherwise generate the sequence comparison and database comparison sections with the B<-c> and B<-d> options.

I<release_letter.pl MANDATORY arguments:> B<NONE>

I<script_template.pl  OPTIONAL arguments:>


B<-d> create the database comparison file.

B<-c> create the sequence composition comparison file.

B<-l> actually write the letter out.

=back
=over 4

=head1 REQUIREMENTS

=over 4

=item None.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
