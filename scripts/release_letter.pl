#!/usr/local/bin/perl5.8.0 -w
#
# release_letter.pl                            
# 
# by Anthony Rogers                             
#
# Last updated by: $Author: pad $               
# Last updated on: $Date: 2010-07-08 11:43:48 $

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

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);



##############
# variables  #                                                                   
##############

my $basedir         = $wormbase->basedir;     # BASE DIR
my $ace_dir         = $wormbase->autoace;     # AUTOACE DATABASE DIR
my $reports_dir     = $wormbase->reports;     # AUTOACE REPORTS
my $wormpep_dir     = $wormbase->wormpep;     # CURRENT WORMPEP


my $ver     = $wormbase->get_wormbase_version;
my $old_ver = $ver -1;

my $date        = `date`;

my $webdir = "/nfs/WWWdev/SANGER_docs/htdocs/Projects/C_elegans/WORMBASE";

$wormbase->release_composition($log) if defined($opt_c);

# make the release letter
if( defined($opt_l)) {
  my $release_letter = "$reports_dir/letter.WS$ver";
  open (RL,">$release_letter");
  print RL "New release of WormBase WS$ver, Wormpep$ver and Wormrna$ver $date\n\n";
  print RL "WS$ver was built by [INSERT NAME HERE]\n";
  print RL "======================================================================\n\n";
  print RL "This directory includes:\n";
  print RL "i)   database.WS$ver.*.tar.gz    -   compressed data for new release\n";
  print RL "ii)  models.wrm.WS$ver           -   the latest database schema (also in above database files)\n";
  print RL "iii) CHROMOSOMES/subdir         -   contains 3 files (DNA, GFF & AGP per chromosome)\n";
  print RL "iv)  WS$ver-WS$old_ver.dbcomp         -   log file reporting difference from last release\n";
  print RL "v)   wormpep$ver.tar.gz          -   full Wormpep distribution corresponding to WS$ver\n";
  print RL "vi)   wormrna$ver.tar.gz          -   latest WormRNA release containing non-coding RNA's in the genome\n";
  print RL "vii)  confirmed_genes.WS$ver.gz   -   DNA sequences of all genes confirmed by EST &/or cDNA\n";
  print RL "viii) cDNA2orf.WS$ver.gz           -   Latest set of ORF connections to each cDNA (EST, OST, mRNA)\n";
  print RL "ix)   gene_interpolated_map_positions.WS$ver.gz    - Interpolated map positions for each coding/RNA gene\n";
  print RL "x)    clone_interpolated_map_positions.WS$ver.gz   - Interpolated map positions for each clone\n";
  print RL "xi)   best_blastp_hits.WS$ver.gz  - for each C. elegans WormPep protein, lists Best blastp match to
                            human, fly, yeast, C. briggsae, and SwissProt & TrEMBL proteins.\n";

  print RL "xii)  best_blastp_hits_brigprot.WS$ver.gz   - for each C. briggsae protein, lists Best blastp match to
                                     human, fly, yeast, C. elegans, and SwissProt & TrEMBL proteins.\n";  

  print RL "xiii) geneIDs.WS$ver.gz   - list of all current gene identifiers with CGC & molecular names (when known)\n";
  print RL "xiv)  PCR_product2gene.WS$ver.gz   - Mappings between PCR products and overlapping Genes\n";

  print RL "\n\n";
  print RL "Release notes on the web:\n-------------------------\n";
  print RL "http://www.wormbase.org/wiki/index.php/Release_Schedule\n\n\n\n";
  
  # make the chromosomal sequence changes file
  open (CC, "> $reports_dir/chromosome_changes") || die "Can't open file $reports_dir/chromosome_changes\n";
  my @mapping_data = Remap_Sequence_Change::read_mapping_data($ver-1, $ver, $wormbase->species);
  my $text = Remap_Sequence_Change::write_changes($wormbase, $ver, @mapping_data);
  print CC $text;
  close(CC);

  my @release_files = ("$reports_dir/composition","$reports_dir/chromosome_changes","$reports_dir/genedata","$reports_dir/wormpep");
  
  #include all the pre-generated reports
  my $file = shift(@release_files);
  while (defined($file)) {
    open (READIN, "<$file") || die "cant open $file\n";
    while(<READIN>) {
      print RL "$_";
    }
    close READIN;
    print RL "\n\n";
    $file = shift(@release_files);
  }


  # Find out Gene->CDS, Transcript, Pseudogene connections
  my $tace = $wormbase->tace;
  my $db = Ace->connect(-path  => $ace_dir,
                        -program =>$tace) || $log->log_and_die("Connection failure: ",Ace->error);
  my $query = "Find Gene WHERE (Corresponding_CDS OR Corresponding_transcript OR Corresponding_pseudogene) AND CGC_name AND Species = \"*elegans\"";
  my $gene_seq_count = $db->fetch(-query=> "$query");

  # wormpep status overview
  my %wp_status;
  my $wormpep_datafile = "$basedir/WORMPEP/wormpep$ver/wormpep_current";
  
  $wp_status{Confirmed}  = `grep Confirmed    $wormpep_datafile | wc -l`;
  $wp_status{Supported}  = `grep confirmed    $wormpep_datafile | wc -l`;
  $wp_status{Predicted}  = `grep Predicted    $wormpep_datafile | wc -l`;
  $wp_status{Gene}       = $gene_seq_count;
  $wp_status{Uniprot}    = `grep 'UniProt:'        $wormpep_datafile | wc -l`;
  $wp_status{Protein_ID} = `grep 'protein_id' $wormpep_datafile | wc -l`;
  $wp_status{Total}      = $wp_status{Confirmed} + $wp_status{Supported} + $wp_status{Predicted}; 
  
  printf RL "\n\n";
  printf RL "Status of entries: Confidence level of prediction (based on the amount of transcript evidence)\n";
  printf RL "-------------------------------------------------\n";
  printf RL "Confirmed            %6d (%2.1f%%)\tEvery base of every exon has transcription evidence (mRNA, EST etc.)\n", $wp_status{Confirmed}, (($wp_status{Confirmed}/$wp_status{Total}) * 100);
  printf RL "Partially_confirmed  %6d (%2.1f%%)\tSome, but not all exon bases are covered by transcript evidence\n", $wp_status{Supported}, (($wp_status{Supported}/$wp_status{Total}) * 100);
  printf RL "Predicted            %6d (%2.1f%%)\tNo transcriptional evidence at all\n", $wp_status{Predicted}, (($wp_status{Predicted}/$wp_status{Total}) * 100);

  printf RL "\n\n\n";
  printf RL "Status of entries: Protein Accessions\n";
  printf RL "-------------------------------------\n";
  printf RL "UniProtKB accessions %6d (%2.1f%%)\n", $wp_status{Uniprot}, (($wp_status{Uniprot}/$wp_status{Total}) * 100);
  printf RL "\n\n\n";
  printf RL "Status of entries: Protein_ID's in EMBL\n";
  printf RL "---------------------------------------\n";
  printf RL "Protein_id           %6d (%2.1f%%)\n", $wp_status{Protein_ID}, (($wp_status{Protein_ID}/$wp_status{Total}) * 100);
  printf RL "\n\n\n";
  printf RL "Gene <-> CDS,Transcript,Pseudogene connections\n";
  printf RL "----------------------------------------------\n";
  printf RL "Caenorhabditis elegans entries with WormBase-approved Gene name %6d\n", $wp_status{Gene};
  printf RL "\n\n";
  
  ######################################
  #  Get the GO annotations  - Ranjana #
  ######################################

  #Table maker definition
  #SCRIPT:GeneGO_codes.def
  # WBgene00000001 GO:0000001 IEA
  if (-e "$basedir/autoace/wquery/SCRIPT:GeneGO_codes.def") {
    my $def = "$basedir/autoace/wquery/SCRIPT:GeneGO_codes.def";
    my $command = "Table-maker -p $def\nquit\n";
    $log->write_to("Retrieving GO::Gene data, using Table-maker...\n");

    my %Gene2code;
    my %GO2code;

    open (TACE, "echo '$command' | $tace $ace_dir | ") || die "Cannot query acedb. $command  $tace\n";
    while (<TACE>) {
      chomp;
      s/"//g;
      s/GO://g;
      unless (/(WBGene\d+)\s+(\d+)\s+(\w+)/) {next;}
      if (/(WBGene\d+)\s+(\d+)\s+(\w+)/){
	#query with WBGene and GO_code is returned
	$Gene2code{$1} = $3;
	#query with GO term and GO_code is returned
	$GO2code{$2} = $3;
      } 
    }
    close TACE;

    # of genes with GO terms
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

    #and the same for the annotations
    # of GO annotations
    my $query2  = "Find GO_term";
    my $GO_annotations = $db->fetch(-query=> "$query2");
    # GO2Gene
    my $query3  = "Find GO_term WHERE Gene";
    my $GOwithgene = $db->fetch(-query=> "$query3");
    # of IEA GO annotations
    my $IEAno = grep {/IEA/} values %GO2code; 

    # need to know how many GO_terms keys for next $non_IEAno calc.
    my $GOkeys = keys %GO2code;

    # of non-IEA GO annotations
    my $nonIEAno = $GOkeys - $IEAno; 
    
    
    print RL "GO Annotation Stats WS$ver\n----------------------------------------------\n\n";
    print RL "		+----------------------------------------+--------+\n";
    print RL "		| Genes with GO_term connections         | $gc  |\n";
    print RL "		|             From RNAi mappings         | $rnai  |\n";
    print RL "		|                    From citace         |  $citace  |\n";
    print RL "		|                      Inherited         | $inherit  |\n";
    print RL "		|            IEA GO_code present         | $IEA_gc  |\n";
    print RL "		|        non-IEA GO_code present         |  $non_IEA_gc  |\n";
    print RL "		+----------------------------------------+--------+\n";
    print RL "		| Total No. GO_terms                     | $GO_annotations  |\n";
    print RL "		| GO_terms connected to Genes            |  $GOwithgene  |\n";
    print RL "		| GO annotations connected with IEA      |  $IEAno  |\n";
    print RL "		| GO annotations connected with non-IEA  |  $nonIEAno  |\n";
    print RL "		+----------------------------------------+--------+\n";
  }
  else {
    $log->write_to("\nERROR - GeneGO_codes.def abscent from autoace/wquery\nThese stats will be missing from the release letter\n\n");
  }
  # Close the database connection now we have finished with it
  $db->close;


#################################################################################################
  
  # Synchronisation with GenBank / EMBL
  my @chromosomes = ("I","II","III","IV","V","X");
  my $csome = shift @chromosomes;
  print RL "\nSynchronisation with GenBank / EMBL:\n------------------------------------\n\n";
  my $check = 0;
  while ($csome) {
    my $errors = `grep ERROR $ace_dir/yellow_brick_road/CHROMOSOME_$csome.agp_seq.log`;
    while( $errors =~ m/for\s(\p{IsUpper}\w+)/g ) {
      print RL "CHROMOSOME_$csome\tsequence $1\n";
      $check = 1;
    }
    $csome = shift @chromosomes;
  }
  if ($check == 0) {
    print RL "No synchronisation issues\n\n";
  }
  
  print RL "\n";
  
  # Gap summary (hard coded at the moment)
  print RL "There are no gaps remaining in the genome sequence\n";
  print RL "---------------\n";
  print RL "For more info mail worm\@sanger.ac.uk\n";
  print RL "-===================================================================================-\n";
  
  # User filled sections
  print RL "\n\n\n";
  print RL "New Data:\n---------\n\n\n";
  print RL "Genome sequence updates:\n-----------------------\n\n\n";
  print RL "New Fixes:\n----------\n\n\n";
  print RL "Known Problems:\n---------------\n\n\n";
  print RL "Other Changes:\n--------------\n\n";
  print RL "Proposed Changes / Forthcoming Data:\n-------------------------------------\n\n\n";
  print RL "Model Changes:\n------------------------------------\n\n\n";
  
  # Installation guide
  print RL "-===================================================================================-\n";
  print RL "\n\n";
  print RL "Quick installation guide for UNIX/Linux systems\n-----------------------------------------------\n\n";
  print RL "1. Create a new directory to contain your copy of WormBase,\n\te.g. /users/yourname/wormbase\n\n";
  print RL "2. Unpack and untar all of the database.*.tar.gz files into\n\tthis directory. You will need approximately 2-3 Gb of disk space.\n\n";
  print RL "3. Obtain and install a suitable acedb binary for your system\n\t(available from www.acedb.org).\n\n";
  print RL "4. Use the acedb 'xace' program to open your database, e.g.\n\ttype 'xace /users/yourname/wormbase' at the command prompt.\n\n";
  print RL "5. See the acedb website for more information about acedb and\n\tusing xace.\n\n";
  
  
  print RL "____________  END _____________\n";
  close(RL);


  print "DONT FORGET TO FILL IN THE LAST FEW FIELDS IN THE LETTER\n found at $release_letter\n";
  

##################
# Check the files
##################

  $wormbase->check_file($release_letter, $log,
                  minsize => 5500);
}



# say goodbye
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
