#!/software/bin/perl -w
#
# script_template.pl                           
# 
# by Gary Williams                         
#
# This script looks for overlaps of a gene span against the modENCODE
# tiling array Transcriptionally Active Regions (TARs) and works out a
# simple expression level value for the gene. The resulting file of
# expression values is intended to be used by the Caltech expression
# people in SPELL.
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2013-05-14 12:59:47 $      


# new version:
# don't try to pull information out of the acedb database
# do as much as possible from the GFF dumps.
# do Overlap comparison to get overlaps of each TranscriptionallyActiveRegion against Transcript exons.
# (check we can match with the same TAR region against several overlapping exons)
# find the proportion of the TAR that overlaps with the exon and multiply the TAR score by that proportion
# sum the scores for each transcript
# for each gene: divide the sum of the transcript scores by the number of transcripts from that gene











use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
use Modules::Overlap;
#use Sequence_extract;
#use Coords_converter;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $database, $species);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "database:s" => \$database,
	    "species:s"  => \$species,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species,
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

if (!defined $database) {$database = $wormbase->database("autoace")}

my $ovlp = Overlap->new($database, $wormbase);


##########################
# MAIN BODY OF SCRIPT
##########################


# hash of the analysis names keyed by the modENCODE DCC_ID
my %analysis = (

		"modENCODE_476" => "TAR_early_embryo_20dC_0_4hrs_post_fertilization_N2", 
		"modENCODE_3407" => "TAR_EE_Z1_Z4_male",# not in mSTAD data
		"modENCODE_3172" => "TAR_embryo_0hr_reference", 
		"modENCODE_654" => "TAR_embryo_A_class_motor_neurons", 
		"modENCODE_456" => "TAR_embryo_all_cells_reference", 
		"modENCODE_459" => "TAR_embryo_AVA_neurons", # 
		"modENCODE_3173" => "TAR_embryo_AVE_neurons", # 
		"modENCODE_2499" => "TAR_embryo_BAG_neurons", # 
		"modENCODE_470" => "TAR_embryo_body_wall_muscle", # 
		"modENCODE_458" => "TAR_embryo_coelomocytes", # 
		"modENCODE_467" => "TAR_embryo_dopaminergic_neurons", # 
		"modENCODE_468" => "TAR_embryo_GABA_motor_neurons", # 
		"modENCODE_661" => "TAR_embryo_germline_precursor_cells", # 
		"modENCODE_662" => "TAR_embryo_hypodermal_cells", # 
		"modENCODE_457" => "TAR_embryo_intestine", # 
		"modENCODE_455" => "TAR_embryo_panneural", # 
		"modENCODE_2548" => "TAR_embryo_pharyngeal_muscle", # 
		"modENCODE_2500" => "TAR_embryo_PVC_neurons", # not in mSTAD data 
		"modENCODE_3408" => "TAR_emb_Z1_Z4", # not in mSTAD data 
		"modENCODE_481" => "TAR_gonad_from_young_adult_20dC_42hrs_post_L1_N2", # not in mSTAD data 
		"modENCODE_484" => "TAR_L1_20dC_0hrs_post_L1_N2", # 
		"modENCODE_472" => "TAR_L2_25dC_14hrs_post_L1_N2", # 
		"modENCODE_469" => "TAR_L2_A_class_neuron", # 
		"modENCODE_3405" => "TAR_L2_AFD", # not in mSTAD data 
		"modENCODE_465" => "TAR_L2_body_wall_muscle", # 
		"modENCODE_657" => "TAR_L2_coelomocytes", # 
		"modENCODE_464" => "TAR_L2_excretory_cell", # 
		"modENCODE_3406" => "TAR_L2_GABA_alr_1", # not in mSTAD data 
		"modENCODE_466" => "TAR_L2_GABA_neurons", # 
		"modENCODE_658" => "TAR_L2_glutamate_receptor_expressing_neurons", # 
		"modENCODE_463" => "TAR_L2_intestine", # 
		"modENCODE_462" => "TAR_L2_panneural", # 
		"modENCODE_477" => "TAR_L2_polyA_enriched_20dC_14hrs_post_L1_N2", # not in mSTAD data 
		"modENCODE_461" => "TAR_L2_reference__mockIP", # 
		"modENCODE_474" => "TAR_L3_25dC_25hrs_post_L1_N2", # 
		"modENCODE_655" => "TAR_L3_L4_dopaminergic_neuron", # 
		"modENCODE_2454" => "TAR_L3_L4_hypodermal_cells", # 
		"modENCODE_460" => "TAR_L3_L4_PVD___OLL_neurons", # 
		"modENCODE_3440" => "TAR_L3_L4_rectal_epithelial_cells", # not in mSTAD data 
		"modENCODE_659" => "TAR_L3_L4_reference__mockIP", # 
		"modENCODE_473" => "TAR_L4_25dC_36hrs_post_L1_N2", # not in mSTAD data 
		"modENCODE_479" => "TAR_late_embryo_20dC_6_12hrs_post_fertilization_N2", # 
		"modENCODE_478" => "TAR_male_L4_25dC_36hrs_post_L1_CB4689", # 
		"modENCODE_485" => "TAR_soma_only_mid_L4_25dC_36hrs_post_L1_JK1107", # not in mSTAD data 
		"modENCODE_475" => "TAR_young_adult_25dC_42hrs_post_L1_N2", # 
		"modENCODE_660" => "TAR_Young_Adult_Cephalic_sheath__CEPsh", # 
		"modENCODE_656" => "TAR_Young_Adult_reference__mockIP", # 
		"modENCODE_491" => "TAR_pathogen_control_OP50_25dC_24hr_exposure_post_adulthood_N2", # not in mSTAD data 
		"modENCODE_490" => "TAR_pathogen_control_OP50_25dC_48hr_exposure_post_adulthood_N2", # not in mSTAD data 
		"modENCODE_487" => "TAR_pathogen_Efaecalis_25dC_24hr_exposure_post_adulthood_N2", # not in mSTAD data 
		"modENCODE_486" => "TAR_pathogen_Pluminscens_25dC_24hr_exposure_post_adulthood_N2", # not in mSTAD data 
		"modENCODE_489" => "TAR_pathogen_Smarcescens_25dC_24hr_exposure_post_adulthood_N2", # not in mSTAD data 
		"modENCODE_488" => "TAR_pathogen_Smarcescens_25dC_48hr_exposure_post_adulthood_N2", # not in mSTAD data 
	       );



# output directory
my $outdir = $wormbase->autoace."/TARS";
mkdir $outdir, 0777;
if (-e glob("$outdir/*.out")) {
  $wormbase->run_command("rm $outdir/*.out", $log);
}

my @chroms = $wormbase->get_chromosome_names(-prefix => 1, mito => 0);


foreach my $chrom (@chroms) {
  gene_expression($chrom);
}

$log->write_to("Making final expr.tar.gz file ...\n");
chdir ($outdir) || $log->log_and_die("Couldn't chdir to $outdir\n");
unlink ("expr.tar") || $log->log_and_die("ERROR: Cannot delete file expr.tar :\t$!\n");;
my $status = $wormbase->run_command("tar cf expr.tar *.out", $log);
$log->write_to("status of tar command: $status\n");

unlink "expr.tar.gz";
$status = $wormbase->run_command("gzip -f expr.tar", $log);
$log->write_to("status of gzip command: $status\n");


$log->mail();
print "Finished.\n" if ($verbose);

exit(0);





##############################################################
#
# Subroutines
#
##############################################################


sub gene_expression {
  my ($chromosome) = @_;
  
  print "Sequence: ", $chromosome, "\n";

  my $common_dir      = "$database/COMMON_DATA";       # GFF directory  
  my %dummy;
  my %worm_gene2geneID_name = $wormbase->FetchData('worm_gene2geneID_name', \%dummy, $common_dir);

  my @exons = $ovlp->get_Coding_transcript_exons($chromosome);
  
  # get the TARs
  my %GFF_data = 
    (
     method			=> "TranscriptionallyActiveRegion", # for looking for a GFF_SPLITS file which won't exist
     gff_source			=> "TranscriptionallyActiveRegion",
     gff_type			=> "transcribed_fragment",
     ID_after			=> "Note\\s+",
     score                      => 1,
   );
#CHROMOSOME_III  TranscriptionallyActiveRegion   transcribed_fragment    21292   21449   8.95931 .       .       Note "modENCODE_3173"

  my @TARs = $ovlp->read_GFF_file($chromosome, \%GFF_data);


  #print @genes_list;
  #print @tars_list;

  # OK - we have the genes and the TARs in lists
  # now look for overlaps

  my $exons_overlap = $ovlp->compare(\@exons);

  my %transcript_score;
  my %transcript_len;

  foreach my $tar (@TARs) { 
    my $tar_id = $tar->[0];
    #print "next TAR to check = $tar_id\n";
    my @exon_matches = $exons_overlap->match($tar);
    #print "Have ", scalar @exon_matches, " overlaps to TAR\n";
    my @exon_ids = $exons_overlap->matching_IDs;
    #print "Matching exon IDs: @exon_ids\n";
    foreach my $exon_match (@exon_matches) {

      my ($proportion1, $proportion2) = $exons_overlap->matching_proportions($exon_match);

      my $score = $tar->[5];

      my $transcript_name = $exon_match->[0];

      my $current_score = $transcript_score{$transcript_name}{$tar_id};

      # get the proportion of the exon that is overlapped by the TAR
      # if this is more than 50% then store the highest score for this transcript and TAR
      if ($proportion2 > 0.5) {
	if (!defined $current_score || $current_score <  $score) {
	  $transcript_score{$transcript_name}{$tar_id} = $score;
	}
      }
    }
  }




  # now write out the transcript expression values 
  foreach my $transcript_name (keys %transcript_score) {
    foreach my $tar_id (keys %{$transcript_score{$transcript_name}}) {

      my $file = $analysis{$tar_id};
      if (!defined $file) {$log->log_and_die("Can't find an analysis object name for the analysis from: $tar_id\n");}

      my $transcript_score = $transcript_score{$transcript_name}{$tar_id};
      if (!defined $transcript_score || $transcript_score == 0) {$transcript_score = 0.0000000001} # Wen wants this small value used when the value is zero

      open (TRANS_TAR, ">> $outdir/transcripts_${file}.out") || $log->log_and_die("Can't open file $outdir/transcripts_${file}.out\n");
      print TRANS_TAR "$transcript_name\t$transcript_score\n";
      close(TRANS_TAR);
    }
  }




  # now get the highest score for the genes for each TAR experiment

  my %gene_score;
  my $cds_regex = $wormbase->cds_regex_noend;

  foreach my $transcript_name (keys %transcript_score) {
    foreach my $tar_id (keys %{$transcript_score{$transcript_name}}) {
      my $transcript_score = $transcript_score{$transcript_name}{$tar_id};

      my ($gene_name) = ($transcript_name =~ /($cds_regex)/);
      if (!defined $gene_name) {$gene_name = $transcript_name}# the tRNAs (e.g. C06G1.t2) are not changed correctly

      my $gene_score = $gene_score{$gene_name}{$tar_id};

      if (!defined $gene_score || $gene_score <  $transcript_score) {
	$gene_score{$gene_name}{$tar_id} = $transcript_score;
      }
    }
  }


  # then output the best scores for each gene

  foreach my $gene_name (keys %gene_score) {
    my $gene_id = $worm_gene2geneID_name{$gene_name};
    if (!defined $gene_id) {$log->log_and_die ("Can't find gene_id $gene_name\n")}

    foreach my $tar_id (keys %{$gene_score{$gene_name}}) {
      my $gene_score = $gene_score{$gene_name}{$tar_id};
      if ($gene_score == 0) {$gene_score = 0.0000000001} # Wen wants this small value used when the value is zero

      # now for each TAR, output the gene scores and translating
      # the modENCODE DCC_ID number to an analysis name

      my $file = $analysis{$tar_id};
      if (!defined $file) {$log->log_and_die("Can't find an analysis object name for the analysis from: $tar_id\n");}

      open (GENE_TAR, ">> $outdir/genes_${file}.out") || $log->log_and_die("Can't open file $outdir/genes_${file}.out\n");
      print GENE_TAR "$gene_id\t$gene_score\n";
      close(GENE_TAR);
    }
  }

}


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

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

=back

This script does...blah blah blah

script_template.pl MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

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

=item xxx (xxx@sanger.ac.uk)

=back

=cut
