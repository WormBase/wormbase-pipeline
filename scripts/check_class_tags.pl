#!/software/bin/perl -w
# check_class_tags.pl
#
# by Gary Williams
#
# Get the counts of various tags in selected objects in various classes.
# Compare these counts to those from previous releases at this stage in the Build.
#

use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Ace;
use File::Compare;
use List::MoreUtils qw(uniq);

my ($wormbase, $verbose, $store, $help, $debug, $test, $species, $database);

GetOptions(
	   "help"       => \$help,
	   "debug=s"    => \$debug,
	   'store=s'    => \$store,
	   'species:s'  => \$species,
	   'database:s' => \$database, # used for loading the data file with previous versions when setting up for the first time
	   'test'       => \$test,
	   'verbose:s'  => \$verbose,
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


############################################
# Get paths and versions
############################################

my $tace = $wormbase->tace;


$species = $wormbase->species;
my $version = $wormbase->get_wormbase_version;
my $prev_version = $version-1;
my $prev_prev_version = $version-2;


$database ||= $wormbase->autoace;

my $dbname_0    = "WS${prev_prev_version}";
my $dbname_1    = "WS${prev_version}";
my $dbname_2    = "WS${version}";

my $count = 0;
my $errors = 0;


my $file = $wormbase->build_data . "/COMPARE/${species}_class_tags_count.dat"; # file holding pparse details from previous Builds


$log->write_to("Checking ".$wormbase->full_name.": - ".$database."\n\n");


my %check = &get_classes_and_ids_to_check();

&print_header(%check);

&check_classes_and_ids(%check);

$log->write_to("\n$count class counts checked, $errors potential errors found\n");

$log->mail;
exit;




##############################################################
#
# Subroutines
#
##############################################################

sub usage {
    my $error = shift;

    if ( $error eq "Help" ) {

        # Normal help menu
        system( 'perldoc', $0 );
        exit(0);
    }
}

##################################################################
# these are the classes that we wish to check for differences in the numbers of tags
# the example IDs are pretty much chosen at random - if you find a better one to check, feel free to change it.

sub get_classes_and_ids_to_check {


  my %classes = (
		 'elegans' => {
			       '2_point_data' => '190',
			       Accession_number => '1FHO_A',
			       Ace2SO => 'transposable_element_ace2so',
			       Analysis => 'Million_Mutation_Project',
			       Anatomy_function => 'WBbtf0550',
			       Anatomy_name => 'AB.prpaaapa',
			       Anatomy_term => 'WBbt:0001554',
			       Antibody => '[cgc4594]:snt-1',
			       Author => 'Accili D',
			       CDS => 'AC8.12',
			       Cell => 'ABprpaapap',
			       Cell_group => 'hypodermis',
			       Clone => 'C03B5',
			       Condition => 'cgc4489_ref',
			       Contig => 'ctg899',
			       Database => 'AceView',
			       Expression_cluster => 'WBPaper00024671:AFD_vs_AWB_downregulated',
			       Expression_pattern => 'Chronogram273',
			       Expr_profile => 'B0213.14',
			       Feature => 'WBsf000351',
			       Feature_data => 'yk1244g06.5:polyA',
			       Gene => 'WBGene00000273',
			       Genetic_code => 'Selenocysteine',
			       Gene_class => 'eva',
			       Gene_cluster => 'rDNA_cluster',
			       Gene_name => '6R55.2',
			       Gene_regulation => 'cgc6998_egl-15',
			       GO_code => 'TAS',
			       GO_term => 'GO:0000351',
			       Grid => 'Y87-96',
			       Homology_group => 'COG0471',
			       Homol_data => ['AC3:Expr','AC3:Mass-spec','AC3:RepeatMasker','AC3:RNAi','AC3:SAGE','AC3:TEC_RED','AC3:wublastx_brenneri','AC3:wublastx_briggsae','AC3:wublastx_brugia','AC3:wublastx_fly','AC3:wublastx_human','AC3:wublastx_japonica','AC3:wublastx_ovolvulus','AC3:wublastx_pristionchus','AC3:wublastx_remanei','AC3:wublastx_slimSwissProt','AC3:wublastx_worm','AC3:wublastx_yeast','AC3:wublastx_sratti',
'Oligo_set:CHROMOSOME_V:WUSTL_1', 'Oligo_set:CHROMOSOME_V:UCSF', 'Oligo_set:CHROMOSOME_V:NIMBLEGEN_1', 'Oligo_set:CHROMOSOME_V:AGILENT_1', 'Oligo_set:CHROMOSOME_V:AFFY_1', 'CHROMOSOME_V:RNAi', 'BLAT_mRNA:CHROMOSOME_V', 'BLAT_mRNA:CHROMOSOME_V_8', 'BLAT_EST:CHROMOSOME_V', 'BLAT_EST:CHROMOSOME_V_1', 'BLAT_Trinity:CHROMOSOME_V', 'BLAT_Trinity:CHROMOSOME_V_1', 'BLAT_OST:CHROMOSOME_V', 'BLAT_OST:CHROMOSOME_V_1', 'BLAT_tc1:CHROMOSOME_V', 'BLAT_tc1:CHROMOSOME_V_1', 'BLAT_RST:CHROMOSOME_V', 'BLAT_RST:CHROMOSOME_V_1', 'CHROMOSOME_V:SAGE', 'CHROMOSOME_V:Mass-spec', 'CHROMOSOME_V:Expr', 'CHROMOSOME_V:wublastx_brenneri', 'CHROMOSOME_V:wublastx_briggsae', 'CHROMOSOME_V:wublastx_brugia', 'CHROMOSOME_V:wublastx_fly', 'CHROMOSOME_V:wublastx_human', 'CHROMOSOME_V:wublastx_japonica', 'CHROMOSOME_V:wublastx_ovolvulus', 'CHROMOSOME_V:wublastx_pristionchus', 'CHROMOSOME_V:wublastx_remanei', 'CHROMOSOME_V:wublastx_slimSwissProt', 'CHROMOSOME_V:wublastx_sratti', 'CHROMOSOME_V:wublastx_worm', 'CHROMOSOME_V:wublastx_yeast'],
			       Interaction => 'WBInteraction000000162',
			       Laboratory => 'RW',
			       Library => 'Vancouver_fosmid',
			       Life_stage => 'WBls:0000072',
			       Locus => 'syP8',
			       Map => 'X',
			       Mass_spec_experiment => 'Zhu_1',
			       Mass_spec_peptide => 'MSP:AAAEEYPVDIVDLSDDFK',
			       Method => 'miRNA_primary_transcript',
			       Microarray => 'WashU_GSC_C.elegans_Genome_Array',
			       Microarray_experiment => 'WBPaper00013462:14_days_N2_5',
			       Microarray_results => '172031_x_at',
			       Molecule => 'WBMol:00000194',
			       Motif => 'BlnI',
			       Movie => '012.C12.i2.z7.mov',
			       Multi_pt_data => '913',
			       Oligo => 'cenix:11-b10_T7',
			       Oligo_set => '172031_x_at',
			       Operon => 'CEOP2352',
			       Paper => 'WBPaper00000277',
			       PCR_product => 'cenix:12-h8',
			       Person => 'WBPerson10000',
			       Person_name => 'A Gottschalk',
			       Phenotype => 'WBPhenotype:0000195',
			       Picture => '295_BC10719.png',
			       Position_Matrix => 'WBPmat00000273',
			       Pos_neg_data => '1417',
			       Protein => 'ENSEMBL:ENSMUSP00000042619',
			       Protein => 'WP:CE10000',
			       Pseudogene => 'C44C10.2',
			       Rearrangement => 'meDf5',
			       RNAi => 'WBRNAi00000273',
			       SAGE_tag => 'SAGE:aaaaaaaaatccacgtt',
			       Sequence => ['yk786f06.5', 'C25A1', 'F56A3', 'C04H5', 'B0432', 'C07A9', 'F30H5', 'C10C6', 'B0545', 'C12D8', 'K04F1', 'C02C6', 'AH9'],					    
			       Sequence_collection => 'Genome:C_elegans-WBcel235',
			       SK_map => 'AH10.1:Sep2001',
			       SO_term => 'SO:0000458',
			       Species => 'Achromobacter cycloclastes',
			       Strain => 'BC2420',
			       Structure_data => 'WBStructure000191',
			       Transcript => 'B0205.9',
			       Transcription_factor => 'WBTranscriptionFactor000119',
			       Transgene => 'eIs2137',
			       Transposon => 'WBTransposon00000195',
			       Transposon_family => 'TURMOIL2',
			       Tree => 'Z4.aaa post-emb lineage male vers 2',
			       TreeNode => 'ABarppppp',
			       Variation => 'WBVar00000273',
			      },
		 'briggsae' => {
			       Ace2SO => 'transposable_element_ace2so',
			       Analysis => 'RNASeq.briggsae.AF16.WBls:0000038.Hermaphrodite.WBbt:0007833.PRJNA75295.SRX052081',
			       CDS => 'CBG00033',
			       Clone => 'CBG29121',
			       Condition => 'RNASeq.briggsae.L4_larva',
			       Feature => 'WBsf028129',
			       Feature_data => 'AF520619:TSL',
			       Gene => 'WBGene00086998', # misses Ortholog Ortholog_other Other_name
			       Gene_name => 'Cbr-glb-26',
			       Homol_data => 'cb25.fpc2220b:wublastx_slimSwissProt',
			       Method => 'cDNA_for_RNAi',
			       Protein => 'ENSEMBL:ENSMUSP00000042619',
			       #Pseudogene => '',
			       Sequence => ['cb25.fpc0002', 'cb25.fpc0011c', 'cb25.fpc0081', 'cb25.fpc0143a', 'chrI', 'chrII'],
			       Transcript => 'CBG00122',
			       Variation => 'WBVar00000752',

			       },
		 'brenneri' => {
			       Ace2SO => 'transposable_element_ace2so',
			       Analysis => 'RNASeq.brenneri.DF5081.WBls:0000038.Unknown.WBbt:0007833.PRJNA75295.SRX100771',
			       CDS => 'CBN00033',
			       Clone => 'CBN32995',
			       Condition => 'RNASeq.brenneri.L4_larva.Replicate2',
			       Feature_data => 'Cbre_Contig0:RNASeq_forward_reads',
			       Gene => 'WBGene00158496', # misses Ortholog Ortholog_other Other_name
			       Gene_name => 'Cbn-acly-2',
			       Homol_data => 'BLAT_EST:Cbre_Contig0_18',
			       Method => 'cDNA_for_RNAi',
			       Pseudogene => 'CBN09775',
			       Sequence => ['Cbre_Contig1', 'Cbre_Contig10', 'Cbre_Contig20', 'Cbre_Contig50', 'Cbre_Contig100', 'Cbre_Contig200', 'Cbre_Contig400', 'Cbre_Contig600', 'Cbre_Contig800'],
			       Transcript => 'CBN00079',
			       },
		 'remanei' => {
			       Ace2SO => 'transposable_element_ace2so',
			       Analysis => 'RNASeq.remanei.SB146.WBls:0000027.Unknown.WBbt:0007833.PRJNA75295.SRX101885',
			       CDS => 'CRE00076',
			       Clone => 'CRE32638',
			       Condition => 'RNASeq.remanei.L4_larva.Replicate2',
			       Feature_data => 'AY589598:TSL',
			       Gene => 'WBGene00051012', # misses Ortholog Ortholog_other Other_name
			       Gene_name => 'Cre-acd-1',
			       Sequence => ['Crem_Contig0', 'Crem_Contig10', 'Crem_Contig15', 'Crem_Contig30', 'Crem_Contig100', 'Crem_Contig200', 'Crem_Contig300', 'Crem_Contig500', 'Crem_Contig800'],
			      },
		 'japonica' => {
			       Ace2SO => 'transposable_element_ace2so',
			       Analysis => 'RNASeq.japonica.DF5081.WBls:0000038.Unknown.WBbt:0007833.PRJNA75295.SRX100091',
			       CDS => 'CJA00088',
			       Condition => 'RNASeq_Hillier.japonica.L4_larva_Replicate',
			       Feature_data => 'CA758779:low',
			       Gene => 'WBGene00119208', # misses Ortholog Ortholog_other Other_name
			       Gene_name => 'Cjp-acd-5',
			       Sequence => ['Cjap.Contig2', 'Cjap.Contig7087', 'Cjap.Contig772', 'Cjap.Contig7751', 'Cjap.Contig8547', 'Cjap.Contig91', 'Cjap.Contig9449', 'Cjap.Contig981', 'Cjap.Contig9854'],
			       Transcript => 'CJA00088',
			       },
		 'pristionchus' => {
			       CDS => 'PPA00099',
			       Feature_data => 'AA191935:polyA_site',
			       Gene => 'WBGene00089610', # misses Ortholog Ortholog_other Other_name
			       Gene_name => 'Ppa-abcf-1',
			       Sequence => ['Ppa_Contig0', 'Ppa_Contig10', 'Ppa_Contig15', 'Ppa_Contig30', 'Ppa_Contig100', 'Ppa_Contig200'],
				   },
		 'brugia' => {
			       Ace2SO => 'coding_transcript_ace2so',
			       CDS => 'Bm7483',
			       Condition => 'RNASeq.brugia.ERP000948.adult_female',
			       Feature => 'WBsf899556',
			       Feature_data => 'Bm_v4_Chr4_scaffold_001:Dust',
			       Gene => 'WBGene00220262',
			       Gene_name => 'Bma-aagr-4',
			       Homol_data => 'Bm_v4_Chr4_scaffold_001:wublastx_brenneri',
			       Method => 'BLAT_EST_BEST',
			       Protein => 'BM:BM32546',
			       Pseudogene => 'Bm477',
			       Sequence => ['Bm_024', 'Bm_013', 'Bm_007', 'Bm_008', 'Bm_014', 'Bm_v4_Chr4_scaffold_001'],
			       Transcript => 'Bm1',
			     },
		 'ovolvulus' => {
			       Ace2SO => 'coding_transcript_ace2so',
			       CDS => 'OVOC9637',
			       Feature_data => 'OVOC_OO_000024:TRF',
			       Gene => 'WBGene00246446',
			       Gene_name => 'Ovo-eat-4',
			       Homol_data => 'OVOC_OO_000024:wublastx_brenneri',
			       Method => 'BLAT_EST_BEST',
			       Protein => 'OV:OVP01471',
			       Sequence => ['OVOC_OO_000001', 'OVOC_OO_000008', 'OVOC_OO_000054', 'OVOC_OO_000132', 'OVOC_OO_000629', 'OVOC_OO_000690'],
			       Transcript => 'OVOC8637',
			     },
                  'sratti'  => { 
			       Ace2SO => 'coding_transcript_ace2so',
			       CDS => 'SRAE_0000000800',
			       Feature_data => 'SRAE_scaffold23:TRF',
			       Gene => 'XXX',
			       Gene_name => 'SRAE_0000000800',
			       Homol_data => 'SRAE_0000000800:wublastx_brenneri',
			       Method => 'BLAT_EST_BEST',
			       Protein => 'SRP:SRP01471',
			       Sequence => ['SRAE_chr2', 'SRAE_scaffold1'],
			       Transcript => 'SRAE_0000000800',
                             },
		);

  
  return %{$classes{$species}};
}


##################################################################
# print the header

sub print_header {

  my (%check) = @_; # key=class, value=id


  my %results_0; 			# results from the last but one release
  my %results_1;			# results from the last release
  my $got_prev_results=0;
  my $got_prev_prev_results=0;
  for (my $version_count = 10; $version_count; $version_count--) { # go back up to 10 releases
    %results_0 = ();
    if (!$got_prev_results) {%results_1 = ()}
    foreach my $class (keys %check) {
      my $id = $check{$class};
      %results_0 = &get_prev_count($species, $prev_prev_version, $class, $id);
      if (!$got_prev_results) {%results_1 = &get_prev_count($species, $prev_version, $class, $id)}
      if (keys %results_0) {last}
    }
    
    # check to see if we have results from the previous release number we are currently checking
    # (elegans is the only species done every release)
    if (keys %results_1) {
      $got_prev_results = 1;
    }
    if (keys %results_0) {
	$got_prev_prev_results = 1;
    }
    
    if ($got_prev_prev_results) {last}
    
    # go back another version
    $prev_prev_version--;
    $dbname_0    = "WS${prev_prev_version}";
    if (!$got_prev_results) {
      $prev_version--;
      $dbname_1    = "WS${prev_version}";
    }
  }
  
  # display message if there are no previous results
  if (!$got_prev_results) {
    $log->write_to("\n\nNo results have been found for this species for the last 10 releases\n\n");
    $dbname_0 = '';
    $dbname_1 = '';
  }
  
  # display header
  $log->write_to(sprintf("%-22s %-22s %-22s %7s %7s %7s %10s\n", "CLASS","ID", "TAG","($dbname_0)",$dbname_1,$dbname_2,"Difference"));
}
##################################################################

sub check_classes_and_ids {

  my (%check) = @_;

  foreach my $class (keys %check) {
    my $id = $check{$class};
    # I've been lazy and hold the Sequence clone IDs as a array-ref, and the IDs of the other classes as scalars. So check which this is.
    if (ref($id)) {
      foreach my $i (sort @{$id}) {
	print "\nChecking tags in $class : \"$id\"\n" if ($verbose);
	&check_for_missing_tags($class, $i);
      }
    } else {
      print "\nChecking tags in $class : \"$id\"\n" if ($verbose);
      &check_for_missing_tags($class, $id);
    }
  }
}

##################################################################
# check for differences in the numbers of all tags in an object
# the tace command show -a gives a output that is easy to parse, like:

#Sequence : "AC3"
#DNA	 "AC3" 38951
#Gene_child	 "WBGene00195634" 30835 30983
#Gene_child	 "WBGene00199275" 25396 25250
#CDS_child	 "AC3.1:wp98" 5575 7013
#CDS_child	 "AC3.3" 18739 17380
#Transcript	 "AC3.14" 30835 30983
#Transcript	 "AC3.15" 25396 25250
#Transcript	 "AC3.16" 17679 17828
#Pseudogene	 "AC3.9" 11805 13838
#Pseudogene	 "AC3.13" 17199 16871
#Pseudogene	 "AC3.12:wp227" 16245 15940
#Genomic_non_canonical	 "WRM0623cH01" 1614 36865
#PCR_product	 "cenix:116-g2" 9726 10433
#PCR_product	 "cenix:137-a1" 27145 28341
#PCR_product	 "cenix:18-a1" 30698 31254
#Allele	 "WBVar00001390" 615 615
#Allele	 "WBVar00001395" 3025 3025
#Allele	 "WBVar00001400" 4819 4819
#Oligo_set	 "Aff_AC3.1" 6262 7013
#Oligo_set	 "Aff_AC3.2" 12439 13179
#Oligo_set	 "Aff_AC3.3" 17979 17380
#Feature_object	 "WBsf016834" 25468 25469
#Feature_object	 "WBsf017331" 30436 30441
#Feature_object	 "WBsf017332" 30455 30456
#Homol_data	 "AC3:RNAi" 1 38951
#Homol_data	 "AC3:SAGE" 1 38951
#Source	 "CHROMOSOME_V"
#Overlap_right	 "F15H10" 38848
#Overlap_left	 "K07C5"
#Clone_left_end	 "AC3" 1
#Clone_left_end	 "F15H10" 24902
#Clone_right_end	 "K07C5" 5515
#Clone_right_end	 "AC3" 38951
#Database	 "EMBL" "NDB_AC" "Z71177"
#Database	 "EMBL" "NDB_SV" "Z71177.3"
#DB_remark	 "[121025] Sequence correction: SNP 0 bases  @ 34693"
#Keyword	 "HTG"
#EMBL_dump_info	 EMBL_dump_method "worm_EMBL-dump"
#From_author	 "McMurray AA"
#From_laboratory	 "HX"
#Date_directory	 "030414"
#Species	 "Caenorhabditis elegans"
#Strain	 "N2"
#Clone	 "AC3"
#Remark	 "[041026 pad] Genome sequencing error was corrected, removed a single G from position 29071 within the clone." Paper_evidence "WBPaper00024276"
#Remark	 "[041026 pad] Genome sequencing error was corrected, removed a single G from position 29071 within the clone." Accession_evidence "NDB" "BJ109865"
#Genomic_canonical	
#MD5	 "e4ead5658016c4defc883e270a20638d"
#Finished	 1995-12-21


sub check_for_missing_tags {

  my ($class, $id) = @_;

  my @db_slurp = get_tace($class, $id, $database);
  @db_slurp = grep {!/^(acedb\>|\/\/)$/} @db_slurp; # remove tace stuff from output
  map {$_ =~ s/^(\S+).*/$1/} @db_slurp; # replace each element of the array with the first word of each element
  map {$_ =~ s/\n//} @db_slurp; # remove newline


  my %db_count;
  
  foreach my $tag (@db_slurp) {$db_count{$tag}++} # count the unique tags

  my %prev_prev_count = &get_prev_count($species, $prev_prev_version, $class, $id);
  my %prev_count      = &get_prev_count($species, $prev_version,      $class, $id);

  
  ##################################################
  # Calculate difference between databases         #
  ##################################################

  my @tags = uniq (@db_slurp, keys %prev_count, keys %prev_prev_count); # make sure we don't miss any tags (see use List::MoreUtils qw(uniq);)

  foreach my $tag (@tags) {
    my $count0 = $prev_prev_count{$tag} || 0;
    my $count1 = $prev_count{$tag} || 0;
    my $count2 = $db_count{$tag} || 0;
    &store_count($species, $version, $class, $id, $tag, $count2);
    my $diff = $count2 - $count1;
    my $err = "";
    if ($count2 == 0) {
      $err = "***** POSSIBLE ERROR *****";
      $log->error;
      $errors++;
    } elsif (
	     ($count2 < $count1 * 0.9 || 
	      $count2 > $count1 * 1.1)) {
      $err = "***** POSSIBLE ERROR *****";
      $log->error;
      $errors++;
    }
    $count++;
    
    $log->write_to(sprintf("%-22s %-22s %-22s %7d %7d %7d %10d %s\n", $class, $id, $tag, $count0, $count1, $count2, $diff, $err)) if ($err || $verbose);
  }

}

##################################################################
# open tace connection to get the object and slurp up the contents

sub get_tace {
  my ($class, $id, $database) = @_;
  my @slurp;

  my $cmd = "find $class $id\nshow -a\nquit\n";
  open (TACE, "echo '$cmd' | $tace $database |");
  while (<TACE>) {
    chomp;
    next if ($_ eq "");
    next if (/^acedb\>/);
    next if (/^\/\//);
    push @slurp, $_;
  }
  close TACE;

  return @slurp;
}

##################################################################
# get the count of this class found in the previous Build
# retuens a hash of the numbers of each tag in this object
sub get_prev_count {
  my ($species, $version, $class, $id) = @_;

  my %results;

  my $last_count = -1;		# -1 is a flag value indicating no results were found
  if (open (TAG_COUNT, "< $file")) {
    while (my $line = <TAG_COUNT>) {
      chomp $line;
      my ($cc_species, $cc_version, $cc_class, $cc_id, $cc_tag, $cc_count) = split /\t+/, $line;
      # we don't want to get the count from any details that may have
      # been stored by an earlier run of this script in this Build,
      # but we want the most recent version's count apart from that,
      # so get the most recent result that isn't from this Build.
      if ($cc_version == $version && $cc_class eq $class && $cc_species eq $species && $cc_id eq $id) {
	# store to get the last one in the previous build
	$results{$cc_tag} = $cc_count;
      } 
    }
    close (TAG_COUNT);
  }
  return %results;
}

##################################################################
# now store the details for this Build
sub store_count {

  my ($species, $version, $class, $id, $tag, $count) = @_;


  if (open (TAG_COUNT, ">> $file")) {
    if ($version && $class && $species && $count && $tag) {
      print TAG_COUNT "$species\t$version\t$class\t$id\t$tag\t$count\n";
    } else {
      if (!$count) {
	$log->write_to("\n*** ERROR: There are zero $class=$id tag $tag in the database!\n\n");
	$log->error;
      } else {
      $log->log_and_die("*** ERROR: Couldn't write to $file because some of the following is blank\nspecies=$species, version=$version, class=$class, id=$id, tag=$tag, count=$count\n\n");
    }
      $log->error;
    }
    close (TAG_COUNT);
  } else {
    $log->write_to("WARNING: Couldn't write to $file\n\n");
  }
}




##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################

__END__







###############################################

=pod

=head2   NAME - check_class_tags.pl

=head1 USAGE

=over 4

=item check_class_tags.pl [-options]

=back

Get the counts of various tags in selected objects in various classes.
Compare these counts to those from previous releases at this stage in the Build.

1) The following 12 clones are representative of the whole genome in that
they include one Sanger and one St. Louis clone for each chromosome.  Check
each clone to ensure that it contains BLAT data (EST and mRNA), BLAST data,
waba data, gene models, UTRs etc.  Also check for presence of tandem and inverted
repeats which have gone missing in the past

i)    C25A1
ii)   F56A3
iii)  C04H5
iv)   B0432
v)    C07A9
vi)   F30H5
vii)  C10C6
viii) B0545
ix)   C12D8
x)    K04F1
xi)   C02C6
xii)  AH9

2) If the genome sequence has changed, you should inspect the clones containing
those changes to see if there are any strange errors (e.g. duplicate sets of data
which are slightly out of sync.)

3) Check ~wormpub/BUILD/autoace/CHROMOSOMES/composition.all - are there any non-ATCGN
characters

4a) Check that the latest WormPep proteins have proper protein and motif homologies
This has been a problem in some builds where all new WormPep proteins have not got any
BLAST analyses.  Pick a few random Wormpep proteins and especially check that all of
the various blastp homologies are there (human, fly, worm, yeast etc.) and try to
check at least one protein from the ~wormpub/BUILD/WORMPEP/wormpepXXX/new_entries.WSXXX file

4b) Now that we have a curated set of brigpep, should do this periodically for
C. briggase protein objects too...these now have their own set of blastp hits

5) Check PFAM Motif objects have a title tag. It is a problem if there are more than about 20.

6) Run: 
  ls ~wormpub/BUILD/autoace/CHROMOSOMES/*.dna | grep -v masked |grep -v Mt| xargs composition
Make sure this is the same as it was at the start of the build:
  cat ~wormpub/BUILD/autoace/CHROMOSOMES/composition.all
Bad Homol objects can lead to errors esp when chromosome length has been reduced

Thats all...for now!  If you are satisfied the build is ok, please inform the person
building the database. Please continue to add to this list as appropriate.


=over 4

=item MANDATORY arguments: none

=back

=over 4

=item OPTIONAL arguments: -debug, -help, -database


-debug and -help are standard Wormbase script options.

-species specifies the species BUILD database to use.
=back

=head1 AUTHOR - Gary Williams



=cut

