#!/use/bin/env perl
#
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


my $log = Log_files->make_build_log($wormbase);
my $tace = $wormbase->tace;


$species = $wormbase->species;
my $version = $wormbase->get_wormbase_version;

$database ||= $wormbase->autoace;

my $count = 0;
my $errors = 0;

my $file = $wormbase->build_data . "/COMPARE/${species}_class_tags_count.dat"; # file holding pparse details from previous Builds
$log->write_to("Checking ".$wormbase->full_name.": - ".$database."\n\n");

my $DATA = {};

&read_counts($file);
my ($prev_version, $prev_prev_version) = &get_previous_two_versions();

my $dbname_0    = "WS${prev_prev_version}";
my $dbname_1    = "WS${prev_version}";
my $dbname_2    = "WS${version}";

my $check = &get_classes_and_ids_to_check();

$log->write_to(sprintf("%-22s %-22s %-22s %7s %7s %7s %10s\n", "CLASS","ID", "TAG","($dbname_0)",$dbname_1,$dbname_2,"Difference"));

&check_classes_and_ids($check);
$log->write_to("\n$count class counts checked, $errors potential errors found\n");

&store_counts($file);

$log->mail;
exit;

    
##############################################################
##############################################################
sub get_previous_two_versions {
  
  my ($prev_v, $prev_prev_v);

  for(my $v = $version-1; $v >= $version - 10; $v--) {
    if (exists $DATA->{$species}->{$v}) {
      if (not defined $prev_v) {
        $prev_v = $v;
      } elsif (not defined $prev_prev_v) {
        $prev_prev_v = $v;
      }
    }
  }

  $prev_v = 0 if not defined $prev_v;
  $prev_prev_v = 0 if not defined $prev_prev_v;
  
  return ($prev_v, $prev_prev_v);
}



#####################################################3
sub check_classes_and_ids {
  my ($check) = @_;

  delete $DATA->{$species}->{$version} if exists $DATA->{$species}->{$version};

  foreach my $class (keys %$check) {
    my $id = $check->{$class};
    # I've been lazy and hold the Sequence clone IDs as a array-ref, and the IDs of the other classes as scalars. So check which this is.
    if (ref($id)) {
      foreach my $i (sort @{$id}) {
	$log->write_to("Checking tags in $class : \"$i\"\n") if $debug;
	&check_for_missing_tags($class, $i);
      }
    } else {
      $log->write_to("Checking tags in $class : \"$id\"\n") if $debug;
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
####################################################################
sub check_for_missing_tags {

  my ($class, $id) = @_;

  my @db_slurp = get_tace($class, $id, $database);
  @db_slurp = grep {!/^(acedb\>|\/\/)$/} @db_slurp; # remove tace stuff from output
  map {$_ =~ s/^(\S+).*/$1/} @db_slurp; # replace each element of the array with the first word of each element
  map {$_ =~ s/\n//} @db_slurp; # remove newline

  foreach my $tag (@db_slurp) {
    $DATA->{$species}->{$version}->{$class}->{$id}->{$tag}++;
  }

  my $this_counts      = $DATA->{$species}->{$version}->{$class}->{$id};
  my $prev_counts      = (exists $DATA->{$species}->{$prev_version}) ?  $DATA->{$species}->{$prev_version}->{$class}->{$id} : {};
  my $prev_prev_counts = (exists $DATA->{$species}->{$prev_prev_version}) ?  $DATA->{$species}->{$prev_prev_version}->{$class}->{$id} : {};  

  ##################################################
  # Calculate difference between databases         #
  ##################################################

  my @tags = uniq (keys %$this_counts, keys %$prev_counts, keys %$prev_prev_counts); 

  foreach my $tag (sort @tags) {
    my $count2 = $this_counts->{$tag} || 0;
    my $count1 = (exists $prev_counts->{$tag}) ? $prev_counts->{$tag} : 0;
    my $count0 = (exists $prev_prev_counts->{$tag}) ? $prev_prev_counts->{$tag} : 0;

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
    
    $log->write_to(sprintf("%-22s %-22s %-22s %7d %7d %7d %10d %s\n", $class, $id, $tag, $count0, $count1, $count2, $diff, $err)) if ($err || $verbose);
  }

}

##################################################################
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
sub read_counts {
  my ($file) = @_;
  
  if (open (my $fh, $file)) {
    while (<$fh>) {
      next if not /\S/;
      chomp;
      my ($cc_species, $cc_version, $cc_class, $cc_id, $cc_tag, $cc_count) = split(/\s+/, $_);
      $DATA->{$cc_species}->{$cc_version}->{$cc_class}->{$cc_id}->{$cc_tag} = $cc_count;
    }
  }
}  



##################################################################
sub store_counts {
  my ($file) = @_;
  
  open(my $fh, ">$file") or $log->log_and_die("Could not open $file for writing\n");
  
  foreach my $sp (sort keys %$DATA) {
    foreach my $v (sort { $a <=> $b } keys %{$DATA->{$sp}}) {
      foreach my $cl (sort keys %{$DATA->{$sp}->{$v}}) {
        foreach my $id (sort keys %{$DATA->{$sp}->{$v}->{$cl}}) {
          foreach my $tg (sort keys %{$DATA->{$sp}->{$v}->{$cl}->{$id}}) {
            my $cnt = $DATA->{$sp}->{$v}->{$cl}->{$id}->{$tg};

            print $fh join("\t", $sp, $v, $cl, $id, $tg, $cnt), "\n";
          } 
        }
      }
    }                   
  }

  close($fh);
}




##################################################################
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
			       Protein => 'CE10000',
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
			       Transcript => 'B0205.9.1',
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
			       Pseudogene => 'CBG25555',
			       Sequence => ['cb25.fpc0002', 'cb25.fpc0011c', 'cb25.fpc0081', 'cb25.fpc0143a', 'chrI', 'chrII'],
			       Transcript => 'CBG00122.1',
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
			       Transcript => 'CBN00079.1',
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
			       Transcript => 'CJA00088.1',
			       },
		 'pristionchus' => {
			       CDS => 'PPA00099',
			       Feature_data => 'AA191935:polyA_site',
			       Gene => 'WBGene00089610', # misses Ortholog Ortholog_other Other_name
			       Gene_name => 'Ppa-abcf-1',
			       Sequence => ['PPA_ChrI', 'PPA_ChrII', 'PPA_ChrIV', 'PPA_ChrX'],
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
			       Protein => 'BM32546',
			       Pseudogene => 'Bm477',
			       Sequence => ['Bm_024', 'Bm_013', 'Bm_007', 'Bm_008', 'Bm_014', 'Bm_v4_Chr4_scaffold_001'],
			       Transcript => 'Bm1.1',
			     },
		 'ovolvulus' => {
			       Ace2SO => 'coding_transcript_ace2so',
			       CDS => 'OVOC9637',
			       Feature_data => 'OVOC_OO_000024:TRF',
			       Gene => 'WBGene00246446',
			       Gene_name => 'Ovo-eat-4',
			       Homol_data => 'OVOC_OO_000024:wublastx_brenneri',
			       Method => 'BLAT_EST_BEST',
			       Protein => 'OVP01471',
			       Sequence => ['OVOC_OO_000001', 'OVOC_OO_000008', 'OVOC_OO_000054', 'OVOC_OO_000132', 'OVOC_OO_000629', 'OVOC_OO_000690'],
			       Transcript => 'OVOC8637.1',
			     },
                  'sratti'  => { 
			       Ace2SO => 'coding_transcript_ace2so',
			       CDS => 'SRAE_0000000800',
			       Feature_data => 'SRAE_scaffold23:TRF',
			       Gene => 'XXX',
			       Gene_name => 'SRAE_0000000800',
			       Homol_data => 'SRAE_0000000800:wublastx_brenneri',
			       Method => 'BLAT_EST_BEST',
			       Protein => 'SRP01471',
			       Sequence => ['SRAE_chr2', 'SRAE_scaffold1'],
			       Transcript => 'SRAE_0000000800.1',
                             },

		);                  'tmuris'  => { 
			       Ace2SO => 'coding_transcript_ace2so',
			       CDS => 'TMUE_0000001900',
			       Feature_data => 'TMUE_scaffold351:TRF',
			       Gene => 'XXX',
			       Gene_name => 'TMUE_0000002363',
			       Homol_data => 'TMUE_0000001503:wublastx_brenneri',
			       Method => 'BLAT_EST_BEST',
			       Protein => 'TMP00002',
			       Sequence => ['TMUE_LG3', 'TMUE_scaffold273'],
			       Transcript => 'TMUE_0000001726.1',
                             },

  return $classes{$species};
}

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

