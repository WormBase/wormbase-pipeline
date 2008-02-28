#!/software/bin/perl -w
#
# find_anomalies.pl                           
# 
# by Gary Williams                        
#
# This looks for anomalous things such as protein homologies not
# matching a CDS and stores the results in the mysql database
# 'worm_anomaly'
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2008-02-28 09:22:15 $      

# Changes required by Ant: 2008-02-19
# 
# Improve signal to noise ration generally by changing thresholds and
# heavily weighting realy good anomalies.
# 
# Confirmed_intron - weight 10 (Done 2008-02-21)
# 
# Frameshifted_protein - improve.  (not done yet)
#  - Change from linear scoring to
# log10(Blast score) with min threshold of 50 giving an anomaly score
# of 1.x to 3 (Done 2008-02-20)
# 
# long_exon - remove (Done 2008-02-20)
# 
# merge_genes_by_proteins - want to use coordinate data so that false
# positives are filtered out. (not done yet)
# - Score as for merge_genes_by_proteins (Done 2008-02-20)
# 
# mismatched_est - output to a separate file for curators to look at (not done yet)
# 
# overlapping_exons - score increased to 5 (Done 2008-02-20)
# 
# repeat_overlaps_exons - ignore overlaps of circa <20 bases.  (Done 2008-02-22)
# - Report the types and frequency of overlaps. (not done yet)
# 
# short_exons - remove (Done 2008-02-20)
# 
# split_genes_by_protein - score as for merge_genes_by_proteins (Done 2008-02-20)
# 
# split_gene_by_protein_groups - score as for merge_genes_by_proteins (Done 2008-02-20)
# 
# unattached_tsl - remove (Done 2008-02-20)
# 
# unmatched_protein - score as for merge_genes_by_proteins (Done 2008-02-20)
# 
# unmatched_sage - make use of the score which is now in the GFF file (Done 2008-02-21)
# 
# unmatched_tsl - score increased to 5.   (Done 2008-02-20)
# - Add TSL inside genes but not near start of isoform. (not done yet)
# 
# unmatched_waba - get anomaly score 1 to 2 as a log10 transform of
# the waba score.  (Done 2008-02-20)
# 
# weak_intron_splice_site - Check splice site against confirmed introns.  (Done 2008-02-26)
# - Make use of splice score GFF files output by giface. (not done yet)
# 
# unmatched_genefinder - don't report if it overlaps with a repeat (Done 2008-02-22)
# 
# unmatched_twinscan - don't report if it overlaps with a repeat (Done 2008-02-22)
# 
# 
# Add:
# 
# incomplete_pfam_motif - speak to Rob Finn (pfam project leader)
# about the utility of these. (Mail sent: 2008-02-20)
# 
# Compare SignalP in known proteins to the number I predicted and send to Ant.
# Have as anomaly. (not done yet)

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
use Coords_converter;
use Modules::Overlap;
use Modules::PWM;
use DBI;
use POSIX qw(log10);

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($database, $species, $supplementary, $nodb);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "database:s" => \$database,	    # use the specified database instead of currentdb
	    "species:s"  => \$species,
	    "supplementary" => \$supplementary, # add GFF files to the SUPPLEMENTARY_DATA directory of the specified database
	    "nodb"       => \$nodb,
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

my $currentdb = $wormbase->database('current');		# set currentdb explicitly as the default
$database = $currentdb unless $database;


#################################
# Set up some useful paths      #
#################################

my $tace = $wormbase->tace;        # TACE PATH

##########################
# MAIN BODY OF SCRIPT
##########################


my $mysql;
my $db_key_id;
my $go_faster_by_ignoring_db_checks;
$species = $wormbase->species;
if (! $nodb) {
  # mysql database parameters
  my $sqldb = "worm_anomaly";
  if ($species ne 'elegans') {$sqldb = $sqldb . "_" . lc $species;}
  my $dbsn = "DBI:mysql:database=$sqldb;host=ia64b";
  my $dbuser = "wormadmin";
  my $dbpass = "worms";

  $mysql = DBI -> connect($dbsn, $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to database $sqldb, $DBI::errstr";

  # get the last used anomaly_id key value
  my $array_ref = $mysql->selectcol_arrayref("select max(anomaly_id) from anomaly;");
  $db_key_id = $array_ref->[0]; # 
  $go_faster_by_ignoring_db_checks = 0;
  if (defined $db_key_id) {
    print "db_key_id=$db_key_id\n";
  } else {
    # reset the database key value
    $db_key_id = 0; 
    print "db_key_id has been reset to 0\n";
    $go_faster_by_ignoring_db_checks = 1;	# don't need to check records in the database because there are none
  }
}

# output file of data to write to database, primarily for St. Louis to read in
my $datafile = $wormbase->wormpub . "/CURATION_DATA/anomalies.dat";
open (DAT, "> $datafile") || die "Can't open $datafile\n";

# and output the species line for the St. Louis datafile
print DAT "SPECIES\t$species\n\n";

my $coords = Coords_converter->invoke($database, 0, $wormbase);

# open an ACE connection to parse details for mapping to genome
print "Connecting to Ace\n";
my $ace = Ace->connect (-path => $database,
                       -program => $tace) || die "cannot connect to database at $database\n";

# get the Overlap object
my $ovlp = Overlap->new($database, $wormbase);

my $pwm = PWM->new;

my $chromosome_prefix=$wormbase->chromosome_prefix;

my %clonesize;
$wormbase->FetchData('clonesize', \%clonesize, "$database/COMMON_DATA/");

my %clonelab;
$wormbase->FetchData('clone2centre', \%clonelab, "$database/COMMON_DATA/");


# now delete things that have not been updated in this run that you
# would expect to have been updated like protein-homology-based
# anomalies that might have gone away.  This also deletes anomalies
# that we are no longer putting into the database and which can be
# removed.
&delete_anomalies("UNMATCHED_PROTEIN");
&delete_anomalies("SPLIT_GENES_BY_PROTEIN");
&delete_anomalies("SPLIT_GENE_BY_PROTEIN_GROUPS");
&delete_anomalies("SPLIT_GENES_BY_EST");
&delete_anomalies("MERGE_GENES_BY_EST");
&delete_anomalies("UNMATCHED_EST");
&delete_anomalies("UNATTACHED_EST");
&delete_anomalies("UNATTACHED_TSL");
&delete_anomalies("UNMATCHED_TSL");
&delete_anomalies("FRAMESHIFTED_PROTEIN");
&delete_anomalies("OVERLAPPING_EXONS");
&delete_anomalies("LONG_EXON");
&delete_anomalies("SHORT_EXON");
&delete_anomalies("UNMATCHED_WABA");
&delete_anomalies("UNMATCHED_SAGE");
&delete_anomalies("REPEAT_OVERLAPS_EXON");
&delete_anomalies("MERGE_GENES_BY_PROTEIN");
&delete_anomalies("WEAK_INTRON_SPLICE_SITE");
&delete_anomalies("UNMATCHED_TWINSCAN");
&delete_anomalies("UNMATCHED_GENEFINDER");
&delete_anomalies("CONFIRMED_INTRON");
&delete_anomalies("CONFIRMED_EST_INTRON");
&delete_anomalies("CONFIRMED_cDNA_INTRON");
&delete_anomalies("INTRONS_IN_UTR");
&delete_anomalies("SPLIT_GENE_BY_TWINSCAN");
&delete_anomalies("MERGE_GENE_BY_TWINSCAN");


my $ace_output = $wormbase->wormpub . "/CURATION_DATA/anomalies.ace";
if ($species ne 'elegans') { $ace_output = $wormbase->wormpub . "/CURATION_DATA/anomalies_" . lc $species . ".ace";}
open (OUT, "> $ace_output") || die "Can't open $ace_output to write the Method\n";
print OUT "\n\n";
print OUT "Method : \"curation_anomaly\"\n";
print OUT "Colour   RED\n";
print OUT "Strand_sensitive Show_up_strand\n";
print OUT "Right_priority   1.45\n";
print OUT "Width   2\n";
print OUT "Score_by_width\n";
print OUT "Score_bounds 0.01 1.0\n";
print OUT "Remark \"This method is used by acedb to display curation anomaly regions\"\n";
print OUT "\n\n";
close(OUT);

# loop through the chromosomes
my @chromosomes = $wormbase->get_chromosome_names(-mito => 0, -prefix => 0);

foreach my $chromosome (@chromosomes) {

  # if we want the anomalies GFF file
  if ($supplementary) {
    my $gff_file = "$database/CHROMOSOMES/SUPPLEMENTARY_GFF/${chromosome_prefix}${chromosome}_curation_anomalies.gff";
    open (OUTPUT_GFF, ">$gff_file") || die "Can't open $gff_file";
  }

  $log->write_to("Processing chromosome $chromosome\n");

  print "reading GFF data\n";
  my @est_hsp = $ovlp->get_EST_BEST($chromosome);
  #my @est_paired_span = $ovlp->get_paired_span(@est_hsp); # change the ESTs from HSPs to start-to-end span of paired reads
  my @cds_exons  = $ovlp->get_curated_CDS_exons($chromosome);
  my @pseudogenes = $ovlp->get_pseudogene($chromosome);
  my @rRNA = $ovlp->get_rRNA_transcripts($chromosome);
  #my @coding_transcript_exons = $ovlp->get_Coding_transcript_exons($chromosome);
  my @coding_transcripts = $ovlp->get_Coding_transcripts($chromosome);
  my @transposons = $ovlp->get_transposons($chromosome);
  my @transposon_exons = $ovlp->get_transposon_exons($chromosome);
  my @noncoding_transcript_exons = $ovlp->get_Non_coding_transcript_exons($chromosome);
  my @TSL_SL1 = $ovlp->get_TSL_SL1($chromosome);
  my @TSL_SL2 = $ovlp->get_TSL_SL2($chromosome);
  my @homologies = $ovlp->get_blastx_homologies($chromosome);
  my @twinscan_exons = $ovlp->get_twinscan_exons($chromosome);
  my @twinscan_transcripts = $ovlp->get_twinscan_transcripts($chromosome);
  my @genefinder = $ovlp->get_genefinder_exons($chromosome);
  my @waba_coding = $ovlp->get_waba_coding($chromosome);
  my @repeatmasked = $ovlp->get_repeatmasked($chromosome);
  my @UTRs_5 = $ovlp->get_5_UTRs($chromosome);
  my @UTRs_3 = $ovlp->get_3_UTRs($chromosome);
  my @CDS = $ovlp->get_curated_CDS($chromosome);
  my @SAGE_tags = $ovlp->get_SAGE_tags($chromosome);
  my @miRNA = $ovlp->get_miRNA($chromosome);
  my @ncRNA = $ovlp->get_ncRNA($chromosome);
  my @scRNA = $ovlp->get_scRNA($chromosome);
  my @snRNA = $ovlp->get_snRNA($chromosome);
  my @snoRNA = $ovlp->get_snoRNA($chromosome);
  my @stRNA = $ovlp->get_stRNA($chromosome);
  my @tRNA = $ovlp->get_tRNA($chromosome);
  my @tRNAscan_SE_1_23 = $ovlp->get_tRNAscan_SE_1_23RNA($chromosome);
  my @repeatmasked_complex = $ovlp->get_repeatmasked_complex(@repeatmasked);
  my @CDS_introns = $ovlp->get_CDS_introns($chromosome);
  my @check_introns_EST  = $ovlp->get_check_introns_EST($chromosome);
  my @check_introns_cDNA = $ovlp->get_check_introns_cDNA($chromosome);


######################################################################

  print "finding anomalies\n";

  print "finding protein homologies not overlapping CDS exons\n";
  my $matched_protein_aref = &get_protein_differences(\@cds_exons, \@pseudogenes, \@homologies, \@transposons, \@transposon_exons, \@noncoding_transcript_exons, \@rRNA, $chromosome);

  print "finding frameshifts\n";
  &get_frameshifts($matched_protein_aref, $chromosome);

  print "finding genes to be split/merged based on protein homology\n";
  &get_protein_split_merged($matched_protein_aref, $chromosome);

  print "finding genes to be split based on protein homology groups\n";
  &get_protein_split($matched_protein_aref, $chromosome);

  # this finds TEC-RED TSL sites more than 100 bases upstream that are not mentioned in the remarks or evidence
  print "finding isolated TSL sites\n";
  &get_isolated_TSL(\@TSL_SL1, \@TSL_SL2, \@coding_transcripts, \@pseudogenes, \@transposons, \@transposon_exons, \@noncoding_transcript_exons, \@rRNA, $chromosome);

  # get SAGE tags that don't match a gene with score based on frequency
  print "finding non-overlapping SAGE_tags\n";
  &get_unmatched_SAGE(\@coding_transcripts, \@pseudogenes, \@SAGE_tags, \@transposons, \@transposon_exons, \@noncoding_transcript_exons, \@rRNA, \@miRNA, \@ncRNA, \@scRNA, \@snRNA, \@snoRNA, \@stRNA, \@tRNA, \@tRNAscan_SE_1_23, $chromosome);

  print "finding twinscan exons not overlapping CDS exons\n";
  &get_unmatched_twinscan_exons(\@twinscan_exons, \@cds_exons, \@pseudogenes, \@transposons, \@transposon_exons, \@noncoding_transcript_exons, \@rRNA, \@repeatmasked_complex, $chromosome);

  print "finding genefinder exons not overlapping CDS exons\n";
  &get_unmatched_genefinder_exons(\@genefinder, \@cds_exons, \@pseudogenes, \@transposons, \@transposon_exons, \@noncoding_transcript_exons, \@rRNA, \@repeatmasked_complex, $chromosome);

  print "finding WABA coding regions not overlapping CDS exons\n";
  &get_unmatched_waba_coding(\@waba_coding, \@cds_exons, \@pseudogenes, \@transposons, \@transposon_exons, \@noncoding_transcript_exons, \@rRNA, \@repeatmasked_complex, $chromosome);

  print "finding CDS exons overlapping repeatmasker regions\n";
  &get_matched_repeatmasker(\@cds_exons, \@repeatmasked_complex, $chromosome);

  print "finding CDS exons overlapping other genes\n";
  &get_matched_exons(\@cds_exons, \@pseudogenes, \@transposons, \@transposon_exons, \@noncoding_transcript_exons, \@rRNA, $chromosome);

  print "finding EST/genome mismatches\n";
  &get_est_mismatches(\@est_hsp, $chromosome);


  print "finding weak CDS exon splice sites\n";
  &get_weak_exon_splice_sites(\@CDS_introns, $chromosome);

  print "finding multiple UTR introns\n";
  &get_multiple_utr_introns(\@UTRs_5, $chromosome);
  &get_multiple_utr_introns(\@UTRs_3, $chromosome);

  print "read confirmed checked introns\n";
  &get_checked_confirmed_introns(\@check_introns_EST, \@check_introns_cDNA, $chromosome);

  print "finding genes to be split/merged based on twinscan\n";
  &get_twinscan_split_merged(\@twinscan_transcripts, \@CDS, $chromosome);


#################################################
# these don't work very well - don't use

##  this looks at the ESTs  finds those not attached to a transcript
##  print "finding ESTs sites not attached to a transcript\n";
##  &get_unattached_EST(\@est_hsp, $chromosome);

##  # ESTs not matching exons
##  # and those with a match to a coding exon
##  print "finding EST homologies not overlapping exons\n";
##  my $matched_EST_aref = &get_EST_differences(\@exons, \@pseudogenes, \@est_hsp, \@transposons, \@transposon_exons, \@noncoding_transcript_exons, \@rRNA, $chromosome);

##  print "finding split/merged genes based on EST homology\n";
##  &get_EST_split_merged($matched_EST_aref, $chromosome, %Show_in_reverse_orientation);

#################################################
# we don't value these anomalies very much - don't run them for now

##  print "finding short CDS exons\n";
##  &get_short_exons(\@coding_transcripts, $chromosome);

#  # this looks at the EST TSL sites and finds those not attached to genes
#  print "finding TSL sites not attached to genes\n";
#  &get_unattached_TSL(\@TSL_SL1, \@TSL_SL2, $chromosome);


#################################################





  # close the ouput GFF file
  if ($supplementary) {
    close (OUTPUT_GFF);
  }
}

# close the output datafile for St. Louis
if ($datafile) {
  close(DAT);			
}

# disconnect from the mysql database
if (! $nodb) {
  $mysql->disconnect || die "error disconnecting from database", $DBI::errstr;
}

# close the ACE connection
$ace->close;


##################
# Check the files
##################

if ($database eq $wormbase->{'autoace'}) {
  foreach my $chromosome (@chromosomes) {
    my $gff_file = $wormbase->{'chromosomes'} . "/SUPPLEMENTARY_GFF/${\$wormbase->chromosome_prefix}${chromosome}_curation_anomalies.gff";
    $wormbase->check_file($gff_file, $log,
			  minsize => 700000,
			  maxsize => 4000000,
			  lines => ['^##',
				    "^${\$wormbase->chromosome_prefix}${chromosome}\\s+curation_anomaly\\s+\\S+\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+\\s+Evidence\\s+\\S+"],
			  );
  }
}

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################
##########################################
# get the homologies with no matching exons or pseudogenes or transposons
# and those which do match exons or transposons (but not pseudogenes)

#  my @matched_homologies = get_differences(\@transcripts, \@pseudogenes, \@protein_coverage);

sub get_protein_differences {
  my ($exons_aref, $pseudogenes_aref, $homologies_aref, $transposons_aref, $transposon_exons_aref, $noncoding_transcript_exons_aref, $rRNA_aref, $chromosome) = @_;

  my @homologies = @{$homologies_aref};

  my @not_matched = ();		# the resulting list of hashes of homologies with no matching exons/transposons/pseudogenes
  my @matched = ();		# the resulting list of hashes of homologies which match a coding exon

  my $SMALL_OVERLAP = -20;	# amount of small overlap to ignore
      
  # we allow the protein to be in either sense compared to the coding exons
  my $exons_match = $ovlp->compare($exons_aref, near => $SMALL_OVERLAP, same_sense => 0);
  my $pseud_match = $ovlp->compare($pseudogenes_aref);
  my $trans_match = $ovlp->compare($transposons_aref);
  my $trane_match = $ovlp->compare($transposon_exons_aref);
  my $nonco_match = $ovlp->compare($noncoding_transcript_exons_aref);
  my $rrna_match  = $ovlp->compare($rRNA_aref);

  foreach my $homology (@homologies) { # $protein_id, $chrom_start, $chrom_end, $chrom_strand

    # don't want the low-scoring trembl proteins - they are not informative
    if ($homology->[0] =~ /TR:/ &&  $homology->[6] < 200) {next;}

    my $got_a_match = 0;	        # not yet seen a match to anything
    my $got_a_match_to_coding_exon = 0;	# not yet seen a match to a coding exon
    my $matching_exon = "";	        # the name of the exon that matches;

    if (my @matching_exons = $exons_match->match($homology)) {               #&match($homology, $exons_aref, \%exons_match)) 
      $got_a_match = 1;
      $got_a_match_to_coding_exon = 1;
      $matching_exon = $matching_exons[0][0];
    }

    if ( $pseud_match->match($homology)) {              #&match($homology, $pseudogenes_aref, \%pseud_match)) {
      $got_a_match = 1;
    }

    if (my @matching_exons = $trans_match->match($homology)) {               #&match($homology, $transposons_aref, \%trans_match)) {
      $got_a_match = 1;
      $got_a_match_to_coding_exon = 1;
      $matching_exon = $matching_exons[0][0];
    }

    if ($trane_match->match($homology)) {               #&match($homology, $transposon_exons_aref, \%trane_match)) {
      $got_a_match = 1;
    }

    if ($nonco_match->match($homology)) {               #&match($homology, $noncoding_transcript_exons_aref, \%nonco_match)) {
      $got_a_match = 1;
    }

    if ($rrna_match->match($homology)) {                #&match($homology, $rRNA_aref, \%rrna_match)) {
      $got_a_match = 1;
    }

    # output unmatched homologies to the database
    if (! $got_a_match) {
      my $protein_id = $homology->[0];
      my $chrom_start = $homology->[1];
      my $chrom_end = $homology->[2];
      my $chrom_strand = $homology->[3];
      my $protein_score = $homology->[6];

      # reject any hits where either of the proteins have a Blast score < 50
      if ($protein_score < 50) {next;}

      # make the anomaly score based on the protein alignment score normalised between 1 and 3
      # using the log10 of the blast score
      # the BLAST scores seem to be between 0 and 1000
      # use the average of the two alignment's scores
      my $anomaly_score = POSIX::log10($protein_score);
      if ($anomaly_score > 3) {$anomaly_score = 3;}
      if ($anomaly_score < 0) {$anomaly_score = 0;}

      #print "NOT got a match ANOMALY: $protein_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score\n";
      &output_to_database("UNMATCHED_PROTEIN", $chromosome, $protein_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, '');
    }

    # add to the list of homologies to a coding regions
    if ($got_a_match_to_coding_exon) {
      push @{$homology}, $matching_exon; # make a note of the exon we matched
      push @matched, $homology;

    }
  }

  return (\@matched);
}

##########################################
# get the EST homologies with no matching exons or pseudogenes or transposons
# and those which do match exons or transposons (but not pseudogenes)
# my $matched_EST_aref = &get_EST_differences(\@exons, \@pseudogenes, \@est, \@transposons, \@transposon_exons, \@noncoding_transcript_exons, $chromosome);

sub get_EST_differences {
  my ($exons_aref, $pseudogenes_aref, $est_aref, $transposons_aref, $transposon_exons_aref, $noncoding_transcript_exons_aref, $rRNA_aref, $chromosome) = @_;

  my @est = @{$est_aref};

  my @not_matched = ();		# the resulting list of hashes of est with no matching exons/transposons/pseudogenes
  my @matched = ();		# the resulting list of hashes of est which match a coding exon


  my $exons_match = $ovlp->compare($exons_aref, same_sense => 1);
  my $pseud_match = $ovlp->compare($pseudogenes_aref);
  my $trans_match = $ovlp->compare($transposons_aref);
  my $trane_match = $ovlp->compare($transposon_exons_aref);
  my $nonco_match = $ovlp->compare($noncoding_transcript_exons_aref);
  my $rrna_match  = $ovlp->compare($rRNA_aref);


  foreach my $homology (@est) { # $est_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_match = 0;	        # not yet seen a match to anything
    my $got_a_match_to_coding_exon = 0;	# not yet seen a match to a coding exon
    my $matching_exon = "";	        # the name of the exon that matches;

    if (my @matching_exons = $exons_match->match($homology)) {               #&match($homology, $exons_aref, \%exons_match)) 
      $got_a_match = 1;
      $got_a_match_to_coding_exon = 1;
      $matching_exon = $matching_exons[0][0];
    }

    if ( $pseud_match->match($homology)) {              #&match($homology, $pseudogenes_aref, \%pseud_match)) {
      $got_a_match = 1;
    }

    if (my @matching_exons = $trans_match->match($homology)) {               #&match($homology, $transposons_aref, \%trans_match)) {
      $got_a_match = 1;
      $got_a_match_to_coding_exon = 1;
      $matching_exon = $matching_exons[0][0];
    }

    if ($trane_match->match($homology)) {               #&match($homology, $transposon_exons_aref, \%trane_match)) {
      $got_a_match = 1;
    }

    if ($nonco_match->match($homology)) {               #&match($homology, $noncoding_transcript_exons_aref, \%nonco_match)) {
      $got_a_match = 1;
    }

    if ($rrna_match->match($homology)) {                #&match($homology, $rRNA_aref, \%rrna_match)) {
      $got_a_match = 1;
    }

    # output unmatched est to the database
    if (! $got_a_match) {
      my $est_id = $homology->[0];
      my $chrom_start = $homology->[1];
      my $chrom_end = $homology->[2];
      my $chrom_strand = $homology->[3];
      my $est_score = $homology->[6];
      # make the anomaly score based on the protein alignment score normalised between 0 and 1
      # BLAT scores are between 0 and 100
      my $anomaly_score = $est_score/100;
      if ($anomaly_score > 1) {$anomaly_score = 1;}
      if ($anomaly_score < 0) {$anomaly_score = 0;}
      #print "NOT got an EST match ANOMALY: $est_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score\n";
      &output_to_database("UNMATCHED_EST", $chromosome, $est_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, '');
    }

    # add to the list of homologies to a coding regions
    if ($got_a_match_to_coding_exon) {
      push @{$homology}, $matching_exon; # make a note of the exon we matched
      push @matched, $homology;

    }
  }

  return (\@matched);
}

##########################################

# get the frameshifts evident in adjacent homologies of a group of HSPs for proteins

#  my @frameshifts = get_frameshifts(\@matched);

sub get_frameshifts {
  my ($matched_aref, $chromosome) = @_;

  my @matched = @{$matched_aref};

  my $prev_protein_id = "";	# use to collect alignments for a protein's homology when looking for frameshifts
  my $prev_chrom_strand = "";	# the previous HSP's sense
  my $prev_chrom_end = -1;
  my $prev_hit_start = -1;
  my $prev_hit_end = -1;
  my $prev_score = -1;

  my @list_of_aligments = ();	# previous protein ID

  # sort the homologies grouped by protein ID and then chromosomal position
  my @homologies = sort {$a->[0] cmp $b->[0] or $a->[1] <=> $b->[1]} @matched;

  foreach my $homology (@homologies) { # $protein_database, $protein_id, $chrom_start, $chrom_end, $chrom_strand, $hit_start, $hit_end, $clone, $protein_score

    my $protein_id = $homology->[0];
    my $chrom_start = $homology->[1];
    my $chrom_end = $homology->[2];
    my $chrom_strand = $homology->[3];
    my $hit_start = $homology->[4];
    my $hit_end = $homology->[5];
    my $protein_score = $homology->[6];

    #print "Matched: $protein_id, $chrom_start, $chrom_end, $chrom_strand, $hit_start, $hit_end, $protein_score\n" if ($verbose);
    # look for frameshifts against the other alignment blocks of this homology
    if ($protein_id eq $prev_protein_id && $chrom_strand eq $prev_chrom_strand) {
      #print "Looking at next group\n" if ($verbose);

      #print $homology->[0] ."\t". $homology->[1] ."\t". $homology->[2] ."\t". $homology->[3] ."\t". $homology->[4] ."\t". $homology->[5]."\n";

      # look for lack of intron
      my $chrom_diff = $chrom_start - $prev_chrom_end;
      my $prot_diff = $hit_start - $prev_hit_end;
	  
      # output any frameshifts found to the database
      # want the frame to have changed, so look at the chrom_start to prev_chrom_end difference mod 3
      if (( (abs($chrom_diff+1)) % 3 != 0) && $chrom_diff > -30 && $chrom_diff < 15 && $prot_diff > -15 && $prot_diff < 15) {

	# get the region to display
	my $anomaly_start = $prev_chrom_end;
	my $anomaly_end = $chrom_start;
	# get the start and end in the right order
	if ($anomaly_end < $anomaly_start) {
	  my $tmp = $anomaly_end;
	  $anomaly_end = $anomaly_start;
	  $anomaly_start = $tmp; 
	}
	# expand the region to display a little
	$anomaly_start -= 10;
	$anomaly_end += 10;

	# reject any hits where either of the proteins have a Blast score < 50
	if ($protein_score < 50 || $prev_score < 50) {next;}

	# make the anomaly score based on the protein alignment score normalised between 1 and 3
	# using the log10 of the blast score
	# the BLAST scores seem to be between 0 and 1000
	# use the average of the two alignment's scores
	my $average_blast_score = ($protein_score+$prev_score)/2;
	my $anomaly_score = POSIX::log10($average_blast_score);
	if ($anomaly_score > 3) {$anomaly_score = 3;}
	if ($anomaly_score < 0) {$anomaly_score = 0;}

	#print "FRAMESHIFTED_PROTEIN ANOMALY: $protein_id, $anomaly_start, $anomaly_end, $chrom_strand, $anomaly_score\n";
	&output_to_database("FRAMESHIFTED_PROTEIN", $chromosome, $protein_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, '');
      }

    }      

    # update the previous HSP details
    $prev_chrom_end = $chrom_end;
    $prev_hit_start = $hit_start;
    $prev_hit_end = $hit_end;
    $prev_score = $protein_score;
    $prev_protein_id = $protein_id;
    $prev_chrom_strand = $chrom_strand;

  }

}



##########################################

#  my @split_merged = &get_split_merged($matched_aref);

# find groups of homology that indicate that the genes they match should be split or merged
# look for non-overlapping HSPs that jump back down the position in the protein that is aligned as you move along the chromosome -> split
# look for matches to exons with HSPs that do not jump back down the position in the protein that is aligned as you move along the chromosome -> merge

sub get_protein_split_merged {
  my ($matched_aref, $chromosome) = @_;

  my @matched = @{$matched_aref};

  my $prev_protein_id = "";	# use to collect alignments for a protein's homology when looking for frameshifts
  my $prev_exon = "";
  my $prev_hit_start = -1;
  my $prev_hit_end = -1;
  my $prev_chrom_start = -1;
  my $prev_chrom_end = -1;
  my $prev_score = 0;
  my $prev_chrom_strand = "";

  my $got_a_new_exon;
  my $got_a_big_decrease_in_HSP_start;
  my $got_a_continuation_of_the_HSPs;


  # sort the homologies grouped by protein ID and then chromosomal position
  my @homologies = sort {$a->[0] cmp $b->[0] or $a->[1] <=> $b->[1]} @matched;

  foreach my $homology (@homologies) { # $protein_id, $chrom_start, $chrom_end, $chrom_strand, $hit_start, $hit_end, $protein_score, $matching_exon

    my $protein_id = $homology->[0];
    my $chrom_start = $homology->[1];
    my $chrom_end = $homology->[2];
    my $chrom_strand = $homology->[3];
    my $hit_start = $homology->[4];
    my $hit_end = $homology->[5];
    my $protein_score = $homology->[6];
    # $homology->[7] can hold the the other_data field from the GFF line, but is usually undef
    my $matching_exon = $homology->[8];

    # we don't use the briggsae or remanei proteins in this analysis
    # because their predictions are heavily based on twinscan
    # predictions and so we would be simply confirming possibly
    # erroneous twinscan data.
    if ($protein_id =~ /^BP:/ || $protein_id =~ /^RP:/) {next;}

    #  check if these are isoforms of the same gene and if so then treat them as the same gene
    if ($matching_exon =~ /(\S+\.\d+)[a-z]/) {
      $matching_exon = $1;
    }
    # and change the strange transcript names like H25P06.1.1, H25P06.1.2 to H25P06.1
    if ($matching_exon =~ /(\S+\.\d+)\.\d+/) {
      $matching_exon = $1;
    }

    #print "Matched: $protein_id, $chrom_start, $chrom_end, $chrom_strand, $hit_start, $hit_end, $protein_score, $matching_exon\n";

    # look for homologies which cover two genes or genes which cover a repeated homology
    if ($protein_id eq $prev_protein_id && $chrom_strand eq $prev_chrom_strand) {

      # see if exon has changed
      if ($prev_exon ne "" && $prev_exon ne $matching_exon) {
	$got_a_new_exon = 1;
	#print "Merge? prev $prev_exon this $matching_exon\n";
      } else {
	$got_a_new_exon = 0;
      }

      # see if the protein HSP order is unchanged - a possible case for merging if continued over into the next gene
      # same sense, non-overlapping on the chromosome and jump in HSP start
      #
      # We want to avoid picking up the case where there are two duplicated single-exon genes
      # and the homologies have a hit_start that differs by only a few bases - this should be ignored.
      # So we look for cases where the start continues as greater than
      # the previous mid point for evidence of a possible merge.
      #
      my $prev_midpoint = $prev_hit_start + ($prev_hit_end - $prev_hit_start)/2;
      my $midpoint = $hit_start + ($hit_end - $hit_start)/2;
      if (($chrom_strand eq '+' && $prev_chrom_end < $chrom_start && $hit_start > $prev_midpoint) ||
	  ($chrom_strand eq '-' && $prev_chrom_end < $chrom_start && $hit_start < $prev_midpoint && $prev_hit_start != -1)) {
	$got_a_continuation_of_the_HSPs = 1;
      } else {
	$got_a_continuation_of_the_HSPs = 0;
      }
      
      # see if we have a big decrease in the start of the HSP - a possible case for splitting
      if (($chrom_strand eq '+' && $prev_chrom_end < $chrom_start && $prev_hit_start > $hit_start + 100) ||
	  ($chrom_strand eq '-' && $prev_chrom_end < $chrom_start && $prev_hit_start < $hit_start - 100 && $prev_hit_start != -1)) {
	$got_a_big_decrease_in_HSP_start = 1;
	#print "Split? $protein_id chrom_start $chrom_end prev_chrom_start $prev_chrom_start $chrom_strand prev_hit_start $prev_hit_start hit_start $hit_start\n";
      } else {
	$got_a_big_decrease_in_HSP_start = 0;
      }


      # we have a merge if the pattern of HSPs shows a continuation of
      # the protein matching in order and this is a new gene and the
      # distance between the alignments is less than 5 kb
      if ($got_a_new_exon && $got_a_continuation_of_the_HSPs && $chrom_start - $prev_chrom_end < 5000) {
	# reject any hits where either of the proteins have a Blast score < 50
	if ($protein_score < 50 || $prev_score < 50) {next;}

	# make the anomaly score based on the protein alignment score normalised between 1 and 3
	# using the log10 of the blast score
	# the BLAST scores seem to be between 0 and 1000
	# use the average of the two alignment's scores
	my $average_blast_score = ($protein_score+$prev_score)/2;
	my $anomaly_score = POSIX::log10($average_blast_score);
	if ($anomaly_score > 3) {$anomaly_score = 3;}
	if ($anomaly_score < 0) {$anomaly_score = 0;}

	if ($anomaly_score > 0.1) {
	  &output_to_database("MERGE_GENES_BY_PROTEIN", $chromosome, $protein_id, $prev_chrom_start, $chrom_end, $chrom_strand, $anomaly_score, "merge $prev_exon and $matching_exon");
	}

      }

      if ($got_a_big_decrease_in_HSP_start && ! $got_a_new_exon) {
	# reject any hits where either of the proteins have a Blast score < 50
	if ($protein_score < 50 || $prev_score < 50) {next;}

	# make the anomaly score based on the protein alignment score normalised between 1 and 3
	# using the log10 of the blast score
	# the BLAST scores seem to be between 0 and 1000
	# use the average of the two alignment's scores
	my $average_blast_score = ($protein_score+$prev_score)/2;
	my $anomaly_score = POSIX::log10($average_blast_score);
	if ($anomaly_score > 3) {$anomaly_score = 3;}
	if ($anomaly_score < 0) {$anomaly_score = 0;}

	&output_to_database("SPLIT_GENES_BY_PROTEIN", $chromosome, $protein_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, "split $matching_exon");
      }
    }

    $prev_exon = $matching_exon;
    $prev_hit_start = $hit_start;
    $prev_hit_end = $hit_end;
    $prev_score = $protein_score;
    $prev_chrom_start = $chrom_start;
    $prev_chrom_end = $chrom_end;
    $prev_protein_id = $protein_id;
    $prev_chrom_strand = $chrom_strand;
  }

}


##########################################

#  my &get_protein_split($matched_aref, $chromosome);

# find groups of homology that indicate that the genes they match should be split
# look for genes that have two or more non-overlapping groups of homology from different sets of genes

sub get_protein_split {
  my ($matched_aref, $chromosome) = @_;

  my @matched = @{$matched_aref};

  my $prev_transcript = "";

  my @brigpep = ();
  my @rempep = ();
  my @wormpep = ();

  # sort the homologies grouped by matching transcript and then protein_id
  my @homologies = sort {$a->[8] cmp $b->[8] or $a->[0] cmp $b->[0]} @matched;

  foreach my $homology (@homologies) { # $protein_id, $chrom_start, $chrom_end, $chrom_strand, $hit_start, $hit_end, $protein_score, $matching_transcript

    my $protein_id = $homology->[0];
    my $chrom_start = $homology->[1];
    my $chrom_end = $homology->[2];
    my $chrom_strand = $homology->[3];
    my $hit_start = $homology->[4];
    my $hit_end = $homology->[5];
    my $protein_score = $homology->[6];
    # $homology->[7] can hold the the other_data field from the GFF line, but is usually undef
    my $matching_transcript = $homology->[8];

    # don't look at protein homologies with a score below 30
    if ($protein_score <= 30) {next;}

    #  check if these are isoforms of the same gene and if so then treat them as the same gene
    if ($matching_transcript =~ /(\S+\.\d+)[a-z]/) {
      $matching_transcript = $1;
    }
    # and change the strange transcript names like H25P06.1.1, H25P06.1.2 to H25P06.1
    if ($matching_transcript =~ /(\S+\.\d+)\.\d+/) {
      $matching_transcript = $1;
    }

    # do we have a new matching transcript? If so, see if the transcript we were doing should be split
    if ($matching_transcript ne $prev_transcript && $prev_transcript ne "") {

      #print "Got all homologies for $prev_transcript\n";

      # check if there are two or more different boxes with more than one protein in them and score over the box cutoff
      # if so, then output a database record for each contributing protein in the boxes
      my @box_results;
      if (@box_results = split_check(@brigpep)) {
	foreach my $pep (@box_results) {
	  my $pep_protein_id = $pep->[0];
	  my $pep_chrom_start = $pep->[1];
	  my $pep_chrom_end = $pep->[2];
	  my $pep_chrom_strand = $pep->[3];
	  my $pep_protein_score = $pep->[6];
	  my $pep_matching_transcript = $pep->[8];

	  # reject any hits where either of the proteins have a Blast score < 50
	  if ($pep_protein_score < 50) {next;}
	  # make the anomaly score based on the protein alignment score normalised between 1 and 3
	  # using the log10 of the blast score
	  # the BLAST scores seem to be between 0 and 1000
	  # use the average of the two alignment's scores
	  my $anomaly_score = POSIX::log10($pep_protein_score);
	  if ($anomaly_score > 3) {$anomaly_score = 3;}
	  if ($anomaly_score < 0) {$anomaly_score = 0;}

	  output_to_database("SPLIT_GENE_BY_PROTEIN_GROUPS", $chromosome, $pep_protein_id, $pep_chrom_start, $pep_chrom_end, $pep_chrom_strand, $anomaly_score, "split $pep_matching_transcript");

	}
      }

      if (@box_results = split_check(@rempep)) {
	foreach my $pep (@box_results) {
	  my $pep_protein_id = $pep->[0];
	  my $pep_chrom_start = $pep->[1];
	  my $pep_chrom_end = $pep->[2];
	  my $pep_chrom_strand = $pep->[3];
	  my $pep_protein_score = $pep->[6];
	  my $pep_matching_transcript = $pep->[8];

	  # reject any hits where either of the proteins have a Blast score < 50
	  if ($pep_protein_score < 50) {next;}
	  # make the anomaly score based on the protein alignment score normalised between 1 and 3
	  # using the log10 of the blast score
	  # the BLAST scores seem to be between 0 and 1000
	  # use the average of the two alignment's scores
	  my $anomaly_score = POSIX::log10($pep_protein_score);
	  if ($anomaly_score > 3) {$anomaly_score = 3;}
	  if ($anomaly_score < 0) {$anomaly_score = 0;}

	  output_to_database("SPLIT_GENE_BY_PROTEIN_GROUPS", $chromosome, $pep_protein_id, $pep_chrom_start, $pep_chrom_end, $pep_chrom_strand, $anomaly_score, "split $pep_matching_transcript");
	}
      }

      if (@box_results = split_check(@wormpep)) {
	foreach my $pep (@box_results) {
	  my $pep_protein_id = $pep->[0];
	  my $pep_chrom_start = $pep->[1];
	  my $pep_chrom_end = $pep->[2];
	  my $pep_chrom_strand = $pep->[3];
	  my $pep_protein_score = $pep->[6];
	  my $pep_matching_transcript = $pep->[8];

	  # reject any hits where either of the proteins have a Blast score < 50
	  if ($pep_protein_score < 50) {next;}
	  # make the anomaly score based on the protein alignment score normalised between 1 and 3
	  # using the log10 of the blast score
	  # the BLAST scores seem to be between 0 and 1000
	  # use the average of the two alignment's scores
	  my $anomaly_score = POSIX::log10($pep_protein_score);
	  if ($anomaly_score > 3) {$anomaly_score = 3;}
	  if ($anomaly_score < 0) {$anomaly_score = 0;}

	  output_to_database("SPLIT_GENE_BY_PROTEIN_GROUPS", $chromosome, $pep_protein_id, $pep_chrom_start, $pep_chrom_end, $pep_chrom_strand, $anomaly_score, "split $pep_matching_transcript");
	}
      }



      # initialise the new lists to be put into boxes
      @brigpep = ();
      @rempep = ();
      @wormpep = ();
    }

    # only want to store details from briggsae, remanei and elegans
    if ($protein_id =~ /^BP:/) {
	
      # add this protein to the list of homologies to be put into boxes
      push @brigpep, $homology;

    } elsif ($protein_id =~ /^RP:/) {

      # add this protein to the list of homologies to be put into boxes
      push @rempep, $homology;

    } elsif ($protein_id =~ /^WP:/) {
      
      # add this protein to the list of homologies to be put into boxes
      push @wormpep, $homology;
    }
    
    $prev_transcript = $matching_transcript;
  }

}

##########################################
#    @box_results = split_check(@wormpep);
# put the homologies into non-overlapping boxes
# find non-overlapping ranges of proteins that do not share protein IDs
# the boxes of protein homology must have at least two different proteins in them and must have a total Blast score greater than 100
# return the list of proteins that contribute to the two or more boxes

sub split_check {
  my (@peplist) = @_;

  my @box_result = ();


  my @boxes = ();
 
  my $last_sense = '';
 
  foreach my $homology (@peplist) {
    my $protein_id = $homology->[0];
    my $chrom_start = $homology->[1];
    my $chrom_end = $homology->[2];
    my $sense = $homology->[3];
    my $hit_start = $homology->[4];
    my $hit_end = $homology->[5];
    my $protein_score = $homology->[6];
    my $matching_transcript = $homology->[8];

    #print "$protein_id $chrom_start $chrom_end\n";

    $last_sense = $sense;
 
    my $got_a_match = 0;

    #print "Have @boxes boxes in the list\n";

    foreach my $box (@boxes) {
      # find an existing box to merge to
      # merge if it contains a previous block of this protein's homology or it has an overlap to this block of homology
      #print "$protein_id Box   " . $box->{'chrom_start'} . " " . $box->{'chrom_end'} . " IDs:  " , @{$box->{'ID'}} , "\n";
      if (((grep /$protein_id/, @{$box->{'ID'}}) ||
	  $box->{'chrom_start'} <= $chrom_end && $box->{'chrom_end'} >= $chrom_start) &&
          ! exists $box->{'deleted'}) {
                                                                                                                                              
        #print "*** Box $box->{'chrom_start'}..$box->{'chrom_end'}, matches $protein_id, $chrom_start, $chrom_end, $protein_score\n";
                                                                                                                                              
	# add this homology to the box
        $box->{'count'}++;        # add to the count of matches in this box
        $box->{'total_score'} += $protein_score;        # add to the sum of the alignment scores
        push @{$box->{'ID'}}, ($protein_id);        # add to the protein IDs
                                                                                                                                              
        # update start/end
        if ($box->{'chrom_start'} > $chrom_start) {$box->{'chrom_start'} = $chrom_start;}
        if ($box->{'chrom_end'} < $chrom_end) {$box->{'chrom_end'} = $chrom_end;}
                                                                                                                                              
        $got_a_match = 1;
                                                                                                                                              
        # check to see if we now need to merge two overlapping boxes
        my $past_first_box = 0;
        foreach my $other_box (@boxes) {
         if (! $past_first_box) {
            # see if this is our current $box
            if ($other_box->{'chrom_start'} == $box->{'chrom_start'} &&
                $other_box->{'chrom_end'} == $box->{'chrom_end'} &&
                ! exists $other_box->{'deleted'}
                ) {
              $past_first_box = 1;
            }
            next;
          };
                                                                                                                                              
          # now start checks for overlaps of boxes
          if (((grep /$protein_id/, @{$other_box->{'ID'}}) ||
	       $other_box->{'chrom_start'} <= $box->{'chrom_end'} && $other_box->{'chrom_end'} >= $box->{'chrom_start'}) &&
              ! exists $other_box->{'deleted'}
              ) {
            $box->{'count'} += $other_box->{'count'};
            $box->{'total_score'} += $other_box->{'total_score'};
            push @{$box->{'ID'}},  @{$other_box->{'ID'}};
            if ($box->{'chrom_start'} > $other_box->{'chrom_start'}) {$box->{'chrom_start'} = $other_box->{'chrom_start'}};
            if ($box->{'chrom_end'} < $other_box->{'chrom_end'}) {$box->{'chrom_end'} = $other_box->{'chrom_end'}};
            # delete the other box
            #print "deleted $other_box->{'type'} $other_box->{'chrom_start'} $other_box->{'chrom_end'} \n";
            $other_box->{'deleted'} = 1; # mark this box as deleted
          }
        }
      }
    }
                                                                                                                                              
                                                                                                                                              
    # no existing boxes found, so start a new one
    if (! $got_a_match) {
      #print "New box for $protein_id, $chrom_start, $chrom_end, $protein_score\n";
      # want to store: $chromosome_start, $chromosome_end, $sense, $protein_id, $protein_score
      my $new_box = {};             # reference to hash
      $new_box->{'sense'} = $sense;
      $new_box->{'chrom_start'} = $chrom_start;
      $new_box->{'chrom_end'} = $chrom_end;
      $new_box->{'total_score'} = $protein_score;
      push @{$new_box->{'ID'}}, ($protein_id);
      $new_box->{'count'} = 1;
      push @boxes, $new_box;
    }
  }
                                                                                                                                              
# now look at each box to see if it has passed our tests
  my $count_of_boxes_passing_tests = 0;
  foreach my $box (@boxes) {
    if (! exists $box->{'deleted'} &&
	$box->{'total_score'} > 100 &&
	$box->{'count'} > 1
	) {
      $count_of_boxes_passing_tests++;
      $box->{'passed_test'} = 1;
      #print "Box passed test with $box->{'chrom_start'}..$box->{'chrom_end'} IDs @{$box->{'ID'}}\n";
    }
  }

# now look to see if we have more than one box and return the details of the IDs in those boxes
  #print "Found $count_of_boxes_passing_tests boxes passing test\n";
  if ($count_of_boxes_passing_tests > 2) {
    foreach my $box (@boxes) {
      if (exists $box->{'passed_test'}) {
	# foreach of the IDs in this box, put the details in @box_result
	foreach my $id (@{$box->{'ID'}}) {
	  foreach my $homology (@peplist) {
	    if ($homology->[0] eq $id) {
	      push @box_result, $homology;
	    }
	  }
	}
      }
    }
  }

  return @box_result
}

##########################################
#
#  &get_EST_split_merged($matched_EST_aref, $chromosome, %Show_in_reverse_orientation);

# find groups of homology that indicate that the genes they match should be split or merged
# look for non-overlapping HSPs that jump back down the position in the protein that is aligned as you move along the chromosome -> split
# look for matches to transcripts with HSPs that do not jump back down the position in the protein that is aligned as you move along the chromosome -> merge

sub get_EST_split_merged {
  my ($matched_EST_aref, $chromosome, %Show_in_reverse_orientation) = @_;

  my @matched = @{$matched_EST_aref};

  my $prev_EST_id = "";	# use to collect alignments for a EST's homology when looking for frameshifts
  my $prev_transcript = "";
  my $prev_hit_start = -1;
  my $prev_hit_end = -1;
  my $prev_chrom_start = -1;
  my $prev_chrom_end = -1;
  my $prev_score = 0;
  my $prev_chrom_strand = "";

  my $got_a_new_transcript;
  my $got_a_big_decrease_in_HSP_start;
  my $got_a_continuation_of_the_HSPs;
 

  # sort the homologies grouped by EST ID and then chromosomal position
  my @homologies = sort {$a->[0] cmp $b->[0] or $a->[1] <=> $b->[1]} @matched;

  foreach my $homology (@homologies) { # $EST_id, $chrom_start, $chrom_end, $chrom_strand, $hit_start, $hit_end, $EST_score, other_data, $matching_transcript

    my $EST_id = $homology->[0];
    my $chrom_start = $homology->[1];
    my $chrom_end = $homology->[2];
    my $chrom_strand = $homology->[3];
    my $hit_start = $homology->[4];
    my $hit_end = $homology->[5];
    my $EST_score = $homology->[6];
    my $matching_transcript = $homology->[8];

    #  check if these are isoforms of the same gene and if so then treat them as the same gene
    if ($matching_transcript =~ /(\S+\.\d+)[a-z]/) {
      $matching_transcript = $1;
    }
    # and change the strange transcript names like H25P06.1.1, H25P06.1.2 to H25P06.1
    if ($matching_transcript =~ /(\S+\.\d+)\.\d+/) {
      $matching_transcript = $1;
    }

    #print "Matched: $EST_id, $chrom_start, $chrom_end, $chrom_strand, $hit_start, $hit_end, $EST_score, $matching_transcript\n";

    # look for homologies which cover two genes or genes which cover a repeated homology
    if ($EST_id eq $prev_EST_id && $chrom_strand eq $prev_chrom_strand) {

      # see if transcript has changed
      if ($prev_transcript ne "" && $prev_transcript ne $matching_transcript) {
	$got_a_new_transcript = 1;
	#print "Merge? prev $prev_transcript this $matching_transcript\n";
      } else {
	$got_a_new_transcript = 0;
      }

      # see if the EST HSP order is unchanged - a possible case for merging if continued over into the next gene
      # same sense, non-overlapping on the chromosome and jump in HSP start
      #
      # NB. The logic here is slightly different to that used when
      # looking for split/merges by protein homology because the 3'
      # EST reads are in the opposite sense (although the sense has
      # been flipped in the routine that reads them in). The EST
      # hit_start and hit_end positions therefore will be going in the
      # wrong direction in those ESTs, so we need to identify them and
      # account for the hit-positions going the 'wrong' way.
      # 
      # We want to avoid picking up the case where there are two duplicated single-exon genes
      # and the homologies have a hit_start that differs by only a few bases - this should be ignored.
      # So we look for cases where the start continues as greater than
      # the previous mid point for evidence of a possible merge.
      #
      if (exists $Show_in_reverse_orientation{$EST_id}) {	# flip the logic of these around
	my $prev_midpoint = $prev_hit_start + ($prev_hit_end - $prev_hit_start)/2;
	my $midpoint = $hit_start + ($hit_end - $hit_start)/2;
	if (($chrom_strand eq '+' && $prev_chrom_end < $chrom_start && $hit_start < $prev_midpoint) ||
	    ($chrom_strand eq '-' && $prev_chrom_end < $chrom_start && $hit_start > $prev_midpoint && $prev_hit_start != -1)) {
	  $got_a_continuation_of_the_HSPs = 1;
	} else {
	  $got_a_continuation_of_the_HSPs = 0;
	}
      
	# see if we have a big increase in the start of the HSP - a possible case for splitting
	if (($chrom_strand eq '+' && $prev_chrom_end < $chrom_start && $prev_hit_start < $hit_start - 100) ||
	    ($chrom_strand eq '-' && $prev_chrom_end < $chrom_start && $prev_hit_start > $hit_start + 100 && $prev_hit_start != -1)) {
	  $got_a_big_decrease_in_HSP_start = 1;
	  #print "Split? $EST_id chrom_start $chrom_end prev_chrom_start $prev_chrom_start $chrom_strand prev_hit_start $prev_hit_start hit_start $hit_start\n";
	} else {
	  $got_a_big_decrease_in_HSP_start = 0;
	}

      } else {			# normal case - just like the proteins
	my $prev_midpoint = $prev_hit_start + ($prev_hit_end - $prev_hit_start)/2;
	my $midpoint = $hit_start + ($hit_end - $hit_start)/2;
	if (($chrom_strand eq '+' && $prev_chrom_end < $chrom_start && $hit_start > $prev_midpoint) ||
	    ($chrom_strand eq '-' && $prev_chrom_end < $chrom_start && $hit_start < $prev_midpoint && $prev_hit_start != -1)) {
	  $got_a_continuation_of_the_HSPs = 1;
	} else {
	  $got_a_continuation_of_the_HSPs = 0;
	}
      
	# see if we have a big decrease in the start of the HSP - a possible case for splitting
	if (($chrom_strand eq '+' && $prev_chrom_end < $chrom_start && $prev_hit_start > $hit_start + 100) ||
	    ($chrom_strand eq '-' && $prev_chrom_end < $chrom_start && $prev_hit_start < $hit_start - 100 && $prev_hit_start != -1)) {
	  $got_a_big_decrease_in_HSP_start = 1;
	  #print "Split? $EST_id chrom_start $chrom_end prev_chrom_start $prev_chrom_start $chrom_strand prev_hit_start $prev_hit_start hit_start $hit_start\n";
	} else {
	  $got_a_big_decrease_in_HSP_start = 0;
	}

      }


      # we have a merge if the pattern of HSPs shows a continuation of
      # the EST matching in order and this is a new gene and the
      # distance between the alignments is less than 10 kb
      if ($got_a_new_transcript && $got_a_continuation_of_the_HSPs && ($chrom_start - $prev_chrom_end) < 10000) {
	# output to database
	# make the anomaly score based on the EST alignment score normalised between 0 and 1
	# the BLAT scores are between 0 and 100
	# use the average of the two alignment's scores to discourage reporting cases where one is very poor
	my $anomaly_score = ($EST_score+$prev_score)/200;
	if ($anomaly_score > 1) {$anomaly_score = 1;}
	if ($anomaly_score < 0) {$anomaly_score = 0;}
	#print "MERGE genes ANOMALY $prev_transcript and $matching_transcript\t$EST_id, $prev_chrom_start, $chrom_end, $chrom_strand, $anomaly_score, 'merge $prev_transcript and $matching_transcript'\n";
	&output_to_database("MERGE_GENES_BY_EST", $chromosome, $EST_id, $prev_chrom_start, $chrom_end, $chrom_strand, $anomaly_score, "merge $prev_transcript and $matching_transcript");


      }

      if ($got_a_big_decrease_in_HSP_start && ! $got_a_new_transcript) {
	# output to database
	# make the anomaly score based on the EST alignment score normalised between 0 and 1
	# the BLAT scores are between 0 and 100
	# use the average of the two alignment's scores to discourage reporting cases where one is very poor
	my $anomaly_score = ($EST_score+$prev_score)/200;
	if ($anomaly_score > 1) {$anomaly_score = 1;}
	if ($anomaly_score < 0) {$anomaly_score = 0;}
	#print "SPLIT gene ANOMALY $matching_transcript\t$EST_id, $prev_chrom_start, $chrom_end, $chrom_strand, $anomaly_score, 'split $matching_transcript'\n";
	&output_to_database("SPLIT_GENES_BY_EST", $chromosome, $EST_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, "split $matching_transcript");
      }
    }

    $prev_transcript = $matching_transcript;
    $prev_hit_start = $hit_start;
    $prev_hit_end = $hit_end;
    $prev_score = $EST_score;
    $prev_chrom_start = $chrom_start;
    $prev_chrom_end = $chrom_end;
    $prev_EST_id = $EST_id;
    $prev_chrom_strand = $chrom_strand;
  }

}

##########################################
# find ESTs that have no association with a transcript in their object
#  &get_unattached_EST(\@est, $chromosome)

sub get_unattached_EST {
  my ($est_aref, $chromosome) = @_;

  my $anomaly_score = 1.0;	# we must look at this with high priority

  foreach my $EST (@$est_aref) { # $id, $start, $end, $sense

    my $EST_id = $EST->[0];
    my $chrom_start = $EST->[1];
    my $chrom_end = $EST->[2];
    my $chrom_strand = $EST->[3];
    #print "EST: $EST_id, $chrom_start, $chrom_end, $chrom_strand\n";

    my $sequence_obj = $ace->fetch(Sequence => $EST_id);

    my @transcripts = $sequence_obj->get('Matching_transcript');
    if ($#transcripts == -1) {
      #print "Not associated with a transcript\n";
      #print "NOT ASSOCIATED: $EST_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score\n";
      &output_to_database("UNATTACHED_EST", $chromosome, $EST_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, 'EST not attached to any transcript');
    }
  }

}



##########################################
# find TSL sites that are defined by an EST but have no association with a transcript in their object
#  &get_unattached_TSL(\@TSL_SL1, \@TSL_SL2, $chromosome)

sub get_unattached_TSL {
  my ($SL1_aref, $SL2_aref, $chromosome) = @_;

  my @SL1 = @$SL1_aref;
  my @SL2 = @$SL2_aref;
  my @SL = sort {$a->[1] <=> $b->[1]} (@SL1, @SL2); # merge and sort the two TSL lists

  my $anomaly_score = 1.0;	# we must look at this with high priority

  foreach my $TSL (@SL) { # $id, $start, $end, $sense

    my $TSL_id = $TSL->[0];
    my $chrom_start = $TSL->[1];
    my $chrom_end = $TSL->[2];
    my $chrom_strand = $TSL->[3];
    #print "TSL: $TSL_id, $chrom_start, $chrom_end, $chrom_strand\n";

    my $feature_obj = $ace->fetch(Feature => $TSL_id);

    my @sequences = $feature_obj->get('Defined_by_sequence');

    # see if we are defined by at least one EST, rather than a TEC-RED
    my $got_an_EST = 0;
    foreach my $sequence (@sequences) {
      #print "> $sequence\n";
      my $seq_obj = $ace->fetch('Sequence' => $sequence);
      if (! $seq_obj->fetch('TSL_tag')) { # ignore TEC-RED sequences
	$got_an_EST = 1;
	last;
      }
    }

    if ($got_an_EST) {
      my @transcripts = $feature_obj->get('Associated_with_transcript');
      if ($#transcripts == -1) {
	#print "Not associated with a transcript\n";
        #print "NOT ASSOCIATED: $TSL_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score\n";
	&output_to_database("UNATTACHED_TSL", $chromosome, $TSL_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, 'Trim CDS to start at the TSL? - also beware of chimeric ESTs.');
      }
    }
  }

}

##########################################
# get TSL sites that do not match a coding transcript or pseudogene
#  &get_isolated_TSL(\@TSL_SL1, \@TSL_SL2, \@coding_transcript, \@transposons, \@transposon_exons, \@noncoding_transcript_exons, $chromosome);

sub get_isolated_TSL {

  my ($SL1_aref, $SL2_aref, $transcripts_aref, $pseudogenes_aref, $transposons_aref, $transposon_exons_aref, $noncoding_transcript_exons_aref, $rRNA_aref, $chromosome) = @_;

  my @SL1 = @$SL1_aref;
  my @SL2 = @$SL2_aref;
  my @SL = sort {$a->[1] <=> $b->[1]} (@SL1, @SL2); # merge and sort the two TSL lists
      
  # allow the TSL to be within 75 bases of the transcript to give a match
  my $transcripts_match = $ovlp->compare($transcripts_aref, near => 75, same_sense => 1);
  my $pseud_match = $ovlp->compare($pseudogenes_aref, same_sense => 1);
  my $trans_match = $ovlp->compare($transposons_aref, same_sense => 1);
  my $trane_match = $ovlp->compare($transposon_exons_aref, same_sense => 1);
  my $nonco_match = $ovlp->compare($noncoding_transcript_exons_aref, same_sense => 1);
  my $rrna_match  = $ovlp->compare($rRNA_aref, same_sense => 1);

  foreach my $tsl (@SL) { # $TSL_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_match = 0;
  
    if ($transcripts_match->match($tsl)) { 
      $got_a_match = 1;
    }

    if ($pseud_match->match($tsl)) { 
      $got_a_match = 1;
    }

    if ($trans_match->match($tsl)) { 
      $got_a_match = 1;
    }

    if ($nonco_match->match($tsl)) { 
      $got_a_match = 1;
    }

    if ($trane_match->match($tsl)) { 
      $got_a_match = 1;
    }

    if ($rrna_match->match($tsl)) { 
      $got_a_match = 1;
    }

    # output unmatched TSL sites to the database
    if (! $got_a_match) {
      my $TSL_id = $tsl->[0];
      my $chrom_start = $tsl->[1];
      my $chrom_end = $tsl->[2];
      my $chrom_strand = $tsl->[3];

      # make the anomaly score 5 because this is very informative
      my $anomaly_score = 5;
      #print "TSL NOT got a match ANOMALY: $TSL_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score\n";
      &output_to_database("UNMATCHED_TSL", $chromosome, $TSL_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, '');
    }
  }

}

##########################################
# get twinscan exons that do not match a coding transcript or pseudogene
# &get_unmatched_twinscan_exons(\@twinscan, \@exons, \@pseudogenes, \@transposons, \@transposon_exons, \@noncoding_transcript_exons, $chromosome)


sub get_unmatched_twinscan_exons {

  my ($twinscan_aref, $exons_aref, $pseudogenes_aref, $transposons_aref, $transposon_exons_aref, $noncoding_transcript_exons_aref, $rRNA_aref, $repeatmasked_aref, $chromosome) = @_;

  my $exons_match = $ovlp->compare($exons_aref, same_sense => 1);
  my $pseud_match = $ovlp->compare($pseudogenes_aref, same_sense => 1);
  my $trans_match = $ovlp->compare($transposons_aref, same_sense => 1);
  my $trane_match = $ovlp->compare($transposon_exons_aref, same_sense => 1);
  my $nonco_match = $ovlp->compare($noncoding_transcript_exons_aref, same_sense => 1);
  my $rrna_match  = $ovlp->compare($rRNA_aref, same_sense => 1);
  my $repeat_match= $ovlp->compare($repeatmasked_aref, near => -20);

  foreach my $twinscan (@{$twinscan_aref}) { # $twinscan_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_match = 0;
  
    if ($exons_match->match($twinscan)) { #&match($twinscan, $exons_aref, \%exons_match)) {
      $got_a_match = 1;
    }

    if ($pseud_match->match($twinscan)) { #&match($twinscan, $pseudogenes_aref, \%pseud_match)) {
      $got_a_match = 1;
    }

    if ($trans_match->match($twinscan)) { #&match($twinscan, $transposons_aref, \%trans_match)) {
      $got_a_match = 1;
    }

    if ($trane_match->match($twinscan)) { #&match($twinscan, $transposon_exons_aref, \%trane_match)) {
      $got_a_match = 1;
    }

    if ($nonco_match->match($twinscan)) { #&match($twinscan, $noncoding_transcript_exons_aref, \%nonco_match)) {
      $got_a_match = 1;
    }

    if ($rrna_match->match($twinscan)) { #&match($twinscan, $rRNA_aref, \%rrna_match)) {
      $got_a_match = 1;
    }

    # don't want to report unmatched twinscan that overlaps with a repeat
    if ($repeat_match->match($twinscan)) { #&match($twinscan, $repeatmasked_aref, \%repeat_match)) { 
      $got_a_match = 1;
    }

    # output unmatched TWINSCAN sites to the database
    if (! $got_a_match) {
      my $TWINSCAN_id = $twinscan->[0];
      my $chrom_start = $twinscan->[1];
      my $chrom_end = $twinscan->[2];
      my $chrom_strand = $twinscan->[3];

      # make the anomaly score 1
      my $anomaly_score = 1.0;
      #print "TWINSCAN NOT got a match ANOMALY: $TWINSCAN_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score\n";
      &output_to_database("UNMATCHED_TWINSCAN", $chromosome, $TWINSCAN_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, '');
    }
  }

}


##########################################
# get genefinder exons that do not match a coding transcript or pseudogene
# &get_unmatched_genefinder_exons(\@genefinder, \@exons, \@pseudogenes, \@transposons, \@transposon_exons, \@noncoding_transcript_exons, $chromosome)


sub get_unmatched_genefinder_exons {

  my ($genefinder_aref, $exons_aref, $pseudogenes_aref, $transposons_aref, $transposon_exons_aref, $noncoding_transcript_exons_aref, $rRNA_aref, $repeatmasked_aref, $chromosome) = @_;

  my $exons_match = $ovlp->compare($exons_aref, same_sense => 1);
  my $pseud_match = $ovlp->compare($pseudogenes_aref, same_sense => 1);
  my $trans_match = $ovlp->compare($transposons_aref, same_sense => 1);
  my $trane_match = $ovlp->compare($transposon_exons_aref, same_sense => 1);
  my $nonco_match = $ovlp->compare($noncoding_transcript_exons_aref, same_sense => 1);
  my $rrna_match  = $ovlp->compare($rRNA_aref, same_sense => 1);
  my $repeat_match= $ovlp->compare($repeatmasked_aref, near => -20);	# allow 20 bases overlap in either sense with a repeat before we call it a hit

  foreach my $genefinder (@{$genefinder_aref}) { # $genefinder_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_match = 0;
  
    if ($exons_match->match($genefinder)) { #&match($genefinder, $exons_aref, \%exons_match)) {
      $got_a_match = 1;
    }

    if ($pseud_match->match($genefinder)) { #&match($genefinder, $pseudogenes_aref, \%pseud_match)) {
      $got_a_match = 1;
    }

    if ($trans_match->match($genefinder)) { #&match($genefinder, $transposons_aref, \%trans_match)) {
      $got_a_match = 1;
    }

    if ($trane_match->match($genefinder)) { #&match($genefinder, $transposon_exons_aref, \%trane_match)) {
      $got_a_match = 1;
    }

    if ($nonco_match->match($genefinder)) { #&match($genefinder, $noncoding_transcript_exons_aref, \%nonco_match)) {
      $got_a_match = 1;
    }

    if ($rrna_match->match($genefinder)) { #&match($genefinder, $rRNA_aref, \%rrna_match)) {
      $got_a_match = 1;
    }

    # don't want to report unmatched genefinder that overlaps with a repeat
    if ($repeat_match->match($genefinder)) { #&match($genefinder, $repeatmasked_aref, \%repeat_match)) {
      $got_a_match = 1;
    }

    # output unmatched GENEFINDER sites to the database
    if (! $got_a_match) {
      my $GENEFINDER_id = $genefinder->[0];
      my $chrom_start = $genefinder->[1];
      my $chrom_end = $genefinder->[2];
      my $chrom_strand = $genefinder->[3];

      # make the anomaly score 1
      my $anomaly_score = 1.0;
      #print "GENEFINDER NOT got a match ANOMALY: $GENEFINDER_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score\n";
      &output_to_database("UNMATCHED_GENEFINDER", $chromosome, $GENEFINDER_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, '');
    }
  }

}

##########################################
# find the SAGE transcripts that don't overlap the coding transcripts and pseudogenes
#  &get_unmatched_SAGE(\@coding_transcripts, \@SAGE_tags);

sub get_unmatched_SAGE {

  my ($coding_transcripts_aref, $pseudogenes_aref, $SAGE_tags_aref, $transposons_aref, $transposon_exons_aref, $noncoding_transcript_exons_aref, $rRNA_aref, $miRNA_aref, $ncRNA_aref, $scRNA_aref, $snRNA_aref, $snoRNA_aref, $stRNA_aref, $tRNA_aref, $tRNAscan_SE_1_23RNA_aref, $chromosome) = @_;
 
  my @SAGE_tags = @{$SAGE_tags_aref};

  # count a hit within 20 bases as a match, so tags have to be more
  # than 20bp away before they are counted as unmatched
  my $NEAR = 20;		

  my $coding_transcripts_match = $ovlp->compare($coding_transcripts_aref, near => $NEAR);
  my $pseud_match = $ovlp->compare($pseudogenes_aref, near => $NEAR);
  my $trans_match = $ovlp->compare($transposons_aref, near => $NEAR);
  my $trane_match = $ovlp->compare($transposon_exons_aref, near => $NEAR);
  my $nonco_match = $ovlp->compare($noncoding_transcript_exons_aref, near => $NEAR);
  my $rrna_match = $ovlp->compare($rRNA_aref, near => $NEAR);
  my $mirna_match = $ovlp->compare($miRNA_aref, near => $NEAR);
  my $ncrna_match = $ovlp->compare($ncRNA_aref, near => $NEAR);
  my $scrna_match = $ovlp->compare($scRNA_aref, near => $NEAR);
  my $snrna_match = $ovlp->compare($snRNA_aref, near => $NEAR);
  my $snorna_match = $ovlp->compare($snoRNA_aref, near => $NEAR);
  my $strna_match = $ovlp->compare($stRNA_aref, near => $NEAR);
  my $trna_match = $ovlp->compare($tRNA_aref, near => $NEAR);
  my $tRNAscan_SE_1_23rna_match = $ovlp->compare($tRNAscan_SE_1_23RNA_aref, near => $NEAR);

  foreach my $sage (@SAGE_tags) { # $SAGE_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_match = 0;
  
    if ($coding_transcripts_match->match($sage)) { #&match($sage, $coding_transcripts_aref, \%coding_transcripts_match)) {
      $got_a_match = 1;
    }

    if ($pseud_match->match($sage)) { #&match($sage, $pseudogenes_aref, \%pseud_match)) {
      $got_a_match = 1;
    }

    if ($trans_match->match($sage)) { #&match($sage, $transposons_aref, \%trans_match)) {
      $got_a_match = 1;
    }

    if ($trane_match->match($sage)) { #&match($sage, $transposon_exons_aref, \%trane_match)) {
      $got_a_match = 1;
    }

    if ($nonco_match->match($sage)) { #&match($sage, $noncoding_transcript_exons_aref, \%nonco_match)) {
      $got_a_match = 1;
    }

    if ($rrna_match->match($sage)) { #&match($sage, $rRNA_aref, \%rrna_match)) {
      $got_a_match = 1;
    }

    if ($mirna_match->match($sage)) { #&match($sage, $miRNA_aref, \%mirna_match)) {
      $got_a_match = 1;
    }

    if ($ncrna_match->match($sage)) { #&match($sage, $ncRNA_aref, \%ncrna_match)) {
      $got_a_match = 1;
    }

    if ($scrna_match->match($sage)) { #&match($sage, $scRNA_aref, \%scrna_match)) {
      $got_a_match = 1;
    }

    if ($snrna_match->match($sage)) { #&match($sage, $snRNA_aref, \%snrna_match)) {
      $got_a_match = 1;
    }

    if ($snorna_match->match($sage)) { #&match($sage, $snoRNA_aref, \%snorna_match)) {
      $got_a_match = 1;
    }

    if ($strna_match->match($sage)) { #&match($sage, $stRNA_aref, \%strna_match)) {
      $got_a_match = 1;
    }

    if ($trna_match->match($sage)) { #&match($sage, $tRNA_aref, \%trna_match)) {
      $got_a_match = 1;
    }

    if ($tRNAscan_SE_1_23rna_match->match($sage)) { #&match($sage, $tRNAscan_SE_1_23RNA_aref, \%tRNAscan_SE_1_23rna_match)) {
      $got_a_match = 1;
    }

    # output unmatched SAGE sites to the database
    if (! $got_a_match) {
      my $SAGE_id = $sage->[0];
      my $chrom_start = $sage->[1];
      my $chrom_end = $sage->[2];
      my $chrom_strand = $sage->[3];
      my $other_data = $sage->[4];
      my ($score) = ($other_data =~ /count\s+(\d+)/); # get the score from the 'count' field of the 'other' data

      my $anomaly_score = $score/50;	# the SAGE_tag scores are between 1 and about 100, so get an anomaly score of from 0 to about 2
      if ($anomaly_score > 2.5) {$anomaly_score = 2.5;}	# cap the score so we don't get extraordinarily large ones

      #print "NOT got a match ANOMALY: $SAGE_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score\n";
      &output_to_database("UNMATCHED_SAGE", $chromosome, $SAGE_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, '');
    }
  }

}



##########################################
# get genefinder exons that do not match a coding transcript or pseudogene
#  &get_unmatched_waba_coding(\@waba_coding, \@exons, \@pseudogenes, \@transposons, \@transposon_exons, \@noncoding_transcript_exons, \@rRNA, \@repeatmasked, $chromosome);



sub get_unmatched_waba_coding {

  my ($waba_aref, $exons_aref, $pseudogenes_aref, $transposons_aref, $transposon_exons_aref, $noncoding_transcript_exons_aref, $rRNA_aref, $repeatmasked_aref, $chromosome) = @_;

  my $exons_match = $ovlp->compare($exons_aref);
  my $pseud_match = $ovlp->compare($pseudogenes_aref);
  my $trans_match = $ovlp->compare($transposons_aref);
  my $trane_match = $ovlp->compare($transposon_exons_aref);
  my $nonco_match = $ovlp->compare($noncoding_transcript_exons_aref);
  my $rrna_match  = $ovlp->compare($rRNA_aref);
  my $repeat_match= $ovlp->compare($repeatmasked_aref);

  foreach my $waba (@{$waba_aref}) { # $waba_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_match = 0;
  
    if ($exons_match->match($waba)) { #&match($waba, $exons_aref, \%exons_match)) {
      $got_a_match = 1;
    }

    if ($pseud_match->match($waba)) { #&match($waba, $pseudogenes_aref, \%pseud_match)) {
      $got_a_match = 1;
    }

    if ($trans_match->match($waba)) { #&match($waba, $transposons_aref, \%trans_match)) {
      $got_a_match = 1;
    }

    if ($trane_match->match($waba)) { #&match($waba, $transposon_exons_aref, \%trane_match)) {
      $got_a_match = 1;
    }

    if ($nonco_match->match($waba)) { #&match($waba, $noncoding_transcript_exons_aref, \%nonco_match)) {
      $got_a_match = 1;
    }

    if ($rrna_match->match($waba)) { #&match($waba, $rRNA_aref, \%rrna_match)) {
      $got_a_match = 1;
    }

    if ($repeat_match->match($waba)) { #&match($waba, $repeatmasked_aref, \%repeat_match)) {
      $got_a_match = 1;
    }

    # output unmatched WABA sites to the database
    if (! $got_a_match) {
      my $WABA_id = $waba->[0];
      my $chrom_start = $waba->[1];
      my $chrom_end = $waba->[2];
      my $chrom_strand = $waba->[3];
      my $waba_score = $waba->[6];

      # reject any hits where either of the proteins have a Blast score < 50
      if ($waba_score < 50) {next;}

      # make the anomaly score based on the protein alignment score normalised between 1 and about 2.5
      # using the log10 of the waba score
      # the waba scores seem to be between 0 and about 250
      my $anomaly_score = POSIX::log10($waba_score);
      if ($anomaly_score > 2.5) {$anomaly_score = 2.5;}
      if ($anomaly_score < 0) {$anomaly_score = 0;}

      #print "WABA NOT got a match ANOMALY: $WABA_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score\n";
      &output_to_database("UNMATCHED_WABA", $chromosome, $WABA_id, $chrom_start, $chrom_end, '.', $anomaly_score, '');
    }
  }

}
##########################################
# get CDS exons that overlap a CDS exon on the opposite sense
#  &get_matched_exons(\@cds_exons, \@pseudogenes, \@transposons, \@transposon_exons, \@noncoding_transcript_exons, \@rRNA, $chromosome);



sub get_matched_exons {

  my ($cds_exons_aref, $pseudogenes_aref, $transposons_aref, $transposon_exons_aref, $noncoding_transcript_exons_aref, $rRNA_aref, $chromosome) = @_;

  my $exons_match = $ovlp->compare($cds_exons_aref, other_sense => 1);
  my $pseud_match = $ovlp->compare($pseudogenes_aref, other_sense => 1);
  my $trans_match = $ovlp->compare($transposons_aref, other_sense => 1);
  my $trane_match = $ovlp->compare($transposon_exons_aref, other_sense => 1);
  my $nonco_match = $ovlp->compare($noncoding_transcript_exons_aref, other_sense => 1);
  my $rrna_match  = $ovlp->compare($rRNA_aref, other_sense => 1);

  foreach my $exon (@{$cds_exons_aref}) { # $exon_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_match = 0;
    my $matching_thing;
    my @result;
    my @names;

    if (@result = $exons_match->match($exon)) { #&match($exon, $cds_exons_aref, \%exons_match)) { # look for matches of exons against exons
      $got_a_match = 1;
      @names = $exons_match->matching_IDs;
      $matching_thing = "Overlaps gene: @names";
    }

    if (@result = $pseud_match->match($exon)) { #&match($exon, $pseudogenes_aref, \%pseud_match)) {
      $got_a_match = 1;
      @names = $pseud_match->matching_IDs;
      $matching_thing = "Overlaps pseudogene: @names";
    }

    if (@result = $trans_match->match($exon)) { #&match($exon, $transposons_aref, \%trans_match)) {
      $got_a_match = 1;
      @names = $trans_match->matching_IDs;
      $matching_thing = "Overlaps transposon: @names";
    }

    if (@result = $trane_match->match($exon)) { #&match($exon, $transposon_exons_aref, \%trane_match)) {
      $got_a_match = 1;
      @names = $trane_match->matching_IDs;
      $matching_thing = "Overlaps transposon: @names";
    }

    if (@result = $nonco_match->match($exon)) { #&match($exon, $noncoding_transcript_exons_aref, \%nonco_match)) {
      $got_a_match = 1;
      @names = $nonco_match->matching_IDs;
      $matching_thing = "Overlaps non-coding transcript: @names";
    }

    if (@result = $rrna_match->match($exon)) { #&match($exon, $rRNA_aref, \%rrna_match)) {
      $got_a_match = 1;
      @names = $rrna_match->matching_IDs;
      $matching_thing = "Overlaps rRNA: @names";
    }

    # output matched EXON sites to the database
    if ($got_a_match) {
      my $exon_id = $exon->[0];
      my $chrom_start = $exon->[1];
      my $chrom_end = $exon->[2];
      my $chrom_strand = $exon->[3];
      my $exon_score = $exon->[6];

      my $anomaly_score = 5;

      #print "EXON overlapping other thing ANOMALY: $exon_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score\n";
      &output_to_database("OVERLAPPING_EXONS", $chromosome, $exon_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, $matching_thing);
    }
  }

}

##########################################
# Finding short CDS exons
# note any exons shorter than 30 bases
#  &get_short_exons(\@exons, $chromosome);

sub get_short_exons {
  my ($exons_aref, $chromosome) = @_;

  foreach my $exon (@{$exons_aref}) { # $exon_id, $chrom_start, $chrom_end, $chrom_strand

    my $exon_id = $exon->[0];
    my $chrom_start = $exon->[1];
    my $chrom_end = $exon->[2];
    my $chrom_strand = $exon->[3];
    my $exon_score = $exon->[6];

    my $anomaly_score = 1;

    if ($chrom_end - $chrom_start + 1< 30) { # note any exons shorter than 30 bases
      &output_to_database("SHORT_EXON", $chromosome, $exon_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, '');
    }
  }
}

##########################################
# Finding EST/genome mismatches
# note any 1 to 3 base misalignment that is seen in two or more ESTs
# this might indicate a genomic sequencing error or some RNA editing
# &get_est_mismatches(\@est, $chromosome);

sub get_est_mismatches {
  my ($est_aref, $chromosome) = @_;
  my @ests = @{$est_aref};

  my $prev_EST_id = "";	
  my $prev_chrom_start = -1;
  my $prev_chrom_end = -1;
  my $prev_score = 0;
  my $prev_chrom_strand = "";
  my $prev_hit_start = -1;
  my $prev_hit_end = -1;


  
  # sort the homologies grouped by EST ID and then chromosomal position
  my @sorted_ests = sort {$a->[0] cmp $b->[0] or $a->[1] <=> $b->[1]} @ests;

  my %mismatch_counts;
  foreach my $alignment (@sorted_ests) { # $EST_id, $chrom_start, $chrom_end, $chrom_strand, $hit_start, $hit_end, $EST_score, $matching_transcript

    my $EST_id = $alignment->[0];
    my $chrom_start = $alignment->[1];
    my $chrom_end = $alignment->[2];
    my $chrom_strand = $alignment->[3];
    my $hit_start = $alignment->[4];
    my $hit_end = $alignment->[5];
    my $EST_score = $alignment->[6];

    if ($prev_EST_id eq $EST_id && # the same EST as the previous hit
	$EST_score > 99.0 &&	# good quality
	$chrom_end - $chrom_start > 30 && # a decent length
	$prev_chrom_end - $prev_chrom_start > 30) { # the previous hit was a decent length
      if ((($hit_start - $prev_hit_end == 1) && ($chrom_start - $prev_chrom_end) < 4) ||
	  (($hit_start - $prev_hit_end < 4) && ($chrom_start - $prev_chrom_end) == 1)) {
	my $key = "$chrom_start-$chrom_end";
	$mismatch_counts{$key}{count}++;
	$mismatch_counts{$key}{EST_id} = $EST_id;
	$mismatch_counts{$key}{chrom_start} = $chrom_start;
	$mismatch_counts{$key}{chrom_end} = $chrom_end;
	$mismatch_counts{$key}{chrom_strand} = $chrom_strand;
	$mismatch_counts{$key}{overlaps} = 0;
      }
    }

    $prev_EST_id = $EST_id;
    $prev_chrom_start = $chrom_start;
    $prev_chrom_end = $chrom_end;
    $prev_chrom_strand = $chrom_strand;
    $prev_hit_start = $hit_start;
    $prev_hit_end = $hit_end;
    $prev_score = $EST_score;
  }

  # we now have the set of small misalignments in %mismatch_counts

  # get rid of the singleton misalignments
  foreach my $key (keys %mismatch_counts) {
    if ($mismatch_counts{$key}{count} == 1) {
      delete $mismatch_counts{$key};
    }
  }


  # now count how many good alignments overlap with them
  foreach my $alignment (@sorted_ests) { # $EST_id, $chrom_start, $chrom_end, $chrom_strand, $hit_start, $hit_end, $EST_score, $matching_transcript

    my $chrom_start = $alignment->[1];
    my $chrom_end = $alignment->[2];
    my $chrom_strand = $alignment->[3];

    foreach my $key (keys %mismatch_counts) {
      if ($chrom_start < $mismatch_counts{$key}{chrom_start} && 
	  $chrom_end > $mismatch_counts{$key}{chrom_start} &&
	  $chrom_strand eq $mismatch_counts{$key}{chrom_strand}) {
	$mismatch_counts{$key}{overlaps}++;
      }
    }
  }

  foreach my $key (keys %mismatch_counts) {
    if ($mismatch_counts{$key}{count} > 1) { # this should now always be true
      # get score as proportion of ESTs there that have the misalignment
      my $anomaly_score = $mismatch_counts{$key}{count}/($mismatch_counts{$key}{overlaps} + $mismatch_counts{$key}{count});
      if ($anomaly_score >= 0.5) { # we don't want the dubious cases
	#print "MISMATCHED_EST ", $chromosome, " ID ", $mismatch_counts{$key}{EST_id}, " start ", $mismatch_counts{$key}{chrom_start}, " end ", $mismatch_counts{$key}{chrom_end}, " strand ", $mismatch_counts{$key}{chrom_strand}, "score $anomaly_score\n";
	&output_to_database("MISMATCHED_EST", $chromosome, $mismatch_counts{$key}{EST_id}, $mismatch_counts{$key}{chrom_start}, $mismatch_counts{$key}{chrom_end}, $mismatch_counts{$key}{chrom_strand}, $anomaly_score, 'possible genomic sequence error or evidence of RNA editing');
      }
    }
  }
}
##########################################
# Finding weak exon splice sites
# output to database any splice sites with a score < 0.1
# output to database any splice sites with a score < 0.5 if there are two or more such splice sites in a gene
# get_weak_exon_splice_sites(\@exons, $chromosome);

sub get_weak_exon_splice_sites {
  my ($CDS_introns_aref, $chromosome) = @_;


# intron splice sites with key = "ID,strand,pos" value as list 
# value[0] = 0 if not confirmed, > 0 if confirmed 
# value[1] is CDS ID
# value[2] is strand
# value[3] is chromosomal pos
# value[4] is 5 or 3 (for the 5' or 3' end)

  my %splice_sites;

  ###########################################
  # find out which splice sites are confirmed
  ###########################################

  foreach my $intron (@{$CDS_introns_aref}) { # $intron_id, $chrom_start, $chrom_end, $chrom_strand, other

    my $CDS_id = $intron->[0];
    my $chrom_start = $intron->[1];
    my $chrom_end = $intron->[2];
    my $chrom_strand = $intron->[3];
    my $other_data = $intron->[4];
    my $confirmed = $other_data =~ /Confirmed_/; # test to see if this intron is confirmed

    # 5' splice site
    if (exists $splice_sites{"$CDS_id$chrom_strand$chrom_start"}) {
      $splice_sites{"$CDS_id$chrom_strand$chrom_start"}[0] += $confirmed;
    } else {
      $splice_sites{"$CDS_id$chrom_strand$chrom_start"} = [$confirmed, $CDS_id, $chrom_strand, $chrom_start-1, 5]; # we want the base before the cut
    }

    # 3' splice site
    if (exists $splice_sites{"$CDS_id$chrom_strand$chrom_end"}) {
      $splice_sites{"$CDS_id$chrom_strand$chrom_end"}[0] += $confirmed;
    } else {
      $splice_sites{"$CDS_id$chrom_strand$chrom_end"} = [$confirmed, $CDS_id, $chrom_strand, $chrom_end, 3];       
    }
  }

  my @splices;			# list of unconfirmed splice sites
   

  #########################################################
  # get a list of unconfirmed splice sites sorted by CDS_id
  #########################################################

  foreach my $splice (keys %splice_sites) {
    my $confirmed  = $splice_sites{$splice}->[0];
    my $CDS_id     = $splice_sites{$splice}->[1];
    my $strand     = $splice_sites{$splice}->[2];
    my $splice_pos = $splice_sites{$splice}->[3];
    my $splice_end = $splice_sites{$splice}->[4];
    if ($confirmed == 0) {
      push @splices, [$CDS_id, $strand, $splice_pos, $splice_end];
    }
  }
  print "Have ", scalar (keys %splice_sites), " splice sites, of which ", scalar @splices," are not confirmed\n";
  # sort splices by CDS_id so we can find out how many weak splices there are in a CDS
  my @splices_sorted = sort {$a->[0] cmp $b->[0] } @splices;


  #################################################################
  # foreach splice site, see if it is below the cutoff or
  # cutoff_single thresholds output it if it is below the
  # cutoff_single threshold and output if two or more are below the
  # cutoff threshold (store any that might need to be output)
  #################################################################

  my $pwm = PWM->new;
  my $seq_file = "$database/CHROMOSOMES/${\$wormbase->chromosome_prefix}$chromosome.dna";
  my $seq = read_file($seq_file);

  my $splice_cutoff = 0.5;	# if we have two or more splice sites below this value, then report them
  my $splice_cutoff_single = 0.1; # if we have only one splice site, it has to be below this value before we report it
  my $count_low_sites = 0;

  my $anomaly_score = 1;

  my $prev_CDS = "";
  my @prev_details = ();	# details to store in database of any previous splice site in this CDS with a score < $splice_cutoff

  foreach my $splice (@splices_sorted) {
    my $CDS_id     = $splice->[0];
    my $strand     = $splice->[1];
    my $splice_pos = $splice->[2];
    my $splice_end = $splice->[3];

    if ($prev_CDS ne $CDS_id) {	# we have a new CDS, so reset the count of low-scoring sites in this CDS
      $prev_CDS = $CDS_id;
      $count_low_sites = 0;
      @prev_details = ();
    }


    my $score;			# score of this splice site
    if ($strand eq '+') {
      if ($splice_end == 5) {	# is it a 5' end of the intron
	$score = $pwm->splice5($seq, $splice_pos-1, '+'); # -1 to convert from acedb sequence pos to perl string coords
      } else {			# or the 3' end of the intron
	$score = $pwm->splice3($seq, $splice_pos-1, '+'); # -1 to convert from acedb sequence pos to perl string coords	
      }

    } else {			# strand eq '-'
      if ($splice_end == 3) {	# is it a 5' end of the intron (yes, we test for 3 here to get a real 5' in the reverse sense!)
	$score = $pwm->splice5($seq, $splice_pos, '-'); # -1 to convert from acedb sequence pos to perl string coords
      } else {                  # or the 3' end of the intron
	$score = $pwm->splice3($seq, $splice_pos, '-'); # -1 to convert from acedb sequence pos to perl string coords
      } 
    }

    if ($score < $splice_cutoff) {
      if ($score < $splice_cutoff_single || $count_low_sites > 0) {
	&output_to_database("WEAK_INTRON_SPLICE_SITE", $chromosome, $CDS_id, $splice_pos-1, $splice_pos+1, $strand, $anomaly_score, '');
	if ($count_low_sites == 1 && @prev_details) { # output the stored details of the first site
	  &output_to_database("WEAK_INTRON_SPLICE_SITE", @prev_details);
	}
      } else {		# store the details in case we find a second low-scoring site in this CDS
	@prev_details = ($chromosome, $CDS_id, $splice_pos-1, $splice_pos+1, $strand, $anomaly_score, '');
      }
      $count_low_sites++;
    }

  }				# foreach $splice
}

##########################################
# read file

sub read_file {
  my ($file) = @_;

  $/ = "";
  open (SEQ, $file) or die "Can't open the dna file for $file : $!\n";
  my $seq = <SEQ>;
  close SEQ;
  $/ = "\n";

  $seq =~ s/^>.*\n//;		# remove one initial title line
  $seq =~ s/\n//g;

  return $seq
}
##########################################
# Finding short CDS introns
# note any introns shorter than 30 bases
#  &get_short_introns(\@exons, $chromosome);

sub get_short_introns {
  my ($exons_aref, $chromosome) = @_;

  my $prev_end = 0;
  my $prev_exon = "";

  foreach my $exon (@{$exons_aref}) { # $exon_id, $chrom_start, $chrom_end, $chrom_strand

    my $exon_id = $exon->[0];
    my $chrom_start = $exon->[1];
    my $chrom_end = $exon->[2];
    my $chrom_strand = $exon->[3];
    my $exon_score = $exon->[6];

    my $anomaly_score = 1;

    # we want: 
    # not the first exon in the chromosome
    # and the previous exon came from this CDS
    # and the intron is < 30 bases
    if ($prev_end != 0 && $prev_exon eq $exon_id && $chrom_start - $prev_end + 1< 30) { # note any introns shorter than 30 bases
      &output_to_database("SHORT_INTRON", $chromosome, $exon_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, '');
    }
    $prev_end = $chrom_end;
    $prev_exon = $exon_id;
  }
}

##########################################
# find coding exons that match repeatmasked regions 
#  &get_matched_repeatmasker(\@exons, \@repeatmasked, $chromosome);


sub get_matched_repeatmasker {

  my ($exons_aref, $repeatmasked_aref, $chromosome) = @_;


  # remove the low-complexity repeat regions from the repeats list

  my $repeat_match = $ovlp->compare($repeatmasked_aref, near => -20);
  my @results;
  my @names;

  foreach my $exon (@{$exons_aref}) { # $exon_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_match = 0;
  
    if (@results = $repeat_match->match($exon)) { #&match($exon, $repeatmasked_aref, \%repeat_match)) {
      $got_a_match = 1;
      @names = $repeat_match->matching_IDs;
      # +++ print "REPEATS: @names\n";
    }

    # output matched exons to the database
    if ($got_a_match) {
      my $exon_id = $exon->[0];
      my $chrom_start = $exon->[1];
      my $chrom_end = $exon->[2];
      my $chrom_strand = $exon->[3];
      my $exon_score = $exon->[6];

      my $anomaly_score = 1;

      #print "REPEAT got a match ANOMALY: $exon_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score\n";
      &output_to_database("REPEAT_OVERLAPS_EXON", $chromosome, $exon_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, '');
    }
  }

}

##########################################
# get the UTRs which have got 3 or more exons
sub get_multiple_utr_introns {

  my ($UTR_aref, $chromosome) = @_;


  my %output;
  my %output_start;
  my %output_end;
  my %output_strand;

  my $count;
  my $prev = "";

  foreach my $utr (@{$UTR_aref}) { # $id, $chrom_start, $chrom_end, $chrom_strand

    # output multiple UTR exons to the database
    my $id = $utr->[0];
    my $chrom_start = $utr->[1];
    my $chrom_end = $utr->[2];
    my $chrom_strand = $utr->[3];

    # we don't want to see UTRs from isoform transcripts - there are too many mangled UTRs in them
    if ($id !~ /^\S+\.\d+$/ || $id =~ /^\S+\.\d+[a-z]\.\d+$/) {next;} 

    if ($prev ne $id) {
      $prev = $id;
      $count = 0;
    }
    $count++;
    if ($count > 2) {
      $output{$id} = $count;
      if ($count == 3) {$output_start{$id} = $chrom_start;}
      $output_end{$id} = $chrom_end;
      $output_strand{$id} = $chrom_strand;
    }
  }

  # now output the ones we have found
  foreach my $id (keys %output) {
    &output_to_database("INTRONS_IN_UTR", $chromosome, $id, $output_start{$id}, $output_end{$id}, $output_strand{$id}, $output{$id}, '');
  }

}

##########################################
# get genes to be split/merged based on twinscan

sub get_twinscan_split_merged {

  my ($twinscan_aref, $CDS_aref, $chromosome) = @_;


  my $twin_match = $ovlp->compare($twinscan_aref, same_sense => 1);
  my @results;
  my @names;

  foreach my $CDS (@{$CDS_aref}) { # $CDS_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_match = 0;
    my $matching_ids;

    if (@results = $twin_match->match($CDS)) { # &match($CDS, $twinscan_aref, \%twin_match)) {
      @names = $twin_match->matching_IDs;
      #print "SPLIT $CDS->[0] @names\n";
      if (@results > 1) {$got_a_match = 1;}
    }

    # output matched CDS to the database
    if ($got_a_match) {
      my $CDS_id = $CDS->[0];
      my $chrom_start = $CDS->[1];
      my $chrom_end = $CDS->[2];
      my $chrom_strand = $CDS->[3];

      my $anomaly_score = 1;

      #print "SPLIT_GENE_BY_TWINSCAN, $chromosome, $CDS_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, Twinscan: @names\n";
      &output_to_database("SPLIT_GENE_BY_TWINSCAN", $chromosome, $CDS_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, "Twinscan: @names");
    }
  }

  my $CDS_match = $ovlp->compare($CDS_aref, same_sense => 1);

  foreach my $twin (@{$twinscan_aref}) { # $twin_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_merge = 0;
    my $matching_ids;

    if (@results = $CDS_match->match($twin)) { #&match($twin, $CDS_aref, \%CDS_match)) {
      @names = $CDS_match->matching_IDs;
      #  check if these are isoforms of the same gene and if so then treat them as the same gene
      my $prev_name = "";

      foreach my $name (@names) {
	my $simple_name = $name;
	if ($name =~ /(\S+\.\d+)[a-z]/) {
	  $simple_name = $1;
	}
	#if (! defined $simple_name) {print "simple name not defined with name = $name\n";}
	if ($prev_name ne "" && $simple_name ne $prev_name) {
	  $got_a_merge = 1;
	}
	$prev_name = $simple_name;
      }

      #print "MERGE $twin->[0] @names\n";
    }

    # output matched twinscans to the database
    if ($got_a_merge) {
      my $twin_id = $twin->[0];
      my $chrom_start = $twin->[1];
      my $chrom_end = $twin->[2];
      my $chrom_strand = $twin->[3];

      my $anomaly_score = 1;

      #print "MERGE_GENES_BY_TWINSCAN, $chromosome, $twin_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, CDS: @names\n";
      &output_to_database("MERGE_GENES_BY_TWINSCAN", $chromosome, $twin_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, "CDS: @names");
    }
  }

}

####################################################################################
sub get_checked_confirmed_introns {
  my ($check_introns_EST_aref, $check_introns_cDNA_aref, $chromosome) = @_;

  foreach my $check (@{$check_introns_EST_aref}, @{$check_introns_cDNA_aref}) {
    my $anomaly_score = 10;	# high priority
    my $check_id     = $check->[0];
    my $chrom_start  = $check->[1];
    my $chrom_end    = $check->[2];
    my $chrom_strand = $check->[3];


    &output_to_database("CONFIRMED_INTRON", $chromosome, $check_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, '');
  }

}

##########################################
# get ESTs not matching transcripts

sub get_unmatched_ests {

  my ($est_aref, $trans_aref, $pseud_aref, $tposon_aref, $nctrans_aref, $rrna_aref, $chromosome) = @_;

  # set up the overlap compare objects for the secondary lists
  my $trans_match = $ovlp->compare($trans_aref);
  my $pseud_match = $ovlp->compare($pseud_aref);
  my $tposon_match = $ovlp->compare($tposon_aref);
  my $nctrans_match = $ovlp->compare($nctrans_aref);
  my $rrna_match = $ovlp->compare($rrna_aref);

  foreach my $est (@{$est_aref}) { # $id, $chrom_start, $chrom_end, $chrom_strand
 
    my @matching = $trans_match->match($est);
    push @matching, $pseud_match->match($est);
    push @matching, $tposon_match->match($est);
    push @matching, $nctrans_match->match($est);
    push @matching, $rrna_match->match($est);

    if (@matching == 0) {

    # output unmatched ests to the database
      my $id = $est->[0];
      my $chrom_start = $est->[1];
      my $chrom_end = $est->[2];
      my $chrom_strand = $est->[3];
      if ($chrom_end - $chrom_start < 50) {next;}

      my $anomaly_score = 1;

      &output_to_database("UNMATCHED_EST", $chromosome, $id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, "");
    }
  }
}


####################################################################################
####################################################################################
####################################################################################
####################################################################################

##########################################
# checks if there is sense = + or - and if not, outputs two anomaly records, one for each sense.
# else it just outputs one record for the specified sense.

sub output_to_database {
  my ($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, $explanation ) = @_;

  # get the clone and lab for this location
  my ($clone, $clone_start, $clone_end) = $coords->LocateSpan("${\$wormbase->chromosome_prefix}$chromosome", $chrom_start, $chrom_end);
  my $lab =  $clonelab{$clone};          # get the lab that sequenced this clone  
  if (! defined $lab) {
    if ($clone =~ /SUPERLINK_CB_/) {
      $lab = 'HX';
    } else {
      $lab = 'RW';
    }
  }
  #print "clone $clone is in lab $lab\n";

  # calculate the window value as blocks of 10 kb
  my $window =  int($chrom_start/10000);

  # if there is no strand information, put one anomaly in for each strand
  if ($chrom_strand ne '+' && $chrom_strand ne '-') {
    &tidy_up_senseless_records($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, '.', $window, $anomaly_score, $explanation, $clone, $clone_start, $clone_end, $lab);
    &put_anomaly_record_in_database($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, '+', $window, $anomaly_score, $explanation, $clone, $clone_start, $clone_end, $lab);
    &put_anomaly_record_in_database($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, '-', $window, $anomaly_score, $explanation, $clone, $clone_start, $clone_end, $lab);

  } else {
    &put_anomaly_record_in_database($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, $chrom_strand, $window, $anomaly_score, $explanation, $clone, $clone_start, $clone_end, $lab);

  }
}

##########################################
# output the record to the database.
# checks if there is a record there already and doesn't change the 'ignore' flag if there is
# already a record in existence.

sub put_anomaly_record_in_database {

  my ($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, $chrom_strand, $window, $anomaly_score, $explanation, $clone, $clone_start, $clone_end, $lab) = @_;


  # ignore this if the score is less than 0.01
  if ($anomaly_score < 0.01) {return;}

  #################################
  # output the data to the GFF file
  #################################

  if ($supplementary) {
    print OUTPUT_GFF "${\$wormbase->chromosome_prefix}$chromosome\tcuration_anomaly\t$anomaly_type\t$chrom_start\t$chrom_end\t$anomaly_score\t$chrom_strand\t.\tEvidence \"$anomaly_id\"\n";
  }

  ####################################
  # output the anomaly to the ace file
  ####################################

  open (OUT, ">> $ace_output") || die "Can't open $ace_output\n";
  print OUT "\n\n";

  # get length of the clone
  my $clone_len = $coords->Superlink_length($clone);

  #Sequence : "AC3"
  #Homol_data "AC3:curation_anomaly" 1 38951

  print OUT "\nSequence : \"$clone\"\n";
  print OUT "Homol_data \"$clone:curation_anomaly\" 1 $clone_len\n\n";

  #Homol_data : "AC3:curation_anomaly"
  #Sequence "AC3"
  #Motif_homol "curation_anomaly" curation_anomaly 15   189   215  1      15

  my $anomaly_len = $chrom_end - $chrom_start + 1;

  print OUT "\nHomol_data : \"$clone:curation_anomaly\"\n";
  print OUT "Sequence \"$clone\"\n";
  if ($chrom_strand eq '+') {
    print OUT "Motif_homol \"$anomaly_type:$anomaly_id:curation_anomaly\" \"curation_anomaly\" $anomaly_score $clone_start $clone_end 1 $anomaly_len\n\n";
  } else {
    print OUT "Motif_homol \"$anomaly_type:$anomaly_id:curation_anomaly\" \"curation_anomaly\" $anomaly_score $clone_start $clone_end $anomaly_len 1\n\n";
  }
      
  #Motif : "signal_peptide_motif"
  #Homol_homol "AC3:signal_peptide"

  print OUT "\nMotif : \"$anomaly_type:$anomaly_id:curation_anomaly\"\n";
  print OUT "Homol_homol \"$clone:curation_anomaly\"\n\n";
    
  close(OUT);

  ########################################
  # write the anomaly data to the database
  ########################################

  # need to do a check for a very similar previous record that may
  # need to be overwritten, preserving the status.

  my $nearest_db_key_id = -1;
  my $nearest = 21;	# this default size of distance will cause a new record to be inserted if it is not changed in the test below to be <= 20

  if (! $nodb) {
    if (! $go_faster_by_ignoring_db_checks) {
      # see if there is something there already
      my $db_query = $mysql->prepare ( qq{ SELECT anomaly_id, chromosome_start, chromosome_end FROM anomaly WHERE type = "$anomaly_type" AND chromosome = "$chromosome" AND sense = "$chrom_strand" AND thing_id = "$anomaly_id" AND window = $window});

      $db_query->execute();
      my $ref_results = $db_query->fetchall_arrayref;

      # find the nearest one to the current data
      #print "\tstart search for $anomaly_id\n";
      foreach my $result_row (@$ref_results) {
	my $dist_start = abs($chrom_start - $result_row->[1]);
	my $dist_end = abs($chrom_end - $result_row->[2]);
	#print "db_id=$result_row->[0] chrom_start=$result_row->[1] chrom_end=$result_row->[2]\n";
	#print "searching distance of start pos = $dist_start end pos = $dist_end, nearest = $nearest\n";
	if ($dist_start + $dist_end < $nearest) {
	  $nearest = $dist_start + $dist_end;
	  $nearest_db_key_id = $result_row->[0];
	  #print "got a new best distance $nearest\n";
	}
      }
    }
  }

  # is the distance in $nearest less than 20 bases, rather than the default size of 21?
  # if it is not zero this is probably a move of the anomaly 
  # as a result of genome sequence changes or
  # changes in the blast database size.
  # so we should update the existing record
  if ($test) {
    print "In test mode - so not updating the mysql database\n";
  } elsif (!$nodb) {
    if ($nearest <= 20) {
      $mysql->do(qq{ UPDATE anomaly SET clone="$clone", clone_start=$clone_start, clone_end=$clone_end, centre="$lab", chromosome_start=$chrom_start, chromosome_end=$chrom_end, thing_score=$anomaly_score, explanation="$explanation"   WHERE anomaly_id = $nearest_db_key_id; });
      # NB we do not write the status record for this anomaly_id

    } else {

      # we want a new record inserted
      # write the data to the database
      $db_key_id++;
      $mysql->do(qq{ insert into anomaly values ($db_key_id, "$anomaly_type", "$clone", $clone_start, $clone_end, "$lab", "$chromosome", $chrom_start, $chrom_end, "$chrom_strand", "$anomaly_id", $anomaly_score, "$explanation", $window, 1, NULL); });
      #print "*** inserting new record\n";
    }
  }


  # and output the line for the St. Louis datafile
  print DAT "INSERT\t$anomaly_type\t$chromosome\t$anomaly_id\t$chrom_start\t$chrom_end\t$chrom_strand\t$anomaly_score\t$explanation\t$clone\t$clone_start\t$clone_end\t$lab\n";


}
##########################################
# there have been a lot of records put in the database that have a sense
# as '.'  we want to change these to be two records with a sense '+
# and a sense '-' (this routine will not be required after the first
# time it has been run and cleaned up the database)

sub tidy_up_senseless_records {

  my ($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, $chrom_strand, $window, $anomaly_score, $explanation, $clone, $clone_start, $clone_end, $lab) = @_;

  my $old_db_key_id;
  my $active;

  if ($test) {
    print "In test mode - so not updating the mysql database\n";
  } elsif (!$nodb) {
    my $db_query = $mysql->prepare ( qq{ SELECT anomaly_id, active FROM anomaly WHERE type = "$anomaly_type" AND chromosome = "$chromosome" AND sense = "." AND chromosome_start = $chrom_start AND chromosome_end = $chrom_end AND thing_id = "$anomaly_id" });

    $db_query->execute();
    my $ref_results = $db_query->fetchall_arrayref;

    # find the nearest one to the current data
    #print "\tstart search for $anomaly_id\n";
    foreach my $result_row (@$ref_results) {
      $old_db_key_id = $result_row->[0];
      $active = $result_row->[1];
    
      # we want two new records inserted
      $db_key_id++;
      #print qq{ INSERT INTO anomaly VALUES ($db_key_id, "$anomaly_type", "$clone", $clone_start, $clone_end, "$lab", "$chromosome", $chrom_start, $chrom_end, "+", "$anomaly_id", $anomaly_score, "$explanation", $window, $active, NULL); \n  };
      
      $mysql->do(qq{ INSERT INTO anomaly VALUES ($db_key_id, "$anomaly_type", "$clone", $clone_start, $clone_end, "$lab", "$chromosome", $chrom_start, $chrom_end, "+", "$anomaly_id", $anomaly_score, "$explanation", $window, $active, NULL); });
      $db_key_id++;
      $mysql->do(qq{ INSERT INTO anomaly VALUES ($db_key_id, "$anomaly_type", "$clone", $clone_start, $clone_end, "$lab", "$chromosome", $chrom_start, $chrom_end, "-", "$anomaly_id", $anomaly_score, "$explanation", $window, $active, NULL); });
      				# 
    
      # and we want to delete the old record with no sense
      #print qq{ DELETE FROM anomaly WHERE anomaly_id = $old_db_key_id; \n  };

      $mysql->do(qq{ DELETE FROM anomaly WHERE anomaly_id = $old_db_key_id; });
    }
  }

}

##########################################
# now delete things that have not been updated in this run that you
# would expect to have been updated like protein-homology-based
# anomalies that might have gone away.  This also deletes anomalies
# that we are no longer putting into the database and which can now be
# removed.
#
# This means that an anomaly based on something like a similarity to a
# protein where the protein is no longer existing will be removed from
# the database

sub delete_anomalies{

  my ($type) = @_;

  # Delete anything that hasn't been marked as to be ignored (still
  # active = 1) that is of the required type
  if ($test) {
    print "In test mode - so not updating the mysql database\n";
  } elsif (!$nodb)  {
    $mysql->do(qq{ DELETE FROM anomaly WHERE type = "$type" AND active = 1 });
  }

  # and output the line for the St. Louis datafile
  print DAT "DELETE\t$type\n";

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

=head2 NAME - find_anomalies.pl

=head1 USAGE

=over 4

=item find_anomalies.pl  [-options]

=back

This script populates the worm_anomaly mysql database with data describing some types of anomalies.

It can be run periodicaly e.g. every build which is especially useful
if there have been corrections made to the genomic sequence

These anomalies can be inspected using the script:

history_maker.pl -anomalies -chromosome X

script_template.pl MANDATORY arguments:

=over 4

=item none

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -database

This specifies an ACeDB database to use to read GFF data etc.

The default is to use autoace.

=over 4

=item -nodb

If this is set then everything is run except for writing to the mysql database.

The default is to write to the mysql database.

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Verbose/Debug mode
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything in the mysql database.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back
                                                                                             

=head1 REQUIREMENTS

=over 4

=item There must be a mysql database server running.

=back

=head1 AUTHOR

=over 4

=item Gary Williams (gw3@sanger.ac.uk)

=back

=cut
