#!/usr/local/bin/perl5.8.0 -w
#
# find_anomalies.pl                           
# 
# by Gary Williams                        
#
# This looks for anomalous thnigs such as protein homologies not
# matching a CDS and stores the results in the mysql database
# 'worm_anomaly'
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2006-06-26 12:31:32 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
#use Sequence_extract;
use Coords_converter;
use DBI;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($database);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "database:s" => \$database,	    # use the specified database instead of autoace
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

#################################
# Set up some useful paths      #
#################################

# Set up top level base directories (these are different if in test mode)
my $basedir         = $wormbase->basedir;     # BASE DIR
my $ace_dir         = $wormbase->autoace;     # AUTOACE DATABASE DIR
my $wormpep_dir     = $wormbase->wormpep;     # CURRENT WORMPEP
my $wormrna_dir     = $wormbase->wormrna;     # CURRENT WORMRNA
my $common_data_dir = $wormbase->common_data; # AUTOACE COMMON_DATA
my $chromosomes_dir = $wormbase->chromosomes; # AUTOACE CHROMSOMES
my $reports_dir     = $wormbase->reports;     # AUTOACE REPORTS
my $gff_dir         = $wormbase->gff;         # AUTOACE GFF
my $gff_splits_dir  = $wormbase->gff_splits;  # AUTOACE GFF SPLIT
my $logs_dir        = $wormbase->logs;        # AUTOACE LOGS

# some database paths
my $geneace   = $wormbase->database('geneace');
my $camace    = $wormbase->database('camace');
my $currentdb = $wormbase->database('current');
my $stlace    = $wormbase->database('stlace');
my $citace    = $wormbase->database('citace');
my $cshace    = $wormbase->database('cshace');
my $brigace   = $wormbase->database('brigace');

# other paths
my $ftp_upload_dir  = $wormbase->ftp_upload;  # "/nfs/ftp_uploads/wormbase"
my $tace            = $wormbase->tace;        # TACE PATH
my $giface          = $wormbase->giface;      # GIFACE PATH



##########################
# MAIN BODY OF SCRIPT
##########################

# mysql database parameters
my $dbsn = "DBI:mysql:database=worm_anomaly;host=ecs1f";
my $dbuser = "wormadmin";
my $dbpass = "worms";

my $mysql = DBI -> connect($dbsn, $dbuser, $dbpass, {RaiseError => 1})
      || die "cannot connect to database, $DBI::errstr";

# get the last used anomaly_id key value
my $array_ref = $mysql->selectcol_arrayref("select max(anomaly_id) from anomaly;");
my $db_key_id = $array_ref->[0];
my $go_faster_by_ignoring_db_checks = 0;
if (defined $db_key_id) {
  print "db_key_id=$db_key_id\n";
} else {
  # reset the daatbase key value
  $db_key_id = 0; 
  print "db_key_id has been reset to 0\n";
  $go_faster_by_ignoring_db_checks = 1;	# don't need to check records in the database because there are none
}

$database = $wormbase->autoace if (!defined $database || $database eq "");

my $coords = Coords_converter->invoke($database, 0, $wormbase);
my %clonelab = $wormbase->FetchData("clone2centre", undef, "$database/COMMON_DATA");

# open an ACE connection to parse details for mapping to genome
print "Connecting to Ace\n";
my $ace = Ace->connect (-path => $database,
                       -program => $tace) || die "cannot connect to database at $database\n";




my @chromosomes = qw( I II III IV V X );

foreach my $chromosome (@chromosomes) {

  $log->write_to("Processing chromosome $chromosome\n");
                                                                                                   
  # get the exon data
  print "reading exons\n";
  my @exons = &get_coding_exons($database, $chromosome);

  print "reading pseudogenes\n";
  my @pseudogenes = &get_pseudogenes($database, $chromosome);

  print "reading transposons\n";
  my @transposons = &get_transposons($database, $chromosome);

  print "reading transposon coding exons\n";
  my @transposon_exons = &get_transposon_exons($database, $chromosome);

  print "reading coding transcripts\n";
  my @coding_transcripts = &get_coding_transcripts($database, $chromosome);

  print "reading noncoding transcripts exons\n";
  my @noncoding_transcripts = &get_noncoding_transcripts($database, $chromosome);

  print "reading rRNA transcripts\n";
  my @rRNA = &get_rRNA($database, $chromosome);

  print "reading SL1 TSL sites\n";
  my @TSL_SL1 = &get_TSL_SL1($database, $chromosome);

  print "reading SL1 TSL sites\n";
  my @TSL_SL2 = &get_TSL_SL2($database, $chromosome);

  print "reading SAGE transcripts\n";
  my @SAGE_transcripts = &get_SAGE_transcripts($database, $chromosome);

  # get the wublastx protein homology data
  print "reading protein homologies\n";
  my @homologies = &get_homologies($database, $chromosome);

  # get the twinscan data
  print "reading twinscan exons\n";
  my @twinscan = &get_twinscan($database, $chromosome);

  # get the genefinder data
  print "reading genefinder exons\n";
  my @genefinder = &get_genefinder($database, $chromosome);

  # get the EST BLAT homology data
  print "reading est data\n";
  my @est = &get_est($database, $chromosome);




  # get the homologies showing no match to any exon or pseudogene
  # and those with a match to a coding exon
  print "finding protein homologies not overlapping exons\n";
  my $matched_protein_aref = &get_protein_differences(\@exons, \@pseudogenes, \@homologies, \@transposons, \@transposon_exons, \@noncoding_transcripts, \@rRNA, $chromosome);

  print "finding frameshifts\n";
  &get_frameshifts($matched_protein_aref, $chromosome);

  print "finding split/merged genes based on protein homology\n";
  &get_protein_split_merged($matched_protein_aref, $chromosome);

  print "read confirmed_introns and other GFF files\n";
  &get_GFF_files($chromosome);

  # this looks at the EST TSL sites and finds those not attached to genes
  print "finding TSL sites not attached to genes\n";
  &get_unattached_TSL(\@TSL_SL1, \@TSL_SL2, $chromosome);

  # this finds TEC-RED TSL sites more than 100 bases upstream that are not mentioned in the remarks or evidence
  print "finding isolated TSL sites\n";
  &get_isolated_TSL(\@TSL_SL1, \@TSL_SL2, \@coding_transcripts, \@pseudogenes, \@transposons, \@transposon_exons, \@noncoding_transcripts, \@rRNA, $chromosome);

  # get SAGE tags that don't match a gene with score based on frequency
  print "finding non-overlapping SAGE_transcripts\n";
  &get_unmatched_SAGE(\@exons, \@pseudogenes, \@SAGE_transcripts, \@transposons, \@transposon_exons, \@noncoding_transcripts, \@rRNA, $chromosome);

  print "finding twinscan exons not overlapping curated regions\n";
  &get_unmatched_twinscan_exons(\@twinscan, \@exons, \@pseudogenes, \@transposons, \@transposon_exons, \@noncoding_transcripts, \@rRNA, $chromosome);

  print "finding genefinder exons not overlapping curated regions\n";
  &get_unmatched_genefinder_exons(\@genefinder, \@exons, \@pseudogenes, \@transposons, \@transposon_exons, \@noncoding_transcripts, \@rRNA, $chromosome);

  # ESTs not matching exons
  # and those with a match to a coding exon
  print "finding EST homologies not overlapping exons\n";
  my $matched_EST_aref = &get_EST_differences(\@exons, \@pseudogenes, \@est, \@transposons, \@transposon_exons, \@noncoding_transcripts, \@rRNA, $chromosome);

  print "finding split/merged genes based on EST homology\n";
  &get_EST_split_merged($matched_EST_aref, $chromosome);

}



# disconnect from the mysql database
$mysql->disconnect || die "error disconnecting from database", $DBI::errstr;

# close the ACE connection
$ace->close;

$log->mail();

print "Finished.\n" if ($verbose);
exit(0);


##############################################################
#
# Subroutines
#
##############################################################

# get the coding exons

sub get_coding_exons {

  my ($database, $chromosome) = @_;

  my %GFF_data = 
   (
     directory			=> "$database/GFF_SPLITS/",
     file			=> "CHROMOSOME_${chromosome}_Coding_transcript.gff",
     gff_source			=> "Coding_transcript",
     gff_type			=> "exon",
     anomaly_type		=> "",
     ID_after			=> "Transcript\\s+",
     action                     => ["return_result"],
   );

  return &read_GFF_file(\%GFF_data);

}

##########################################
# read the pseudogenes
#  my @coding_transcripts = &get_coding_transcripts($database, $chromosome);

sub get_pseudogenes {
  my ($database, $chromosome) = @_;

  my %GFF_data = 
   (
     directory			=> "$database/GFF_SPLITS/",
     file			=> "CHROMOSOME_${chromosome}_Pseudogene.gff",
     gff_source			=> "Pseudogene",
     gff_type			=> "Pseudogene",
     anomaly_type		=> "",
     ID_after			=> "Pseudogene\\s+",
     action                     => ["return_result"],
   );

  return &read_GFF_file(\%GFF_data);

}


##########################################

# get transposons
# homologies to these should be ignored, and we don't need to check the frame


sub get_transposons {
  my ($database, $chromosome) = @_;

  my %GFF_data = 
   (
     directory			=> "$database/GFF_SPLITS/",
     file			=> "CHROMOSOME_${chromosome}_Transposon.gff",
     gff_source			=> "Transposon",
     gff_type			=> "transposable_element",
     anomaly_type		=> "",
     ID_after			=> "Transposon\\s+",
     action                     => ["return_result"],
   );

  return &read_GFF_file(\%GFF_data);

}

##########################################

# get transposon_exons
# homologies to these should be ignored, and we don't need to check the frame
# these are the coding exons of the Transposon_CDS objects

sub get_transposon_exons {
  my ($database, $chromosome) = @_;

  my %GFF_data = 
   (
     directory			=> "$database/GFF_SPLITS/",
     file			=> "CHROMOSOME_${chromosome}_Transposon_CDS.gff",
     gff_source			=> "Transposon_CDS",
     gff_type			=> "coding_exon",
     anomaly_type		=> "",
     ID_after			=> "CDS\\s+",
     action                     => ["return_result"],
   );

  return &read_GFF_file(\%GFF_data);

}

##########################################
# get the coding transcript data
#  my @coding_transcripts = &get_coding_transcripts($database, $chromosome);

sub get_coding_transcripts {
  my ($database, $chromosome) = @_;

  my %GFF_data = 
   (
     directory			=> "$database/GFF_SPLITS/",
     file			=> "CHROMOSOME_${chromosome}_Coding_transcript.gff",
     gff_source			=> "Coding_transcript",
     gff_type			=> "protein_coding_primary_transcript",
     anomaly_type		=> "",
     ID_after			=> "Transcript\\s+",
     action                     => ["return_result"],
   );

  return &read_GFF_file(\%GFF_data);

}


##########################################
# get the noncoding transcript data
#  my @transcripts = get_noncoding_transcripts($database, $chromosome);

sub get_noncoding_transcripts {
  my ($database, $chromosome) = @_;

  my %GFF_data = 
   (
     directory			=> "$database/GFF_SPLITS/",
     file			=> "CHROMOSOME_${chromosome}_Non_coding_transcript.gff",
     gff_source			=> "Non_coding_transcript",
     gff_type			=> "exon",
     anomaly_type		=> "",
     ID_after			=> "Transcript\\s+",
     action                     => ["return_result"],
   );

  return &read_GFF_file(\%GFF_data);

}

##########################################
# get the ncRNAs
#  my @rRNA = &get_rRNA($database, $chromosome);

sub get_rRNA {
  my ($database, $chromosome) = @_;

  my %GFF_data = 
   (
     directory			=> "$database/GFF_SPLITS/",
     file			=> "CHROMOSOME_${chromosome}_Non_coding_transcript.gff",
     gff_source			=> "rRNA",
     gff_type			=> "rRNA_primary_transcript",
     anomaly_type		=> "",
     ID_after			=> "Transcript\\s+",
     action                     => ["return_result"],
   );

  return &read_GFF_file(\%GFF_data);

}

##########################################
# get the wublastx protein homologies
#  my @wublastx = get_homologies($database, $chromosome);

sub get_homologies {
  my ($database, $chromosome) = @_;


  my %GFF_data = 
   (
     directory			=> "$database/CHROMOSOMES/", # NB we are reading the full gff file, not the split ones here
     file			=> "CHROMOSOME_${chromosome}.gff",
     gff_source			=> "wublastx",
     gff_type			=> "protein_match",
     homology			=> "1",	# this is a GFF with homology data that we need to store
     anomaly_type		=> "",
     ID_after			=> "Target\\s+\"Protein:",
     action                     => ["return_result"],
   );

  return &read_GFF_file(\%GFF_data);
  
}

##########################################
# get the SL1 TSL sites
#  my @TSL_SL1 = &get_TSL_SL1($database, $chromosome);


sub get_TSL_SL1 {
  my ($database, $chromosome) = @_;

  my %GFF_data = 
   (
     directory			=> "$database/GFF_SPLITS/", # NB we are reading the full gff file, not the split ones here
     file			=> "CHROMOSOME_${chromosome}_SL1.gff",
     gff_source			=> "SL1",
     gff_type			=> "SL1_acceptor_site",
     anomaly_type		=> "",
     ID_after			=> "Feature\\s+",
     action                     => ["return_result"],
   );

  return &read_GFF_file(\%GFF_data);
  
}

##########################################
# get the SL2 TSL sites
#  my @TSL_SL2 = &get_TSL_SL2($database, $chromosome);


sub get_TSL_SL2 {
  my ($database, $chromosome) = @_;

  my %GFF_data = 
   (
     directory			=> "$database/GFF_SPLITS/", # NB we are reading the full gff file, not the split ones here
     file			=> "CHROMOSOME_${chromosome}_SL2.gff",
     gff_source			=> "SL2",
     gff_type			=> "SL2_acceptor_site",
     anomaly_type		=> "",
     ID_after			=> "Feature\\s+",
     action                     => ["return_result"],
   );

  return &read_GFF_file(\%GFF_data);
  
}

##########################################
# read the SAGE transcripts
#  my @SAGE_transcripts = &get_SAGE_transcripts($database, $chromosome);

sub get_SAGE_transcripts {
  my ($database, $chromosome) = @_;

  my %GFF_data = 
   (
     directory			=> "$database/GFF_SPLITS/",
     file			=> "CHROMOSOME_${chromosome}_SAGE_transcript.gff",
     gff_source			=> "SAGE_transcript",
     gff_type			=> "transcript",
     anomaly_type		=> "",
     ID_after			=> "SAGE_transcript\\s+",
     action                     => ["return_result"],
   );

  return &read_GFF_file(\%GFF_data);

}

##########################################
# get the twinscan exons
#  my @twinscan = get_twinscan($database, $chromosome);

sub get_twinscan {
  my ($database, $chromosome) = @_;


  my %GFF_data = 
   (
     directory			=> "$database/CHROMOSOMES/", # NB we are reading the full gff file, not the split ones here
     file			=> "CHROMOSOME_${chromosome}.gff",
     gff_source			=> "twinscan",
     gff_type			=> "coding_exon",
     anomaly_type		=> "",
     ID_after			=> "CDS\\s+",
     action                     => ["return_result"],
   );

  return &read_GFF_file(\%GFF_data);

}


##########################################
# get the genefinder exons
#  my @genefinder = get_genefinder($database, $chromosome);

sub get_genefinder {
  my ($database, $chromosome) = @_;


  my %GFF_data = 
   (
     directory			=> "$database/CHROMOSOMES/", # NB we are reading the full gff file, not the split ones here
     file			=> "CHROMOSOME_${chromosome}.gff",
     gff_source			=> "Genefinder",
     gff_type			=> "coding_exon",
     anomaly_type		=> "",
     ID_after			=> "CDS\\s+",
     action                     => ["return_result"],
   );

  return &read_GFF_file(\%GFF_data);

}

##########################################
# get the EST BLAT homologies - (only BLAT_EST_BEST, the BLAT_EST_OTHER doesn't have a sense associated with it)
#  my @est = &get_est($database, $chromosome);

sub get_est {
  my ($database, $chromosome) = @_;


  my %GFF_data = 
   (
     directory			=> "$database/CHROMOSOMES/", # NB we are reading the full gff file, not the split ones here
     file			=> "CHROMOSOME_${chromosome}.gff",
     gff_source			=> "BLAT_EST_BEST",
     gff_type			=> "EST_match",
     homology			=> "1",	# this is a GFF with homology data that we need to store
     anomaly_type		=> "",
     ID_after			=> "Target\\s+\"Sequence:",
     action                     => ["return_result"],
   );

  # The 3' reads of most ofthe ESTs are really in the reverse sense to
  # the gene they match but ACeDB has a tag to display them in the
  # opposite sense as this looks much better.
  #
  # The GFF file holds the original sense though, so we do an explicit
  # flip of the sense of any 3' reads so that when we come to display
  # the ESTs in the curation tool, we get the dispaly coming up the
  # right way round.
  #
  # This will not be infallible, because there are some ESTs where the
  # orientation has been produced incorrectly and changed by manually
  # setting some tag or other in the EST object and there are some
  # ESTs where your can't see if they are 3' from their ID name, but
  # it should work more often than not.

  my @ests = &read_GFF_file(\%GFF_data);
  my @result;

  foreach my $est (@ests) {
    #my $protein_id = $est->[0];
    #my $chrom_strand = $est->[3];
    if ($est->[0] =~ /\.3$/) {
      $est->[3] = (($est->[3] eq '+') ? '-' : '+'); # flip the sense
    }
    push @result, $est;
  }


  return @result;

}

##########################################
# get the homologies with no matching exons or pseudogenes or transposons
# and those which do match exons or transposons (but not pseudogenes)

#  my @matched_homologies = get_differences(\@transcripts, \@pseudogenes, \@protein_coverage);

sub get_protein_differences {
  my ($exons_aref, $pseudogenes_aref, $homologies_aref, $transposons_aref, $transposon_exons_aref, $noncoding_transcripts_aref, $rRNA_aref, $chromosome) = @_;

  my @homologies = @{$homologies_aref};

  my @exons = @{$exons_aref};

  my @not_matched = ();		# the resulting list of hashes of homologies with no matching exons/transposons/pseudogenes
  my @matched = ();		# the resulting list of hashes of homologies which match a coding exon

  my $SMALL_OVERLAP = -20;	# amount of small overlap to ignore

  my %exons_match = &init_match(near => $SMALL_OVERLAP, same_sense => 1); # allow the protein to have up to 20 bases overlap without it being a match
  my %pseud_match = &init_match();
  my %trans_match = &init_match();
  my %trane_match = &init_match();
  my %nonco_match = &init_match();
  my %rrna_match = &init_match();

  foreach my $homology (@homologies) { # $protein_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_match = 0;	        # not yet seen a match to anything
    my $got_a_match_to_coding_exon = 0;	# not yet seen a match to a coding exon
    my $matching_exon = "";	        # the name of the exon that matches;

    if (&match($homology, $exons_aref, \%exons_match)) {
      $got_a_match = 1;
      $got_a_match_to_coding_exon = 1;
      $matching_exon = $exons_match{'matching_ids'};
   }

    if (&match($homology, $pseudogenes_aref, \%pseud_match)) {
      $got_a_match = 1;
    }

    if (&match($homology, $transposons_aref, \%trans_match)) {
      $got_a_match = 1;
      $got_a_match_to_coding_exon = 1;
      $matching_exon = $exons_match{'matching_ids'};
    }

    if (&match($homology, $transposon_exons_aref, \%trane_match)) {
      $got_a_match = 1;
    }

    if (&match($homology, $noncoding_transcripts_aref, \%nonco_match)) {
      $got_a_match = 1;
    }

    if (&match($homology, $rRNA_aref, \%rrna_match)) {
      $got_a_match = 1;
    }

    # output unmatched homologies to the database
    if (! $got_a_match) {
      my $protein_id = $homology->[0];
      my $chrom_start = $homology->[1];
      my $chrom_end = $homology->[2];
      my $chrom_strand = $homology->[3];
      my $protein_score = $homology->[6];
      # make the anomaly score based on the protein alignment score normalised between 0 and 1
      # the BLAST scores seem to be between 0 and 1000
      my $anomaly_score = $protein_score/1000;
      if ($anomaly_score > 1) {$anomaly_score = 1;}
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
# my $matched_EST_aref = &get_EST_differences(\@exons, \@pseudogenes, \@est, \@transposons, \@transposon_exons, \@noncoding_transcripts, $chromosome);

sub get_EST_differences {
  my ($exons_aref, $pseudogenes_aref, $est_aref, $transposons_aref, $transposon_exons_aref, $noncoding_transcripts_aref, $rRNA_aref, $chromosome) = @_;

  my @est = @{$est_aref};

  my @exons = @{$exons_aref};

  my @not_matched = ();		# the resulting list of hashes of est with no matching exons/transposons/pseudogenes
  my @matched = ();		# the resulting list of hashes of est which match a coding exon


  my %exons_match = &init_match(same_sense => 1); # allow the protein to have up to 20 bases overlap without it being a match
  my %pseud_match = &init_match();
  my %trans_match = &init_match();
  my %trane_match = &init_match();
  my %nonco_match = &init_match();
  my %rrna_match = &init_match();

  foreach my $homology (@est) { # $protein_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_match = 0;	        # not yet seen a match to anything
    my $got_a_match_to_coding_exon = 0;	# not yet seen a match to a coding exon
    my $matching_exon = "";	        # the name of the exon that matches;

    if (&match($homology, $exons_aref, \%exons_match)) {
      $got_a_match = 1;
      $got_a_match_to_coding_exon = 1;
      $matching_exon = $exons_match{'matching_ids'};
   }

    if (&match($homology, $pseudogenes_aref, \%pseud_match)) {
      $got_a_match = 1;
    }

    if (&match($homology, $transposons_aref, \%trans_match)) {
      $got_a_match = 1;
      $got_a_match_to_coding_exon = 1;
      $matching_exon = $exons_match{'matching_ids'};
    }

    if (&match($homology, $transposon_exons_aref, \%trane_match)) {
      $got_a_match = 1;
    }

    if (&match($homology, $noncoding_transcripts_aref, \%nonco_match)) {
      $got_a_match = 1;
    }

    if (&match($homology, $rRNA_aref, \%rrna_match)) {
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

	# make the anomaly score based on the protein alignment score normalised between 0 and 1
	# the BLAST scores seem to be between 0 and 1000
	# use the average of the two alignment's scores to discourage reporting cases where one is very poor
	my $anomaly_score = ($protein_score+$prev_score)/2000;
	if ($anomaly_score > 1) {$anomaly_score = 1;}
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
    my $matching_exon = $homology->[7];

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

      # see if the protein HSP order is unchanged - a possible case for merging if continued over nito the next gene
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
      # distance between the alignments is less than 10 kb
      if ($got_a_new_exon && $got_a_continuation_of_the_HSPs && $chrom_start - $prev_chrom_end < 10000) {
	# output to database
	# make the anomaly score based on the protein alignment score normalised between 0 and 1
	# the BLAST scores seem to be between 0 and 1000
	# use the average of the two alignment's scores to discourage reporting cases where one is very poor
	my $anomaly_score = ($protein_score+$prev_score)/2000;
	if ($anomaly_score > 1) {$anomaly_score = 1;}
	if ($anomaly_score < 0) {$anomaly_score = 0;}
	#print "MERGE genes ANOMALY $prev_exon and $matching_exon\t$protein_id, $prev_chrom_start, $chrom_end, $chrom_strand, $anomaly_score, 'merge $prev_exon and $matching_exon'\n";
	&output_to_database("MERGE_GENES_BY_PROTEIN", $chromosome, $protein_id, $prev_chrom_start, $chrom_end, $chrom_strand, $anomaly_score, "merge $prev_exon and $matching_exon");


      }

      if ($got_a_big_decrease_in_HSP_start && ! $got_a_new_exon) {
	# output to database
	# make the anomaly score based on the protein alignment score normalised between 0 and 1
	# the BLAST scores seem to be between 0 and 1000
	# use the average of the two alignment's scores to discourage reporting cases where one is very poor
	my $anomaly_score = ($protein_score+$prev_score)/2000;
	if ($anomaly_score > 1) {$anomaly_score = 1;}
	if ($anomaly_score < 0) {$anomaly_score = 0;}
	#print "SPLIT gene ANOMALY $matching_exon\t$protein_id, $prev_chrom_start, $chrom_end, $chrom_strand, $anomaly_score, 'split $matching_exon'\n";
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
#
#  &get_EST_split_merged($matched_EST_aref, $chromosome);

# find groups of homology that indicate that the genes they match should be split or merged
# look for non-overlapping HSPs that jump back down the position in the protein that is aligned as you move along the chromosome -> split
# look for matches to exons with HSPs that do not jump back down the position in the protein that is aligned as you move along the chromosome -> merge

sub get_EST_split_merged {
  my ($matched_EST_aref, $chromosome) = @_;

  my @matched = @{$matched_EST_aref};

  my $prev_EST_id = "";	# use to collect alignments for a EST's homology when looking for frameshifts
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


  # sort the homologies grouped by EST ID and then chromosomal position
  my @homologies = sort {$a->[0] cmp $b->[0] or $a->[1] <=> $b->[1]} @matched;

  foreach my $homology (@homologies) { # $EST_id, $chrom_start, $chrom_end, $chrom_strand, $hit_start, $hit_end, $EST_score, $matching_exon

    my $EST_id = $homology->[0];
    my $chrom_start = $homology->[1];
    my $chrom_end = $homology->[2];
    my $chrom_strand = $homology->[3];
    my $hit_start = $homology->[4];
    my $hit_end = $homology->[5];
    my $EST_score = $homology->[6];
    my $matching_exon = $homology->[7];

    #  check if these are isoforms of the same gene and if so then treat them as the same gene
    if ($matching_exon =~ /(\S+\.\d+)[a-z]/) {
      $matching_exon = $1;
    }
    # and change the strange transcript names like H25P06.1.1, H25P06.1.2 to H25P06.1
    if ($matching_exon =~ /(\S+\.\d+)\.\d+/) {
      $matching_exon = $1;
    }

    #print "Matched: $EST_id, $chrom_start, $chrom_end, $chrom_strand, $hit_start, $hit_end, $EST_score, $matching_exon\n";

    # look for homologies which cover two genes or genes which cover a repeated homology
    if ($EST_id eq $prev_EST_id && $chrom_strand eq $prev_chrom_strand) {

      # see if exon has changed
      if ($prev_exon ne "" && $prev_exon ne $matching_exon) {
	$got_a_new_exon = 1;
	#print "Merge? prev $prev_exon this $matching_exon\n";
      } else {
	$got_a_new_exon = 0;
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
      if ($EST_id =~ /\.3$/) {	# flip the logic of these around
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
      if ($got_a_new_exon && $got_a_continuation_of_the_HSPs && $chrom_start - $prev_chrom_end < 10000) {
	# output to database
	# make the anomaly score based on the EST alignment score normalised between 0 and 1
	# the BLAT scores are between 0 and 100
	# use the average of the two alignment's scores to discourage reporting cases where one is very poor
	my $anomaly_score = ($EST_score+$prev_score)/200;
	if ($anomaly_score > 1) {$anomaly_score = 1;}
	if ($anomaly_score < 0) {$anomaly_score = 0;}
	#print "MERGE genes ANOMALY $prev_exon and $matching_exon\t$EST_id, $prev_chrom_start, $chrom_end, $chrom_strand, $anomaly_score, 'merge $prev_exon and $matching_exon'\n";
	&output_to_database("MERGE_GENES_BY_EST", $chromosome, $EST_id, $prev_chrom_start, $chrom_end, $chrom_strand, $anomaly_score, "merge $prev_exon and $matching_exon");


      }

      if ($got_a_big_decrease_in_HSP_start && ! $got_a_new_exon) {
	# output to database
	# make the anomaly score based on the EST alignment score normalised between 0 and 1
	# the BLAT scores are between 0 and 100
	# use the average of the two alignment's scores to discourage reporting cases where one is very poor
	my $anomaly_score = ($EST_score+$prev_score)/200;
	if ($anomaly_score > 1) {$anomaly_score = 1;}
	if ($anomaly_score < 0) {$anomaly_score = 0;}
	#print "SPLIT gene ANOMALY $matching_exon\t$EST_id, $prev_chrom_start, $chrom_end, $chrom_strand, $anomaly_score, 'split $matching_exon'\n";
	&output_to_database("SPLIT_GENES_BY_EST", $chromosome, $EST_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, "split $matching_exon");
      }
    }

    $prev_exon = $matching_exon;
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
# find TSL sites that are defined by an EST but have no association with a transcript in their object
#  &get_unattached_TSL(\@TSL_SL1, \@TSL_SL2)

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
      if ($sequence =~ /\.5$/ || $sequence =~ /\.3$/) {
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
#  &get_isolated_TSL(\@TSL_SL1, \@TSL_SL2, \@coding_transcript, \@transposons, \@transposon_exons, \@noncoding_transcripts, $chromosome);

sub get_isolated_TSL {

  my ($SL1_aref, $SL2_aref, $transcripts_aref, $pseudogenes_aref, $transposons_aref, $transposon_exons_aref, $noncoding_transcripts_aref, $rRNA_aref, $chromosome) = @_;

  my @SL1 = @$SL1_aref;
  my @SL2 = @$SL2_aref;
  my @SL = sort {$a->[1] <=> $b->[1]} (@SL1, @SL2); # merge and sort the two TSL lists

  my %transcripts_match = &init_match(near => 20, same_sense => 1); # allow the TSL to be within 20 bases of the transcript to give a match
  my %pseud_match = &init_match(same_sense => 1);
  my %trans_match = &init_match(same_sense => 1);
  my %trane_match = &init_match(same_sense => 1);
  my %nonco_match = &init_match(same_sense => 1);
  my %rrna_match = &init_match(same_sense => 1);

  foreach my $tsl (@SL) { # $TSL_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_match = 0;
  
    if (&match($tsl, $transcripts_aref, \%transcripts_match)) {
      $got_a_match = 1;
    }

    if (&match($tsl, $pseudogenes_aref, \%pseud_match)) {
      $got_a_match = 1;
    }

    if (&match($tsl, $transposons_aref, \%trans_match)) {
      $got_a_match = 1;
    }

    if (&match($tsl, $transposon_exons_aref, \%trane_match)) {
      $got_a_match = 1;
    }

    if (&match($tsl, $noncoding_transcripts_aref, \%nonco_match)) {
      $got_a_match = 1;
    }

    if (&match($tsl, $rRNA_aref, \%rrna_match)) {
      $got_a_match = 1;
    }

    # output unmatched TSL sites to the database
    if (! $got_a_match) {
      my $TSL_id = $tsl->[0];
      my $chrom_start = $tsl->[1];
      my $chrom_end = $tsl->[2];
      my $chrom_strand = $tsl->[3];

      # make the anomaly score 1
      my $anomaly_score = 1.0;
      #print "TSL NOT got a match ANOMALY: $TSL_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score\n";
      &output_to_database("UNMATCHED_TSL", $chromosome, $TSL_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, '');
    }
  }

}

##########################################
# get twinscan exons that do not match a coding transcript or pseudogene
# &get_unmatched_twinscan_exons(\@twinscan, \@exons, \@pseudogenes, \@transposons, \@transposon_exons, \@noncoding_transcripts, $chromosome)


sub get_unmatched_twinscan_exons {

  my ($twinscan_aref, $exons_aref, $pseudogenes_aref, $transposons_aref, $transposon_exons_aref, $noncoding_transcripts_aref, $rRNA_aref, $chromosome) = @_;

  my %exons_match = &init_match(same_sense => 1);
  my %pseud_match = &init_match(same_sense => 1);
  my %trans_match = &init_match(same_sense => 1);
  my %trane_match = &init_match(same_sense => 1);
  my %nonco_match = &init_match(same_sense => 1);
  my %rrna_match = &init_match(same_sense => 1);

  foreach my $twinscan (@{$twinscan_aref}) { # $twinscan_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_match = 0;
  
    if (&match($twinscan, $exons_aref, \%exons_match)) {
      $got_a_match = 1;
    }

    if (&match($twinscan, $pseudogenes_aref, \%pseud_match)) {
      $got_a_match = 1;
    }

    if (&match($twinscan, $transposons_aref, \%trans_match)) {
      $got_a_match = 1;
    }

    if (&match($twinscan, $transposon_exons_aref, \%trane_match)) {
      $got_a_match = 1;
    }

    if (&match($twinscan, $noncoding_transcripts_aref, \%nonco_match)) {
      $got_a_match = 1;
    }

    if (&match($twinscan, $rRNA_aref, \%rrna_match)) {
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
# &get_unmatched_genefinder_exons(\@genefinder, \@exons, \@pseudogenes, \@transposons, \@transposon_exons, \@noncoding_transcripts, $chromosome)


sub get_unmatched_genefinder_exons {

  my ($genefinder_aref, $exons_aref, $pseudogenes_aref, $transposons_aref, $transposon_exons_aref, $noncoding_transcripts_aref, $rRNA_aref, $chromosome) = @_;

  my %exons_match = &init_match(same_sense => 1);
  my %pseud_match = &init_match(same_sense => 1);
  my %trans_match = &init_match(same_sense => 1);
  my %trane_match = &init_match(same_sense => 1);
  my %nonco_match = &init_match(same_sense => 1);
  my %rrna_match = &init_match(same_sense => 1);

  foreach my $genefinder (@{$genefinder_aref}) { # $genefinder_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_match = 0;
  
    if (&match($genefinder, $exons_aref, \%exons_match)) {
      $got_a_match = 1;
    }

    if (&match($genefinder, $pseudogenes_aref, \%pseud_match)) {
      $got_a_match = 1;
    }

    if (&match($genefinder, $transposons_aref, \%trans_match)) {
      $got_a_match = 1;
    }

    if (&match($genefinder, $transposon_exons_aref, \%trane_match)) {
      $got_a_match = 1;
    }

    if (&match($genefinder, $noncoding_transcripts_aref, \%nonco_match)) {
      $got_a_match = 1;
    }

    if (&match($genefinder, $rRNA_aref, \%rrna_match)) {
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
# resets the match state before a search
# my %state = &init_match(param => arg, ...);
# Args:
#      near => distance - set distance within which a near match is considered to be a match
#      near_5 => distance - ditto for 5' end of our objects
#      near_3 => distance - ditto for 3' end of our objects
#      same_sense => boolean - say it must have both in the same sense for a match

sub init_match {
  my (@params) = @_;

  my %state = (last_used  => 0, # reset last other-list's line 
	       last_used_forward => 0, # reset last other-list's line for the forward sense
	       last_used_reverse => 0, # reset last other-list's line for the reverse sense
	       near_5     => 0, # assume we don't allow near 5' matches to count as a match
	       near_3     => 0,	# ditto for 3'
	       same_sense => 0,	# assume we want to allow a match to an object in either sense, set to true if want the same sense only
	       );

  while (@params) {
    my $param = shift @params;
    if ($param eq "near") {
      my $near = shift @params;
      $state{'near_5'} = $near;
      $state{'near_3'} = $near;
    } elsif ($param eq "near_5") {
      my $near_5 = shift @params;
      $state{'near_5'} = $near_5;
    } elsif ($param eq "near_3") {
      my $near_3 = shift @params;
      $state{'near_3'} = $near_3;
    } elsif ($param eq "same_sense") {
      my $same_sense = shift @params;
      $state{'same_sense'} = $same_sense;

    }

  }
      
  return %state;

}

##########################################
# routine to do generalised matching of the region in this_line versus all of the regions in @other_list
# $result = &match($this_line, \@other_list, \%state);
# 
# Args:
#      $this_line - ref to array holding ID, chrom_start, chrom_end, chrom_strand and possibly other values on the end
#      $other_list - ref to array of arrays sorted by start position - each of which holds ID, chrom_start, chrom_end, chrom_strand and possibly other values on the end
#      $state - ref to hash holding state of the search and controlling parameters for the search and ID of matching results
#             last_used - must be initialised to 0
#             near_5 - the amount of bases for nearby other regions to count as an overlap at the 5' end (negative if want to ignore small overlaps)
#             near_3 - the amount of bases for nearby other regions to count as an overlap allowed at the 3' end (negative if want to ignore small overlaps)
#             matching_ids - returned space-delimited string of other IDs that match
# 
# Returns: 1 if match found, 0 if no match found
#          the status hash holds the IDs of the exons matched in 'matching_ids'

sub match {

  my ($this_line, $other_list, $state) = @_;
  my $got_a_match = 0;		# result of the match - no match found yet

  my $id = $this_line->[0];
  my $chrom_start = $this_line->[1];
  my $chrom_end = $this_line->[2];
  my $chrom_strand = $this_line->[3];
  #print "THIS LINE: $id $chrom_start $chrom_end $chrom_strand\n";

  my $near_3;
  my $near_5;
  if ($chrom_strand eq '+') {
    $near_3 = $state->{'near_3'};
    $near_5 = $state->{'near_5'};
  } else {			# swap the values around in the reverse sense
    $near_3 = $state->{'near_5'};
    $near_5 = $state->{'near_3'};
  }

  $state->{'matching_ids'} = ""; # no matching IDs found yet

  # if searching for same-sense matches, then set last_used to be the minimum of last_used_forward and last_used_reverse
  if ($state->{'same_sense'}) {
    $state->{'last_used'} = $state->{'last_used_forward'};
    if ($state->{'last_used_reverse'} < $state->{'last_used_forward'}) {$state->{'last_used'} = $state->{'last_used_reverse'};}
  }

  #print "last_used $state->{'last_used'}\n";

  for (my $other_start = 0, my $i = $state->{'last_used'}; $other_start < $chrom_end && defined $other_list->[$i]; $i++) {
    my $other = $other_list->[$i];
    # test if there is an overlap (allowing possible nearby matches) and optionally test if the senses are the same
    #print "OTHER LINE: $other->[0] $other->[1] $other->[2] $other->[3]\n";
    if ($other->[1] <= $chrom_end + $near_3 && $other->[2] >= $chrom_start - $near_5 && ($state->{'same_sense'}?($chrom_strand eq $other->[3]):1)) {

      # see if we are testing for them to be in the same sense
      if ($state->{'same_sense'}) {
	if ($chrom_strand eq '+') {
	  # remember where we got up to in this sense
	  $state->{'last_used_forward'} = $i;
	} else {
	  $state->{'last_used_reverse'} = $i;
	}
      }

      $got_a_match = 1;
      $state->{'matching_ids'} .= $other->[0];
      $state->{'last_used'} = $i;
      #print "got a match with $other->[0] at $other->[1] $other->[2]\n";
      last;


    } else {
      #print "no match\n";
    }
    $other_start = $other->[1];
  }
  #print "out of OTHER loop\n";

  return $got_a_match;

}


##########################################
# find the SAGE transcripts that don't overlap the coding exons and pseudogenes
#  &get_unmatched_SAGE(\@coding_transcripts, \@SAGE_transcripts);

sub get_unmatched_SAGE {

  my ($exons_aref, $pseudogenes_aref, $SAGE_transcripts_aref, $transposons_aref, $transposon_exons_aref, $noncoding_transcripts_aref, $rRNA_aref, $chromosome) = @_;

  my @SAGE_transcripts = @{$SAGE_transcripts_aref};

  my %exons_match = &init_match();
  my %pseud_match = &init_match();
  my %trans_match = &init_match();
  my %trane_match = &init_match();
  my %nonco_match = &init_match();
  my %rrna_match = &init_match();

  foreach my $sage (@SAGE_transcripts) { # $SAGE_id, $chrom_start, $chrom_end, $chrom_strand

    my $got_a_match = 0;
  
    if (&match($sage, $exons_aref, \%exons_match)) {
      $got_a_match = 1;
    }

    if (&match($sage, $pseudogenes_aref, \%pseud_match)) {
      $got_a_match = 1;
    }

    if (&match($sage, $transposons_aref, \%trans_match)) {
      $got_a_match = 1;
    }

    if (&match($sage, $transposon_exons_aref, \%trane_match)) {
      $got_a_match = 1;
    }

    if (&match($sage, $noncoding_transcripts_aref, \%nonco_match)) {
      $got_a_match = 1;
    }

    if (&match($sage, $rRNA_aref, \%rrna_match)) {
      $got_a_match = 1;
    }

    # output unmatched SAGE sites to the database
    if (! $got_a_match) {
      my $SAGE_id = $sage->[0];
      my $chrom_start = $sage->[1];
      my $chrom_end = $sage->[2];
      my $chrom_strand = $sage->[3];

      # find the relative abundance of each observation of the SAGE_tags for this SAGE_transcript - this didn't work!
#      my $anomaly_score = SAGE_score($SAGE_id);
      my $anomaly_score = 0.1;	# ... so just set the score to be this
      #print "NOT got a match ANOMALY: $protein_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score\n";
      &output_to_database("UNMATCHED_SAGE", $chromosome, $SAGE_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, '');
    }
  }

}

##########################################
# this never worked properly - not used
#
# find the relative abundance of each observation of the SAGE_tags for this SAGE_transcript
# and multiple by 10,000 to give a score for the SAGE				
# my $anomaly_score = SAGE_score($SAGE_id);
# 
# the sage_tag data looks like:
#
#SAGE_tag : "SAGE:agagaacagt"
#Method   "SAGE_tag"
#Tag_sequence     "CATGAGAGAACAGT"
#Anchoring_enzyme         "NlaIII"
#Tag_length       14
#SAGE_transcript  "2RSSE:agagaacagt"
#Unambiguously_mapped
#Results  "[cgc4805]:dauer" Frequency 3
#Results  "[cgc4805]:dauer" Total_tag_count 132582
#Results  "[cgc4805]:dauer" Relative_abundance 0.000023
#Results  "[cgc4805]:mixed" Frequency 1
#Results  "[cgc4805]:mixed" Total_tag_count 137502
#Results  "[cgc4805]:mixed" Relative_abundance 0.000007
#
# something like this might work:
# at('Results')->right(3)->down(2) and then maybe ->right()->name
#
# or possibly:
# get('Results', 3)->down(2) and then ->right()->name
#
#



#sub SAGE_score {
#
#  my ($SAGE_id) = @_;
#
#  # get SAGE_tag for this SAGE_transcript
#  #print "Ace fetch->$SAGE_id\n";
#  my $SAGE_TR_obj = $ace->fetch("SAGE_transcript" => $SAGE_id);
#  # debug
#  if (! defined $SAGE_TR_obj) {
#    print "Can't fetch SAGE_transcript object for $SAGE_id\n";
#    return (0.1);		# return a default value
#  }
#
#  # get the SAGE_tag
#  my $SAGE_TAG = $SAGE_TR_obj->SAGE_tag;
#  #print "SAGE_TAG $SAGE_TAG\n";
#
#  my $SAGE_TAG_obj = $ace->fetch("SAGE_tag" => $SAGE_TAG);
#  if (! defined $SAGE_TAG_obj) {
#    print "Can't fetch SAGE_tag object for $SAGE_id\n";
#    return (0.1);		# return a default value
#  }
#
#  # debug
#  my @results = $SAGE_TAG_obj->get('Results', 4);
#  foreach my $result_text (@results) {
#    #print "$result_text\n";
#  }
#
#  return (0.1);
#
#  # get the sum of the relative abundances * 10,000 for each library of SAGE_tags
#  my @relative_abundance = $SAGE_TAG_obj->Relative_abundance;
#
#  my $result = 0;
#  foreach my $relative_abundance (@relative_abundance) {
#    $result += ($relative_abundance * 10000);
#    #print "relative_abundance = $relative_abundance\n";
#  }
#
#  return $result;
#}

##########################################
# get GFF files from various locations

sub get_GFF_files {
  my ($chromosome) = @_;

  # eventually this data shuold probably be read from a data file somewhere
  my @GFF_data = 
  (
   {
     directory			=> "/nfs/WWWdev/SANGER_docs/htdocs/Projects/C_elegans/WORMBASE/development_release/GFF/",
     file			=> "CHROMOSOME_${chromosome}.check_intron_*.gff",
     gff_source			=> "",
     gff_type			=> "intron",
     anomaly_type		=> "CONFIRMED_EST_INTRON",
     ID_after			=> "Confirmed_EST\\s+",
     action                     => ["save_to_database"], # anonymous array of actions to perform on each line
   },

   {
     directory			=> "/nfs/WWWdev/SANGER_docs/htdocs/Projects/C_elegans/WORMBASE/development_release/GFF/",
     file			=> "CHROMOSOME_${chromosome}.check_intron_*.gff",
     gff_source			=> "",
     gff_type			=> "intron",
     anomaly_type		=> "CONFIRMED_cDNA_INTRON",
     ID_after			=> "Confirmed_cDNA\\s+",
     action                     => ["save_to_database"],
   },

#   {
#     directory			=> "/nfs/WWWdev/SANGER_docs/htdocs/Projects/C_elegans/WORMBASE/development_release/Checks/",
#     file			=> "",
#     gff_source	       	=> "",
#     gff_type			=> "",
#     anomaly_type		=> "",
#     ID_after			=> "ID tag\\s+",
#     action                    => ["save_to_database", "return_result", "whatever"],
#   },


  );


  foreach my $GFF_data (@GFF_data) {
    &read_GFF_file($GFF_data);
  }


}

##########################################
# use the data in the hash-ref to read in GFF files
# various actions can be performed on the lines to filter and process them

sub read_GFF_file {
  my ($GFF_data) = @_;

  my @result;			# returned result

  my @files = glob($GFF_data->{'directory'} . $GFF_data->{'file'});

  my $score = 1.0;		# default score
  my ($id, $hit_start, $hit_end);

  foreach my $file (@files) {
    open (GFF, "< $file") || die "Can't open $file\n";
    while (my $line = <GFF>) {
      chomp $line;
      if ($line =~ /^\s*$/) {next;}
      if ($line =~ /^#/) {next;}
	  my @f = split /\t/, $line;
	  my ($chromosome, $source, $type, $start, $end, $sense) = ($f[0], $f[1], $f[2], $f[3], $f[4], $f[6]);
	  if ($GFF_data->{'gff_source'} ne "" && $GFF_data->{'gff_source'} ne $source) {next;}
	  if ($GFF_data->{'gff_type'} ne "" && $GFF_data->{'gff_type'} ne $type) {next;}
	  if (exists $GFF_data->{'homology'}) {	# do we need to store the homology data?
	    ($id, $hit_start, $hit_end) = ($f[8] =~ /$GFF_data->{'ID_after'}(\S+)\s+(\d+)\s+(\d+)/);
	    if ($f[5] =~ /\d+/) {$score = $f[5]}; # if the score is numeric, store it
	    #print "got homology data: $hit_start, $hit_end\n";
	  } else {
	    ($id) = ($f[8] =~ /$GFF_data->{'ID_after'}(\S+)/);
	  }
	  if (! defined $id) {next;}
	  $id =~ s/\"//g;	# remove quotes
	  if ($chromosome =~ /CHROMOSOME_(\S+)/) {$chromosome = $1;} # abbreviate chromosome

	  #print "read GFF: $chromosome, $id, $start, $end, $sense\n";

	  # what should we do with the resulting line?
	  foreach my $action (@{$GFF_data->{'action'}}) {
	    #print "action = $action\n";
	    # save this line directly to the database
	    if ($action eq "save_to_database") {
	      &output_to_database($GFF_data->{'anomaly_type'}, $chromosome, $id, $start, $end, $sense, $score, '');
	    }
	    # push this line onto the array to return from this routine
	    if ($action eq "return_result") {
	      #print "GFF line pushed to results array\n";
	      if (exists $GFF_data->{'homology'}) {	# do we need to store the homology data?
		push @result, [$id, $start, $end, $sense, $hit_start, $hit_end, $score];
	      } else {
		push @result, [$id, $start, $end, $sense];		
	      }
	    }
	  }			# end of actions


    }
    close (GFF);
  }
      
  # return saved result sorted by chromosomal start position
  return sort {$a->[1] <=> $b->[1]} @result;

}

##########################################

sub output_to_database {
  my ($anomaly_type, $chromosome, $anomaly_id, $chrom_start, $chrom_end, $chrom_strand, $anomaly_score, $explanation ) = @_;

  # get the clone and lab for this location
  my ($clone, $clone_start, $clone_end) = $coords->LocateSpan($chromosome, $chrom_start, $chrom_end);
  my $lab = $clonelab{$clone};          # get the lab that sequenced this clone
  if ($clone =~ /CHROMOSOME/) {
    $lab = "HX/RW";
  } elsif ($clone =~ /SUPERLINK_(\S\S)/) {
    if ($1 eq "CB") {
      $lab = "HX";
    } else {
      $lab = "RW";
    }
  }
 
  # calculate the window value as blocks of 10 kb
  my $window =  int($chrom_start/10000);

  # need to do a check for a very similar previous record that may
  # need to be overwritten, preserving the status.

  my $nearest_db_key_id = -1;
  my $nearest = 100;	# this size of distance will cause a new record to be inserted if it is not changed to <= 20

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



  # is the distance less than 20 bases?
  # if it is not zero this is probably a move of the anomaly 
  # as a result of genome sequence changes or
  # changes in the blast database size.
  # so we should update the existing record
  if ($nearest <= 20) {
    $mysql->do(qq{ UPDATE anomaly SET   clone="$clone", clone_start=$clone_start, clone_end=$clone_end, centre="$lab", chromosome_start=$chrom_start, chromosome_end=$chrom_end, thing_score=$anomaly_score, explanation="$explanation"   WHERE anomaly_id = $nearest_db_key_id; });
    # NB we do not write the status record for this anomaly_id

  } else {

    # we want a new record inserted
    # write the data to the database
    $db_key_id++;
    $mysql->do(qq{ insert into anomaly values ($db_key_id, "$anomaly_type", "$clone", $clone_start, $clone_end, "$lab", "$chromosome", $chrom_start, $chrom_end, "$chrom_strand", "$anomaly_id", $anomaly_score, "$explanation", $window, 1, NULL); });
    #print "*** inserting new record\n";

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
# This should expand on your brief description above and add details of any options
# that can be used with the program.  Such documentation can be viewed using the perldoc
# command.


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

=item none

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Verbose/Debug mode
 
=back

=over 4

=item -test, Test mode, generate the acefile but do not upload themrun the script, but don't change anything

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back
                                                                                             

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Gary Williams (gw3@sanger.ac.uk)

=back

=cut
