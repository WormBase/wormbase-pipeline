#!/software/bin/perl -w
#
# script to find BLAT alignments that match to more than one CDS and
# set the Ignore tag in them so that they are not used in the
# transcript_builder script.
# 
# by Gary Williams
#
# Last updated by: $Author: pad $     
# Last updated on: $Date: 2013-08-14 12:19:59 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Modules::Overlap;

#use Ace;
#use Sequence_extract;
#use Coords_converter;
#use Feature_mapper;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $database, $output, $species, $chromosome, $noload);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "database:s" => \$database,
            "species:s"  => \$species,
	    "chromosome:s" => \$chromosome, # specify a single chromosome with prefix ('CHROMOSOME_II') for debugging purposes
	    "noload"     => \$noload, # specify that the resulting ace file should not be loaded into the database for debugging purposes
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

if (! defined $database) {$database = $wormbase->autoace}
if (defined $chromosome) {
  $output = $wormbase->acefiles . "/ignore_bad_blat_alignments_${chromosome}.ace";
} else {
  $output = $wormbase->acefiles . "/ignore_bad_blat_alignments.ace";
}

#################################

# Set up top level base directories

my $chromosomes_dir = $wormbase->chromosomes; # AUTOACE CHROMSOMES
my $gff_dir         = $wormbase->gff;         # AUTOACE GFF
my $gff_splits_dir  = $wormbase->gff_splits;  # AUTOACE GFF SPLIT

my %mol_types = ( 'elegans'          => [qw( EST mRNA OST RST Trinity Nanopore)],
                  'briggsae'         => [qw( mRNA EST Trinity)],
                  'remanei'          => [qw( mRNA EST)],
                  'brenneri'         => [qw( mRNA EST)],
                  'japonica'         => [qw( mRNA EST Trinity)],
                  'brugia'           => [qw( mRNA EST Trinity)],
                  'pristionchus'     => [qw( mRNA EST)],
                  'ovolvulus'        => [qw( mRNA EST Trinity)],
                  'sratti'           => [qw( mRNA EST)],
                  'tmuris'           => [qw( mRNA EST Trinity IsoSeq)],
                );


##########################

$log->write_to("The following BLAT alignments span the introns of two or more CDS structures.\n\n");
$log->write_to("They should be inspected to see if the gene models should be merged or if the transcripts are chimeric or incompletely spliced operon transcripts\n\n");
$log->write_to("The following transcripts have had the Ignore tag set and will not be used for building Coding_transcripts in transcript_builder.pl\n\n");

$species = $wormbase->species;
open (ACE, ">$output") || die "Can't open $output\n";

my $regex = $wormbase->seq_name_regex;

my @chromosomes;

if (defined $chromosome) {
  @chromosomes = ($chromosome);
} else {
  @chromosomes = $wormbase->get_chromosome_names(-mito => 0, -prefix => 1);
}

foreach my $chromosome (@chromosomes) {
  $log->write_to("\n\nChromosome: $chromosome\n");

  my $ovlp = Overlap->new($database, $wormbase);

  my @est_introns;
  my @rst_introns;
  my @ost_introns;
  my @mrn_introns;
  my @tri_introns;
  my @iso_introns;
  my @nan_introns;
  my @nan_introns_confirmed;

  my $est_match;
  my $rst_match;
  my $ost_match;
  my $mrn_match;
  my $tri_match;
  my $iso_match;
  my $nan_match;

  if (grep /^EST$/, @{$mol_types{$species}}) {@est_introns = $ovlp->get_intron_from_exons($ovlp->get_EST_BEST($chromosome))};
  if (grep /^RST$/, @{$mol_types{$species}}) {@rst_introns = $ovlp->get_intron_from_exons($ovlp->get_RST_BEST($chromosome))};
  if (grep /^OST$/, @{$mol_types{$species}}) {@ost_introns = $ovlp->get_intron_from_exons($ovlp->get_OST_BEST($chromosome))};
  if (grep /^mRNA$/, @{$mol_types{$species}}) {@mrn_introns = $ovlp->get_intron_from_exons($ovlp->get_mRNA_BEST($chromosome))};
  if (grep /^Trinity$/, @{$mol_types{$species}}) {@tri_introns = $ovlp->get_intron_from_exons($ovlp->get_Trinity_BEST($chromosome))};
  if (grep /^IsoSeq$/, @{$mol_types{$species}}) {@iso_introns = $ovlp->get_intron_from_exons($ovlp->get_IsoSeq_BEST($chromosome))};
  if (grep /^Nanopore$/, @{$mol_types{$species}}) {
    @nan_introns = $ovlp->get_intron_from_exons($ovlp->get_Nanopore_BEST($chromosome));
    my @RNASeq_introns = $ovlp->get_RNASeq_splice($chromosome);
    # get only those Nanopore introns that are confirmed by the RNASeq introns
    my $RNASeq_confirmed_match = $ovlp->compare(\@RNASeq_introns, exact_match => 1, same_sense => 0); # non-canonical splice site RNASeq introns are sometimes placed on the wrong' strand, so compare to both strands
    my %Nano_not_confirmed;
    foreach my $nano (@nan_introns) {
# debug
#      if ($nano->[0] eq 'male_Nanopore_Roach_126350') {
#      if ($nano->[0] ne 'EM_Nanopore_Li_212889') {
#	next;
#    
#  }

      if (! $RNASeq_confirmed_match->match($nano)) {
	$Nano_not_confirmed{$nano->[0]} = 1; # note the name of any Nanopore intron that does not have any matching RNASeq introns 
      }
    }
    foreach my $nano (@nan_introns)  { # get the introns of Nanaopore transcripts that have all their introns confirmed
      if (!exists $Nano_not_confirmed{$nano->[0]}) {
	push @nan_introns_confirmed, $nano;
      } else {
#	$log->write_to("$nano->[0] has an intron with no RNASeq confirmation\n");
	print ACE "\n\n";
	print ACE "Sequence : $nano->[0]\n";
	print ACE "Ignore Remark \"Has an intron with no RNASeq supporting evidence\"\n";
	print ACE "Ignore Inferred_automatically \"ignore_bad_blat_alignments\"\n";
      }
    }
  }
 
  my @CDS_introns = $ovlp->get_curated_CDS_introns($chromosome);

  if (grep /^EST$/, @{$mol_types{$species}}) {$est_match = $ovlp->compare(\@est_introns, exact_match => 1, same_sense => 0)};  # exact match to either sense
  if (grep /^RST$/, @{$mol_types{$species}}) {$rst_match = $ovlp->compare(\@rst_introns, exact_match => 1, same_sense => 0)};  # exact match to either sense
  if (grep /^OST$/, @{$mol_types{$species}}) {$ost_match = $ovlp->compare(\@ost_introns, exact_match => 1, same_sense => 0)};  # exact match to either sense
  if (grep /^mRNA$/, @{$mol_types{$species}}) {$mrn_match = $ovlp->compare(\@mrn_introns, exact_match => 1, same_sense => 0)};  # exact match to either sense
  if (grep /^Trinity$/, @{$mol_types{$species}}) {$tri_match = $ovlp->compare(\@tri_introns, exact_match => 1, same_sense => 0)};  # exact match to either sense
  if (grep /^IsoSeq$/, @{$mol_types{$species}}) {$iso_match = $ovlp->compare(\@iso_introns, exact_match => 1, same_sense => 0)};  # exact match to either sense
  if (grep /^Nanopore$/, @{$mol_types{$species}}) {$nan_match = $ovlp->compare(\@nan_introns_confirmed, exact_match => 1, same_sense => 0)};  # exact match to either sense


  my %overlapping_hsps = (); # EST/RST/OST/mRNA/Trinity/IsoSeq/Nanopore transcripts that match a CDS, keyed by transcript name, value is array of matching CDSs

  foreach my $cds (@CDS_introns) {

    my ($cds_id) = ($cds->[0] =~ /($regex)/); # get just the sequence name
    

    if ((grep /^EST$/, @{$mol_types{$species}}) && $est_match->match($cds)) {
      my @ids = $est_match->matching_IDs;
      foreach my $id (@ids) {
	$overlapping_hsps{$id}{$cds_id} = 1;
      }
    } elsif ((grep /^RST$/, @{$mol_types{$species}}) && $rst_match->match($cds)) {
      my @ids = $rst_match->matching_IDs;
      foreach my $id (@ids) {
	$overlapping_hsps{$id}{$cds_id} = 1;
      }

    } elsif ((grep /^OST$/, @{$mol_types{$species}}) && $ost_match->match($cds)) {
      my @ids = $ost_match->matching_IDs;
      foreach my $id (@ids) {
	$overlapping_hsps{$id}{$cds_id} = 1;
      }

    } elsif ((grep /^mRNA$/, @{$mol_types{$species}}) && $mrn_match->match($cds)) {
      my @ids = $mrn_match->matching_IDs;
      foreach my $id (@ids) {
	$overlapping_hsps{$id}{$cds_id} = 1;
      }

    } elsif ((grep /^Trinity$/, @{$mol_types{$species}}) && $tri_match->match($cds)) {
      my @ids = $tri_match->matching_IDs;
      foreach my $id (@ids) {
	$overlapping_hsps{$id}{$cds_id} = 1;
      }

    } elsif ((grep /^IsoSeq$/, @{$mol_types{$species}}) && $iso_match->match($cds)) {
      my @ids = $iso_match->matching_IDs;
      foreach my $id (@ids) {
	$overlapping_hsps{$id}{$cds_id} = 1;
      }

    } elsif ((grep /^Nanopore$/, @{$mol_types{$species}}) && $nan_match->match($cds)) {
      my @ids = $nan_match->matching_IDs;
      foreach my $id (@ids) {
	$overlapping_hsps{$id}{$cds_id} = 1;
      }

    }

  }
  
  # now look for transcripts that matched more than one CDS
  foreach my $trans (keys %overlapping_hsps) {
    my @cds = keys %{$overlapping_hsps{$trans}};
    if (scalar @cds > 1) {
      $log->write_to("$trans matches @cds\n");
      print ACE "\n\n";
      print ACE "Sequence : $trans\n";
      print ACE "Ignore Remark \"matches more than one CDS: @cds\"\n";
      print ACE "Ignore Inferred_automatically \"ignore_bad_blat_alignments\"\n";
    }
  }
  
  print ACE "\n\n\n// FINISHED $chromosome\n";
}
print ACE "\n// Finished - Closing file - normal exit\n";
close (ACE);

$wormbase->load_to_database($database, $output, "ignore_bad_blat_alignments", $log) if (!$noload);







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
