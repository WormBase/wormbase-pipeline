#!/usr/bin/env perl
#
# assign_orientation.pl                       
# 
# by Gary Williams                        
#
# This assigns orientation to BLATted sequences which have no
# orientation yet.  It uses best splice site and overlap with
# transcripts to find the most probable orientation.
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2013-10-14 10:16:24 $      

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Modules::Overlap;
use Modules::PWM;



######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase,$noload);
my ($all, $species, $gff_directory, $gff_file, $gff_source, $gff_type, $ID_after);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
            "noload"     => \$noload,
	    "all"        => \$all, # do all EST sequences, not just the ones with no orientation
	    "species:s"  => \$species, # the default is elegans
	    "gff_directory:s" => \$gff_directory,	# stuff to specify a specific GFF file to search
	    "gff_file:s"      => \$gff_file,
	    "gff_source:s"    => \$gff_source,
	    "gff_type:s"      => \$gff_type,
	    "ID_after:s"      => \$ID_after);


$species = 'elegans' unless $species;

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


##########################
# MAIN BODY OF SCRIPT
##########################

my %dna_entry;		# store for dna sequences for when reading remanei-style multi-entry sequence file

my $database = $wormbase->autoace;

# output file name
my $version = $wormbase->get_wormbase_version_name();
my $output = $wormbase->wormpub . "/CURATION_DATA/assign_orientation." . $version . ".ace";
if ($wormbase->species ne 'elegans') {
  $output = $wormbase->acefiles . "/assign_orientation.ace";
}

print "find ignored sequences\n" if ($verbose);
my %ignore = &get_Ignore();
print "find sequences without Database\n" if ($verbose);
my %not_db = &get_NOT_Database();
print "find EST orientation\n" if ($verbose);
my %Show_in_reverse_orientation = $wormbase->FetchData('estorientation');

my $pwm = PWM->new;

# loop through the chromosomes
foreach my $chromosome ($wormbase->get_chromosome_names(-mito => 1, -prefix => 1)) {
  print STDERR "Doing $chromosome\n";

# get the Overlap object
my $ovlp = Overlap->new($database, $wormbase);


  my %splice_ok = ();
  my %splice_reverse = ();
  my %gene_ok = ();
  my %gene_reverse = ();

  # read in the lists of GFF data or the specified GFF file
  print "Reading GFF files\n" if ($verbose);
  my @est_hsp;

  # if a specific GFF file is given, read that in instead of the normal EST etc sequences to check
  if ($gff_directory && $gff_file && $gff_source && $gff_type && $ID_after) {
    print "gff_file = $gff_file\n";
    my %GFF_data = (
		    file=>$gff_directory . "/" . $chromosome . "_" . $gff_file, 
		    gff_source=>$gff_source, 
		    gff_type=>$gff_type,
		    ID_after=>$ID_after, 
		    reverse_orientation=>1,
		    homology=>1
		    );
    print "gff_file in hash=" . $GFF_data{file} . "\n";
    @est_hsp = $ovlp->read_GFF_file($chromosome, \%GFF_data);
  } else {
    @est_hsp = $ovlp->get_EST_BEST($chromosome); # get main list of EST entries to search with
    push @est_hsp, $ovlp->get_mRNA_BEST($chromosome); # and the list of mRNAs ...
    push @est_hsp, $ovlp->get_Trinity_BEST($chromosome); # and the list of Trinity transcripts ...
    push @est_hsp, $ovlp->get_OST_BEST($chromosome) if $species eq 'elegans'; # and add on the list of OSTs ..
    push @est_hsp, $ovlp->get_RST_BEST($chromosome) if $species eq 'elegans'; # and the list of RSTs
  }

  # check we have some ESTs/mRNAs etc. 
  if (! scalar @est_hsp) {
    print "No ESTs found in $chromosome\n" if ($verbose);
    next;
  }

  @est_hsp = sort {$a->[1] <=> $b->[1] or $a->[2] <=> $b->[2]} @est_hsp; # and ensure it is sorted by chromosomal start,end position

  my @est = $ovlp->get_span(@est_hsp); # change the ESTs from HSPs to start-to-end span
  my @cds  = $ovlp->get_curated_CDS_exons($chromosome); # secondary list we search against
  my @pseud = $ovlp->get_pseudogene($chromosome); # secondary list we search against
  my @rrna = $ovlp->get_rRNA_exons($chromosome); # secondary list we search against
  my @ncrna = $ovlp->get_ncRNA($chromosome); # secondary list we search against
  #my @jigsaw = $ovlp->get_jigsaw_exons($chromosome); # secondary list we search against
  my @trans = $ovlp->get_Coding_transcript_exons($chromosome); # secondary list we search against


  # get ESTs with no orientation (or get all ESTs if $all or a GFF file is specified)
  my @no53_est;
  foreach my $est (@est) {

    print "$est->[0] ignored\n" if ($verbose && exists $ignore{$est->[0]});
    if (exists $ignore{$est->[0]}) {next;} # withdrawn, rRNA or chimaeric

    print "$est->[0] no Database tag\n" if ($verbose && exists $not_db{$est->[0]});
    if (exists $not_db{$est->[0]}) {next;} # estorientation hash file will not have this sequence's orientation even if it is set
    print "$est->[0] already in reverse orientation\n" if ($verbose && exists $Show_in_reverse_orientation{$est->[0]});
    if ($all || !exists $Show_in_reverse_orientation{$est->[0]} || $gff_file) {
      print "$est->[0] to be analysed\n" if ($verbose);
      push @no53_est, $est;
    }
  }
  print "Have ", scalar @no53_est, " sequences with no orientation\n" if ($verbose);

  # get the corresponding HSPs of these ESTs
  my %ids;
  foreach my $est (@no53_est) {
    $ids{$est->[0]} = 1;
  }
  my @no53_est_hsp;
  foreach my $est (@est_hsp) {
    if (exists $ids{$est->[0]}) {
      push @no53_est_hsp, $est;
    }
  }

  # sort the HSPs by name and then by chromosomal position
  my @sorted_hsp = sort {$a->[0] cmp $b->[0]
		   or
		 $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2]} @no53_est_hsp;


  # read the chromosmal sequence
  $log->write_to("\n\nChromosome $chromosome\n");
  my $seq = read_chromosome($chromosome);

  # get the 5' and 3' splice sites of the EST HSPs
  my $prev_id = "";
  my ($id, $start, $end, $sense, $hit_start, $hit_end, $score);
  my ($prev_start, $prev_end, $prev_sense, $prev_hit_start, $prev_hit_end, $prev_score);
  my $score_5;
  my $score_3;
  my $forward;
  my $reverse;
  foreach my $hsp (@sorted_hsp) {
    if ($hsp->[0] ne $prev_id) {	# new ID
      ($id, $start, $end, $sense, $hit_start, $hit_end, $score) = @{$hsp};
      $prev_id = $id;
    } else {			# same as the previous ID
      ($prev_start, $prev_end, $prev_sense, $prev_hit_start, $prev_hit_end, $prev_score) = ($start, $end, $sense, $hit_start, $hit_end, $score);
      ($id, $start, $end, $sense, $hit_start, $hit_end, $score) = @{$hsp};
      # see if we have a clean splice with no missed EST bases
      if (abs($hit_start - $prev_hit_end) == 1) {
	print "$id $sense\n" if ($verbose);

	# check the splice sites
	# forward sense
	$score_5 = $pwm->splice5($seq, $prev_end-1, '+'); # -1 to convert from sequence pos to perl string coords
	$score_3 = $pwm->splice3($seq, $start-2, '+'); # -1 to convert from sequence pos to perl string coords
	$forward = $score_5 + $score_3;

	# reverse sense
	$score_5 = $pwm->splice5($seq, $start-1, '-'); # -1 to convert from sequence pos to perl string coords
	$score_3 = $pwm->splice3($seq, $prev_end, '-'); # -1 to convert from sequence pos to perl string coords
	$reverse = $score_5 + $score_3;

	# good difference between the intron splice scores on the two strands
	# and one or the other has good splice scores
	if (abs($reverse - $forward) > 1 && ($reverse > 1 || $forward > 1)) {
	  print "$id $forward $reverse " if ($verbose);
	  if (($sense eq '+' && $reverse > 1) || ($sense eq '-' && $forward > 1)) {
	    print "$id $forward $reverse  ************* REVERSE SENSE REQUIRED (splice)\n" if ($verbose);
	    $splice_reverse{$id} = 1;
	  } elsif (($sense eq '+' && $forward > 1) || ($sense eq '-' && $reverse > 1)) {
	    $splice_ok{$id} = 1;
	  }
	}

      }
    }
  }

  if (!$all) {			# don't want to do all this error-prone stuff on the complete set of sequences

    # set up the overlap compare objects for the secondary lists
    my $cds_obj = $ovlp->compare(\@cds);
    my $pseud_obj = $ovlp->compare(\@pseud);
    my $rrna_obj = $ovlp->compare(\@rrna);
    my $ncrna_obj = $ovlp->compare(\@ncrna);
    #my $jigsaw_obj = $ovlp->compare(\@jigsaw);
    my $trans_obj = $ovlp->compare(\@trans);

    print "Searching\n" if ($verbose);
    foreach my $est (@no53_est) { 
      # $est is a ref to:
      # ($est_id, $chrom_start, $chrom_end, $chrom_strand, $hit_start, $hit_end, $score)

      my $id = $est->[0];
      print "next EST to check = $id\n" if ($verbose);

      # look to see if the EST has any matches to CDS and/or pseudogenes
      my @cds_matches = $cds_obj->match($est);
      print "Have ", scalar @cds_matches, " overlaps to cds\n" if ($verbose);
      my $prop_same = 0;		# proportion matching gene on same sense
      my $prop_other = 0;		# proportion matching gene on other sense
      if (@cds_matches) {
	my @senses = $cds_obj->matching_sense;
	for (my $i=0; $i < @senses; $i++) {
	  my ($p1,$p2) = $cds_obj->matching_proportions($cds_matches[$i]);
	  if ($senses[$i] == 1) {
	    if ($prop_same < $p1) {$prop_same = $p1;}
	  } else {
	    if ($prop_other < $p1) {$prop_other = $p1;}
	  }
	}
      }

      my @pseud_matches = $pseud_obj->match($est);
      print "Have ", scalar @pseud_matches, " overlaps to pseudogenes\n" if ($verbose);
      if (@pseud_matches) {
	my @senses = $pseud_obj->matching_sense;
	for (my $i=0; $i < @senses; $i++) {
	  my ($p1,$p2) = $pseud_obj->matching_proportions($pseud_matches[$i]);
	  if ($senses[$i] == 1) {
	    if ($prop_same < $p1) {$prop_same = $p1;}
	  } else {
	    if ($prop_other < $p1) {$prop_other = $p1;}
	  }
	}
      }

      my @rrna_matches = $rrna_obj->match($est);
      print "Have ", scalar @rrna_matches, " overlaps to rRNAs\n" if ($verbose);
      if (@rrna_matches) {
	my @senses = $rrna_obj->matching_sense;
	for (my $i=0; $i < @senses; $i++) {
	  my ($p1,$p2) = $rrna_obj->matching_proportions($rrna_matches[$i]);
	  if ($senses[$i] == 1) {
	    if ($prop_same < $p1) {$prop_same = $p1;}
	  } else {
	    if ($prop_other < $p1) {$prop_other = $p1;}
	  }
	}
      }

      my @ncrna_matches = $ncrna_obj->match($est);
      print "Have ", scalar @ncrna_matches, " overlaps to ncRNAs\n" if ($verbose);
      if (@ncrna_matches) {
	my @senses = $ncrna_obj->matching_sense;
	for (my $i=0; $i < @senses; $i++) {
	  my ($p1,$p2) = $ncrna_obj->matching_proportions($ncrna_matches[$i]);
	  if ($senses[$i] == 1) {
	    if ($prop_same < $p1) {$prop_same = $p1;}
	  } else {
	    if ($prop_other < $p1) {$prop_other = $p1;}
	  }
	}
      }

      #my @jigsaw_matches = $jigsaw_obj->match($est);
      #print "Have ", scalar @jigsaw_matches, " overlaps to jigsaw\n" if ($verbose);
      #if (@jigsaw_matches) {
	#my @senses = $jigsaw_obj->matching_sense;
	#for (my $i=0; $i < @senses; $i++) {
	#  my ($p1,$p2) = $jigsaw_obj->matching_proportions($jigsaw_matches[$i]);
	#  if ($senses[$i] == 1) {
	#    if ($prop_same < $p1) {$prop_same = $p1;}
	#  } else {
	#    if ($prop_other < $p1) {$prop_other = $p1;}
	#  }
	#}
      #}

      # as a last resort, only if we have no other information, we
      # take note of the transcripts that the EST overlaps
      my $trans_prop_same = 0;
      my $trans_prop_other = 0;
      my @trans_matches = $trans_obj->match($est);
      print "Have ", scalar @trans_matches, " overlaps to Coding transcript exons\n" if ($verbose);
      if ($prop_same == 0 && $prop_other == 0) {
	if (@trans_matches) {
	  my @senses = $trans_obj->matching_sense;
	  for (my $i=0; $i < @senses; $i++) {
	    my ($p1,$p2) = $trans_obj->matching_proportions($trans_matches[$i]);
	    if ($senses[$i] == 1) {
	      if ($trans_prop_same < $p1) {$trans_prop_same = $p1;}
	    } else {
	      if ($trans_prop_other < $p1) {$trans_prop_other = $p1;}
	    }
	  }
	}
      }


      # work out the orientation
      # if there is splice evidence for the orientation, then this
      # overrides any other evidence
      # the overlap to coding exons overrides evidence from 
      # overlap to coding transcript exons
      if (!$splice_ok{$id} && !$splice_reverse{$id}) {
	if ($prop_other > $prop_same) {
	  print "$id $prop_same $prop_other  ************* REVERSE SENSE REQUIRED (CDS)\n" if ($verbose);
	  $gene_reverse{$id} = 1;
	} elsif ($prop_same > $prop_other) {
	  $gene_ok{$id} = 1;
	} else {			# no splice or exon evidence, try transcript overlap
	  if ($trans_prop_other > $trans_prop_same) {
	    print "$id $trans_prop_same $trans_prop_other  ************* REVERSE SENSE REQUIRED (transcript)\n" if ($verbose);
	    $gene_reverse{$id} = 1;
	  } else {
	    print "$id $trans_prop_same $trans_prop_other (transcript)\n" if ($verbose);
	    $gene_ok{$id} = 1;
	  }
	}
      }
    }

  }
 
  # now write the ACE file for this chromosome
  open (ACE, ">> $output") || die "Can't open file $output\n";
  # we are happy that the existing orientation of these is OK
  # if no existing orientation, we make the default: EST_5
  my $count_out = 0;
  $log->write_to("\nThe following are confirmed in their default orientation\n") if (scalar(keys %splice_ok) + scalar(keys %gene_ok));
  foreach my $id (keys %splice_ok, keys %gene_ok) {
    if (!exists $Show_in_reverse_orientation{$id}) {
      print ACE "\nSequence : $id\n";
      print ACE "EST_5\n";
      $log->write_to("$id\t");
      if (($count_out % 5) == 4) {$log->write_to("\n");}
      $count_out++;
    }
  }

  # we need to swap the orientation of these
  # they may already have an orientation that needs to be deleted
  $count_out = 0;
  $log->write_to("\nThe following have their orientation reversed\n") if (scalar(keys %splice_reverse) + scalar(keys %gene_reverse));
  foreach my $id (keys %splice_reverse, keys %gene_reverse) {
    if (!exists $Show_in_reverse_orientation{$id} || $Show_in_reverse_orientation{$id} == 5) {
      print ACE "\nSequence : $id\n";
      print ACE "-D EST_5\n";

      print ACE "\nSequence : $id\n";
      print ACE "EST_3\n";
      print ACE "Show_in_reverse_orientation\n";
    } elsif ($Show_in_reverse_orientation{$id} == 3) {
      print ACE "\nSequence : $id\n";
      print ACE "-D EST_3\n";
      print ACE "-D Show_in_reverse_orientation\n";

      print ACE "\nSequence : $id\n";
      print ACE "EST_5\n";
    }
    $log->write_to("$id\t");
    if (($count_out % 5) == 4) {$log->write_to("\n");}
    $count_out++;
  }
  close(ACE);

}


# The elegans orientation data is read into camace during the camace
# merge-split procedure.

# If the species is not 'elegans' then we need to load the orientation
# ace file into the database in ~wormpub/DATABASES/$species This will
# ensure that the orientations will be in the database ready for the
# next Build procedure. We do NOT want these orientations loaded into
# the current BUILD database because we have already done teh GFF
# dumps and the orientations would be inconsistent.
if ($wormbase->species ne 'elegans' and not $noload) {
  my $database = $wormbase->database($wormbase->species);
  $wormbase->load_to_database($database, $output, 'assign_orientation.pl', $log, undef, 1);
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
# read the chromosmal sequence
#  my $seq = read_chromosome($chromosome);

sub read_chromosome {

  my $chromosome = shift;

  # if we have already read in the sequence entries, return the one for this chromosome
  if (exists $dna_entry{$chromosome}) {return $dna_entry{$chromosome};}

  my $seq_file = $wormbase->chromosomes . "/$chromosome.dna";
  my $seq = &read_file($seq_file);

  if (! defined $seq) {
    $seq = &read_entry($wormbase->genome_seq, $chromosome);
  }

  return $seq;

}

##########################################
# read file

sub read_file {
  my ($file) = @_;

  my $seq;

  # try to open a single-chromosome entry file e.g. for elegans, briggsae
  if (open (SEQ, $file)) {
    $/ = "";
    $seq = <SEQ>;
    close SEQ;
    $/ = "\n";
    
    $seq =~ s/>.*\n//;           # remove one title line
    $seq =~ s/\n//g;
  } else {
    print "Can't open the dna file for $file : $!\n" if ($verbose);
  }

  return $seq
}


##########################################
# read entry in a file if the file consists of lots of entries

sub read_entry {
  my ($file, $chromosome) = @_;

  # if we have already read in the sequence entries, return the one for this chromosome
  if (exists $dna_entry{$chromosome}) {return $dna_entry{$chromosome};}

  print "Trying: $file\n" if ($verbose);
  if (! -e $file) {
    die "Can't find chromosome sequence file for $file, $chromosome\n";
  }

  print "Reading DNA sequence entries\n";
  my $seq;

  if (open (SEQ, $file)) {
    my $entry;
    while (my $line = <SEQ>) {
      chomp $line;
      if ($line =~ /^>(\S+)/) {
	if (defined $entry) {	# store the previous sequence entry
	  $dna_entry{$entry} = $seq;
	}
	$seq = "";
	$entry = $1;
      } else {
	$seq .= $line;
      }
    }
    if (defined $entry) {	# store the last entry
      $dna_entry{$entry} = $seq;
    }
    close SEQ;
  } else {
    die "Can't open the dna file for $file : $!\n";
  }

  if (! exists $dna_entry{$chromosome}) {$log->log_and_die("No DNA sequence was found for $chromosome\n");}
  return $dna_entry{$chromosome};
}


##########################################
# get the Sequence objects with the tag 'Ignore'

sub get_Ignore {

  my %result;

  my $cmd = "query find Sequence\nwhere Ignore\nlist";
  my $tace = $wormbase->tace;
  my $database = $wormbase->autoace;

  open (TACE, "echo '$cmd' | $tace $database |");
  while (my $line = <TACE>) {
    chomp $line;
    next if ($line =~ /acedb\>/);
    next if ($line =~ /\/\//);

    my ($id) = ($line =~ /(\S+)/); # remove any preceding or trailing space
    if (defined $id) {
      $result{$id} = 1;
    }
  }
  close TACE;


  return %result;
}

##########################################
# get the Sequence objects without the tag 'Database'

sub get_NOT_Database {

  my %result;

  my $cmd = "query find Sequence\nwhere NOT Database\nlist";
  my $tace = $wormbase->tace;
  my $database = $wormbase->autoace;

  open (TACE, "echo '$cmd' | $tace $database |");
  while (my $line = <TACE>) {
    chomp $line;
    next if ($line =~ /acedb\>/);
    next if ($line =~ /\/\//);

    my ($id) = ($line =~ /(\S+)/); # remove any preceding or trailing space
    if (defined $id) {
      $result{$id} = 1;
    }
  }
  close TACE;


  return %result;
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

=head2 NAME - assign_orientation.pl

=head1 USAGE

=over 4

=item  assign_orientation.pl [-options]

=back

This script assigns the EST_3 or EST_5 tags to Sequence objects after
they have been BLATed to the genome.

It reads the dumped GFF files for EST, mRNA and OST and looks at the
orientation from the COMMON_DATA/estorientation file.  Optionally, the
data can be read instead from a specified GFF file.

If the -all flag is NOT set and a specific GFF file is NOT specified,
then only the sequences with no EST_5 and no EST_3 tag set are
examined. For these sequences, if the orientation cannot be confirmed
from the splice site score then the greatest overlaps with CDS exons,
rRNA exons, pseudogenes and jigsaw exons are used to assign the
orientation. If the orientation still cannot be assigned, then the
greatest overlap with transcript exons are used. If the evidence is
still lacking and there is no other tag already set then the default
EST_5 tag is set.

If the -all flag is set or if a specific GFF file is specified then it
examines every sequence to see if the alignment to the genome
indicates (using the splice-site score) whether the orientation should
be changed. The more error-prone overlap evidence will not be used to
assign orientation.

If a specific GFF input file is specified, then all sequences in the
file are examined and if there is no evidence of the orientation from
the splice site scores, the overlap with exons will be used to assign
the orientation.

assign_orientation.pl MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

=back

=over 4

=item -all, do all EST sequences, not just the ones with no orientation tags. If this option is set then only splice site scores will be used to assign the orientation, the more error-prone overlap evidence will be ignored.

=back

=over 4

=item -gff_directory, if a specific GFF file is to be read in, this is the directory it is in.

=back

=over 4

=item -gff_file, if a specific GFF file is to be read in, this is the file name (without the chromosome prefix and number eg 'BLAT_mRNA_OTHER.gff').

=back

=over 4

=item -gff_source, if a specific GFF file is to be read in, this is the source name in the GFF lines to use.

=back

=over 4

=item -gff_type, if a specific GFF file is to be read in, this is the type name in the GFF lines to use.

=back

=over 4

=item -ID_after, if a specific GFF file is to be read in, the ID of the sequence can be found after this regular expression, typically something like: 'Target\s+"Sequence:'

=back

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

=item GFF files of the sequence alignments. COMMON_DATA/estorientation hash file.

=back

=head1 AUTHOR

=over 4

=item Gary Williams (gw3@sanger.ac.uk)

=back

=cut
