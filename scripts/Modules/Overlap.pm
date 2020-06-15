#!/software/bin/perl -w
#
# Overlap.pm
# 
# by Gary Williams                        
#
# Do fast overlap matching of positions of two sets of things.
#
# Last updated by: $Author: pad $     
# Last updated on: $Date: 2014-01-23 13:36:07 $      

=pod

=head1 NAME

 Overlap

=head1 SYNOPSIS

    
    # get the Overlap object
    $wormbase = Wormbase->new( -debug => $debug, -test => $debug);

    my $ovlp = Overlap->new($database, $wormbase);

    # loop through the chromosomes
    foreach my $chromosome ($wormbase->get_chromosome_names(-mito => 1, -prefix => 1)) {

      # read in the lists of GFF data
      my @est_hsp = $ovlp->get_EST_BEST($chromosome);     # main list of entries to search with
      my @est = $ovlp->get_span(@est_hsp); # change the ESTs from HSPs to start-to-end span

      my @cds  = $ovlp->get_curated_CDS($chromosome); # secondary list we search against
      my @pseud = $ovlp->get_pseudogene($chromosome); # secondary list we search against

      # set up the overlap compare objects for the secondary lists
      my $cds_obj = $ovlp->compare(\@cds, same_sense => 1);
      my $pseud_obj = $ovlp->compare(\@pseud, same_sense => 1);

      # now search for overlaps to each EST on this chromosome
      foreach my $est (@est) { 
        # $est is a ref to:
        # ($est_id, $chrom_start, $chrom_end, $chrom_strand, $hit_start, $hit_end, $score)

        print "next EST to check = $est->[0]\n";

        # look to see if the EST has any matches to CDS and/or pseudogenes
        my @cds_matches = $cds_obj->match($est);
        print "Have ", scalar @cds_matches, " overlaps to cds\n";
        my @ids = $cds_obj->matching_IDs;
	print "Matching IDs: @ids\n";
        my @pseud_matches = $pseud_obj->match($est);
        print "Have ", scalar @pseud_matches, " overlaps to pseudogenes\n";
      }
    }
    


Initialise the Overlap object using $wormbase object to find the database

my $ovlp = Overlap->new($wormbase);



Routines to read GFF files from the database defined in $wormbase

@list = $ovlp->read_GFF_file(%data)

@list = $ovlp->get_EST_BEST($chromosome)
@list = $ovlp->get_OST_BEST($chromosome)
@list = $ovlp->get_mRNA_BEST($chromosome)
@list = $ovlp->get_Trinity_BEST($chromosome)
@list = $ovlp->get_IsoSeq_BEST($chromosome)
@list = $ovlp->get_Nanopore_BEST($chromosome)
@list = $ovlp->get_ncRNA_BEST($chromosome)
@list = $ovlp->get_curated_CDS($chromosome)
@list = $ovlp->get_curated_CDS_exons($chromosome)
@list = $ovlp->get_curated_CDS_introns($chromosome)
@list = $ovlp->get_Coding_transcripts($chromosome)
@list = $ovlp->get_Coding_transcript_exons($chromosome)
@list = $ovlp->get_Coding_transcript_introns($chromosome)
@list = $ovlp->get_pseudogene($chromosome)
@list = $ovlp->get_rRNA_transcripts($chromosome)
@list = $ovlp->get_rRNA_exons($chromosome)
@list = $ovlp->get_genefinder_transcripts($chromosome)
@list = $ovlp->get_genefinder_exons($chromosome)
@list = $ovlp->get_twinscan_transcripts($chromosome)
@list = $ovlp->get_twinscan_exons($chromosome)
@list = $ovlp->get_genelets($chromosome)
@list = $ovlp->get_RNASeq_splice($chromosome)
@list = $ovlp->get_Aggregated_CDS_exons($chromosome)
@list = $ovlp->get_Aggregated_CDS_introns($chromosome)
@list = $ovlp->get_mgene_transcripts($chromosome)
@list = $ovlp->get_mgene_exons($chromosome)
@list = $ovlp->get_transposons($chromosome)
@list = $ovlp->get_transposon_exons($chromosome)
@list = $ovlp->get_polyA_signal($chromosome)
@list = $ovlp->get_polyA_site($chromosome)
@list = $ovlp->get_TSL_SL1($chromosome)
@list = $ovlp->get_TSL_SL2($chromosome)
@list = $ovlp->get_SAGE_tags($chromosome)
@list = $ovlp->get_blastx_homologies($chromosome)
@list = $ovlp->get_waba_coding($chromosome)
@list = $ovlp->get_repeatmasked($chromosome)
@list = $ovlp->get_5_UTRs($chromosome)
@list = $ovlp->get_3_UTRs($chromosome)
@list = $ovlp->get_Genes($chromosome)
@list = $ovlp->get_Operons($chromosome)
@list = $ovlp->get_miRNA_mature($chromosome);
@list = $ovlp->get_miRNA_primary($chromosome);
@list = $ovlp->get_ncRNA($chromosome);
@list = $ovlp->get_scRNA($chromosome);
@list = $ovlp->get_piRNA($chromosome);
@list = $ovlp->get_asRNA($chromosome);
@list = $ovlp->get_lincRNA($chromosome);
@list = $ovlp->get_snRNA($chromosome);
@list = $ovlp->get_snoRNA($chromosome);
@list = $ovlp->get_stRNA($chromosome);
@list = $ovlp->get_tRNA($chromosome);
@list = $ovlp->get_tRNAscan_SE_1_23RNA($chromosome);
@list = $ovlp->get_CDS_introns($chromosome);
@list = $ovlp->get_check_introns_EST($chromosome);
@list = $ovlp->get_check_introns_cDNA($chromosome);
@list = $ovlp->get_mass_spec_peptides($chromosome);

Routines to modify/filter GFF lists

@reduced_list = $ovlp->get_span(@list)
@only_complex_repeats = $ovlp->get_repeatmasked_complex(@repeatmasked);



Routine for comparing two lists of GFF data

$compare_obj = $ovlp->compare(\@secondary_list, same_sense => 0, other_sense => 0, near_5 => 0, near_3 => 0)
    all the hash parameters default to False

@secondary_entries = match($main_list_entry_ref)



Routines to return information about the latest match

@ids               = $compare_obj->matching_IDs
@secondary_entries = $compare_obj->matching_data
@sense             = $compare_obj->matching_sense
($prop1, $prop2)   = $compare_obj->matching_proportions($secondary_entry)






=head1 DESCRIPTION

    This does an efficient match for overlapping positions between
    two sets of data.

=head1 CONTACT


Gary gw3@sanger.ac.uk

=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are preceded with a _

=cut




package Overlap;

use strict;
use lib $ENV{CVS_DIR};
use Carp;



=head2

    Title   :   new
    Usage   :   my $ovlp = Overlap->new($database, $wormbase);
    Function:   initialises the data for the Overlap object
    Returns :   Overlap object
    Args    :   wormbase object, database path


=cut

sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;

  $self->{database} = shift or croak "no database supplied\n";
  $self->{wormbase} = shift or croak "no wormbase supplied\n";

  return $self;
}

=head2

    Title   :   wormbase
    Usage   :   my $wormbase = Overlap->wormbase;
    Function:   accessor for the wormbase attribute
    Returns :   Wormbase object
    Args    :   none

=cut

sub wormbase {
	my $self=shift;
	return $self->{wormbase};
}

=head2

    Title   :   compare
    Usage   :   my $ovlp = overlap_obj->compare();
    Function:   initialises the data for comparing two lists of GFF data
    Returns :   GFF comparison object;
    Args    :   1) array ref of list of secondary positions to search against
    Args    :	2) optional hash of search modifiers:
                   near => distance - set distance within which a near match is considered to be a match
                   near_5 => distance - ditto for 5' end of our objects
                   near_3 => distance - ditto for 3' end of our objects
                   same_sense => boolean - say it must have both in the same sense for a match
                   other_sense => boolean - set to true if want them both in the opposite sense for a match
    		   exact_match => boolean - set to true if want them to start and end at the same position
=cut

sub compare {

  # make a copy of the Overlap object for this comparison
  my $caller = shift;
  my $class = ref $caller;
  my $self = bless({}, $class);
  $self->{database} = $caller->{database};
  $self->{wormbase} = $caller->{wormbase};

  # add the data for this comparison
  $self->{secondary_list} = shift or croak "no secondary list of positions to search\n";
  $self->{state} = $self->_init_match(@_);

  return $self;
}


=head2

    Title   :   _init_match
    Usage   :   my %state = $self->_init_match(param => arg, ...);
    Function:   internal method that resets the match state before a search
    Returns :   hash of the initialised match state
    Args    :   near => distance - set distance within which a near match is considered to be a match
                near_5 => distance - ditto for 5' end of our objects
                near_3 => distance - ditto for 3' end of our objects
                same_sense => boolean - say it must have both in the same sense for a match
                other_sense => boolean - set to true if want the both in the opposite sense for a match
                exact_match => boolean - set to true if want them to start and end at the same position
=cut

sub _init_match {
  my $self = shift;
  my (@params) = @_;

  my %state = (last_used  => 0, # reset last secondary list's line 
	       last_used_forward => 0, # reset last secondary list's line for the forward sense
	       last_used_reverse => 0, # reset last secondary list's line for the reverse sense
	       near_5     => 0, # assume we don't allow near 5' matches to count as a match
	       near_3     => 0,	# ditto for 3'
	       same_sense => 0,	# assume we want to allow a match to an object in either sense, set to true if want the same sense only
	       other_sense => 0, # assume we want to allow a match to an object in either sense, set to true if want the opposite sense only
	       exact_match => 0, # assume we do not want an exact match, but some sort of sloppy overlap
	       );

  while (@params) {
    my $param = shift @params;
    if ($param eq "near") {
      my $near = shift @params;
      $state{near_5} = $near;
      $state{near_3} = $near;
    } elsif ($param eq "near_5") {
      my $near_5 = shift @params;
      $state{near_5} = $near_5;
    } elsif ($param eq "near_3") {
      my $near_3 = shift @params;
      $state{near_3} = $near_3;
    } elsif ($param eq "same_sense") {
      my $same_sense = shift @params;
      $state{same_sense} = $same_sense;
    } elsif ($param eq "other_sense") {
      my $other_sense = shift @params;
      $state{other_sense} = $other_sense;
    } elsif ($param eq "exact_match") {
      my $exact_match = shift @params;
      $state{exact_match} = $exact_match;
    }
  }

  # sanity check
  if ($state{same_sense} == 1 && $state{other_sense} == 1) {
    croak ("You can't choose for a match to only things on the same sense and only things on the opposite sense!\n");
  }
     
  return \%state;    
}

=head2

    Title   :   read_GFF_file
    Usage   :   read in the GFF file
    Function:   return $self->read_GFF_file(\%GFF_data);
    Returns :   list of lists for GFF data sorted by chromosomal start
    Args    :   hash containing:
                  directory	        => directory
                  file			=> file name
                  gff_source		=> GFF source name to get
                  gff_type		=> GFF type name to get
                  ID_after		=> regular expression after which the ID is expected to be
                  reverse_orientation   => 1, optional - specifies that the orientation should be reversed if this is a 3-prime read
                  homology              => 1, optional - specifies that the $hit_start, $hit_end, $score data from a homology alignment should also be read
                  other_data            => 1, optional - specifies that the 'other' data field should be appended to the result
                  score                 => 1, optional - specifies the GFF score column should be read and stored in element 5

=cut


sub read_GFF_file {
  my $self = shift;
  my $chromosome = shift;
  my ($GFF_data) = @_;

  my @result;			# returned result

  # the GFF data can be stored in one of several different places.

  # species with small numbers of chromosomes (elegans, briggsae) can
  # have the data in GFF_SPLITS/$chromosome_$method.gff or it can be
  # in CHROMOSOMES/$chromosome.gff if the GFF_SPLIT file has not been
  # made

  # species with large numbers of chromosomes (remanei, brugia) can
  # have the data in GFF_SPLITS/$method.gffor it can be in
  # CHROMOSOMES/$species.gff if the GFF_SPLIT file has not been made

  my @files;
  # if a specific file is given, then use it
  if ($GFF_data->{file}) {
    @files = glob($GFF_data->{file});
    if (!@files || ! -e $files[0]) {die "can't find GFF files $GFF_data->{file}\n";}

  # else look for a GFF_SPLITS file
  } else {
    my $gff_dir = $self->wormbase->gff;
    my $splits_dir = $self->wormbase->gff_splits;
    my $meth = $GFF_data->{method};

    my $single_chr_single_meth_file = "$splits_dir/${chromosome}_${meth}.gff";
    my $single_chr_collated_meth_file = "$gff_dir/${chromosome}.gff";
    my $collated_chr_single_meth_file = "$splits_dir/${meth}.gff";
    my $collated_chr_collated_meth_file = "$gff_dir/" . $self->wormbase->species . ".gff";

    if (-e $single_chr_single_meth_file) {
      @files = ($single_chr_single_meth_file);
    } elsif (-e $single_chr_collated_meth_file) {
      @files = ($single_chr_collated_meth_file);
    } elsif (-e $collated_chr_single_meth_file) {
      @files = ($collated_chr_single_meth_file);
    } elsif (-e $collated_chr_collated_meth_file) {
      @files = ($collated_chr_collated_meth_file);
    } else {
      die "Could not find suitable GFF file(s) for $chromosome $meth\n";
    }
  }

  my $score = 1.0;		# default score
  my ($id, $hit_start, $hit_end);

  foreach my $file (@files) {

    # see if we have to read in the files from disk, or can we reuse a previous read?
    my $key = "${chromosome}_${file}";
    if (! exists $self->{files}{$key}) {
      open (GFF, "grep \"^$chromosome\\W\" $file|") || die "Can't open $file\n";
      my @file = <GFF>;
      close (GFF);
      $self->{files}{$key} = \@file;
    }

    foreach my $line (@{$self->{files}{$key}}) {
      chomp $line;
      if ($line =~ /^\s*$/) {next;}
      if ($line =~ /^#/) {next;}  
      my @f = split /\t/, $line;
      my ($chr, $source, $type, $start, $end, $sense) = ($f[0], $f[1], $f[2], $f[3], $f[4], $f[6]);

#      if ($chr ne $chromosome) {next;}	  

      if ($GFF_data->{gff_source} ne "") {
	if ($GFF_data->{gff_source} =~ /\|/) {
	  if ($source !~ /$GFF_data->{gff_source}/) {next;}
	} else {
	  if ($source ne $GFF_data->{gff_source}) {next;} # want an exact match if there are no '|'
	}
      }
      if ($GFF_data->{gff_type} ne "") {
	if ($GFF_data->{gff_type} =~ /\|/) {
	  if ($type !~ /$GFF_data->{gff_type}/) {next;}
	} else {
	  if ($type ne $GFF_data->{gff_type}) {next;} # want an exact match if there are no '|'
	}
      }

      if (exists $GFF_data->{score}) {	# store the score whether or not we are doing homology stuff
	if ($f[5] =~ /\d+/) {$score = $f[5]} # if the score is numeric, store it
      }

      if (exists $GFF_data->{homology}) {	# do we need to store the homology data?
	($id, $hit_start, $hit_end) = ($f[8] =~ /$GFF_data->{ID_after}(\S+)\s+(\d+)\s+(\d+)/);
	if ($f[5] =~ /\d+/) {$score = $f[5]} # if the score is numeric, store it
      } else {
	($id) = (exists $GFF_data->{ID_after}) ? 
	    ($f[8] =~ /$GFF_data->{ID_after}(\S+)/) : 
	    ($f[8]);	# if no 'ID_after' then the ID is just the Other field
      }
#      print "$line\n";  
      if (! defined $id && exists $GFF_data->{ID_after}) {print "ID NOT DEFINED\nwhen reading $line\nin file $file\n";next;}
      $id =~ s/\"//g;	# remove quotes
      $id =~ s/\;\S+//;	# remove anything after a ';' as this is s field separator in the 'other' field
      my $other_data;
      if (exists $GFF_data->{'other_data'}) {$other_data = $f[8]}
      if (exists $GFF_data->{homology}) {	# do we need to store the homology data?
	push @result, [$id, $start, $end, $sense, $hit_start, $hit_end, $score, $other_data];
      } elsif (exists $GFF_data->{score}) {
	push @result, [$id, $start, $end, $sense, $other_data, $score];	
      } else {
	push @result, [$id, $start, $end, $sense, $other_data];	
      }
    }
#    close (GFF);
  }
  if (scalar @result == 0) {
    print "WARNING: no data found for ";print "@files"; print " looking for $GFF_data->{gff_source}";print " $GFF_data->{gff_type}\n";
  }
    
      
  # The 3' reads of the ESTs/OSTs/mRNAs BLAT results are really in the reverse sense to
  # the gene they match but ACeDB has a tag to display them in the
  # opposite sense as this looks much better.
  #
  # The GFF file holds the original sense though, so we do an explicit
  # flip of the sense of any reversed 3' reads so that when we come to test
  # the ESTs we get the sense coming up the
  # right way round.
  #
 
  if ($GFF_data->{reverse_orientation} && -e $self->{wormbase}->common_data."estorientation.dat") {
    my @ests = @result;
    @result = ();

    my %Show_in_reverse_orientation = $self->{wormbase}->FetchData('estorientation');

    foreach my $est (@ests) {
      # see if we need to reverse the orientation
      # the tag Show_in_reverse_orientation is set if the tag EST_3 is set
      if (exists $Show_in_reverse_orientation{$est->[0]} &&
	  $Show_in_reverse_orientation{$est->[0]} eq '3') {  
	if ($est->[3] eq '+') {
	  $est->[3] = '-';
	} elsif ($est->[3] eq '-') {
	  $est->[3] = '+';
	}
      }
      push @result, $est;
    }
  }

  # return saved result sorted by chromosomal start and end position
  return sort {$a->[1] <=> $b->[1] or $a->[2] <=> $b->[2]} @result;

}

=head2

    Title   :   get_span
    Usage   :   my @est_span = $ovlp->get_span(@est);
    Function:   get chromosomal position span from start to finish of all of HSPs from any one sequence ID
    Returns :   list of lists for GFF data of full span of an ID from start to end
    Args    :   list of lists for GFF data of HSP hits

=cut

sub get_span {
  my $self = shift;
  my (@input) = @_;

  my @out;

  # sort by $id then start,end position
  my @input_sorted = sort {$a->[0] cmp $b->[0]
			     or
			       $a->[1] <=> $b->[1]
				 or 
				   $a->[2] <=> $b->[2]
				 } @input;

  my $prev_id = "";
  my ($id, $start, $end, $sense, $hit_start, $hit_end, $score);
  foreach my $in (@input_sorted) {
    if ($in->[0] ne $prev_id) {	# new ID
      if ($prev_id ne "") {
	# output the previous ID
	push @out, [$id, $start, $end, $sense, $hit_start, $hit_end, $score];
      }
      ($id, $start, $end, $sense, $hit_start, $hit_end, $score) = @{$in};
      $prev_id = $id;
    } else {			# same as the previous ID
      if ($in->[2] > $end) {$end = $in->[2];} # chromosomal end
      if ($in->[6] > $score) {$score = $in->[6];} # score
      if ($in->[4] < $in->[5]) {	# forward sense hit
	if ($in->[5] > $hit_end) {$hit_end = $in->[5];} # hit end
      } else {			# reverse sense hit
	if ($in->[5] < $hit_end) {$hit_end = $in->[5];} # hit end
      }
    }
  }

  # output the last ID
  if ($prev_id ne "") {
    push @out, [$id, $start, $end, $sense, $hit_start, $hit_end, $score];
  }

  # return result sorted by chromosomal start,end position
  return sort {$a->[1] <=> $b->[1] or $a->[2] <=> $b->[2]} @out;
}

=head2

    Title   :   get_paired_span
    Usage   :   my @est_paired_span = $ovlp->get_paired_span(@est);
    Function:   get chromosomal position span from start to finish of all of HSPs including from any pairs of sequences
    Returns :   list of lists for GFF data of full span of IDs including paired reads from start to end
    Args    :   list of lists for GFF data of HSP hits

=cut

sub get_paired_span {
  my $self = shift;
  my (@input) = @_;

  my %pairs;

  # load paired read info
  my $pairs = $self->{database} . "/EST_pairs.txt";    
  open ( PAIRS, "<$pairs") or die("Cant open $pairs :\t$!\n");
  while ( <PAIRS> ) {
    chomp;
    s/\"//g;
    s/Sequence://g;
    next if( ( $_ =~ /acedb/) or ($_ =~ /\/\//) or $_ =~ /^\s*$/);
    my @data = split;
    $pairs{$data[0]} =  $data[1];
    $pairs{$data[1]} =  $data[0];
  }
  close PAIRS;

  my %store;
  # go through ests
  # if the pair of this est has its start,end stored then expand the stored start,end
  # else if this est has its start,end stored then expand the stored start,end
  # else store this est
  foreach my $est (@input) {
    my $id = $est->[0];
    my $pair = $pairs{$id};
    if (defined $pair && exists $store{$pair}) {
      my $start = $store{$pair}->[1];
      my $end = $store{$pair}->[2];
      if ($est->[1] > $start) {$store{$pair}->[1] = $est->[1];} # get greatest span
      if ($est->[2] < $end) {$store{$pair}->[2] = $est->[2];}
      
    } elsif (exists $store{$id}) {
      my $start = $store{$id}->[1];
      my $end = $store{$id}->[2];
      if ($est->[1] > $start) {$store{$id}->[1] = $est->[1];} # get greatest span
      if ($est->[2] < $end) {$store{$id}->[2] = $est->[2];}

    } else {
      $store{$id} = $est;
    }
  }

  my @out = values %store;

  # return result sorted by chromosomal start,end position
  return sort {$a->[1] <=> $b->[1] or $a->[2] <=> $b->[2]} @out;

}


=head2

    Title   :   get_intron_from_exons
    Usage   :   my @introns = $ovlp->get_introns_from_exons(@exons);
    Function:   get a list of introns (the gaps) from between a list of exons (things with gaps between them)
    Returns :   list of lists for GFF data of the gaps (intron) between the things (exons)
    Args    :   list of lists for GFF data of introns

=cut

sub get_intron_from_exons {
  my $self = shift;
  my (@input) = @_;
  my @out;

  # sort by $id then start,end position
  my @input_sorted = sort {$a->[0] cmp $b->[0]
			     or
			       $a->[1] <=> $b->[1]
				 or 
				   $a->[2] <=> $b->[2]
				 } @input;

  my $prev_id = "";
  my $prev_end;
  my $prev_hit_end;
  my ($id, $start, $end, $sense, $hit_start, $hit_end, $score);

  foreach my $in (@input_sorted) {
    ($id, $start, $end, $sense, $hit_start, $hit_end, $score) = @{$in};
    if ($in->[0] eq $prev_id) {	# just gone past an intron
      if (@{$in} > 5) {
	# store the hit positions if known
	if (($hit_start < $hit_end && $prev_hit_end <= $hit_start) || # 5' read
	    ($hit_start > $hit_end && $prev_hit_end >= $hit_start) || # 3' read
	    ($id =~ /Nanopore/) # Nanopore sequences do not have orientation sorted out yet -- just accept anything
	   ) { 
	  push @out, [$id, $prev_end+1, $start-1, $sense, $prev_hit_end, $hit_start, $score];
	}
      } else {
	# otherwise just store the chromosome positions
	push @out, [$id, $prev_end+1, $start-1, $sense];
      }
    }
    $prev_id = $id;
    $prev_end = $end;
    $prev_hit_end = $hit_end if (defined $hit_end);
  }

  # return result sorted by chromosomal start,end position
  return sort {$a->[1] <=> $b->[1] or $a->[2] <=> $b->[2]} @out;

}



=head2

    Title   :   get_EST_BEST
    Usage   :   my @gff = $ovlp->get_EST_BEST($chromosome)
    Function:   reads the GFF data for EST BEST
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_EST_BEST {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "BLAT_EST_BEST",
     gff_source			=> "BLAT_EST_BEST",
     gff_type			=> "expressed_sequence_match",
     ID_after			=> "Target\\s+\"Sequence:",
     reverse_orientation        => 1,
     homology                   => 1,
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_EST_NEMATODE
    Usage   :   my @gff = $ovlp->get_EST_NEMATODE($chromosome)
    Function:   reads the GFF data for EST NEMATODE
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_EST_NEMATODE {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "BLAT_NEMATODE",
     gff_source			=> "EMBL_nematode_cDNAs-BLAT",
     gff_type			=> "expressed_sequence_match",
     ID_after			=> "Target\\s+\"Sequence:",
     reverse_orientation        => 1,
     homology                   => 1,
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_EST_NEMBASE
    Usage   :   my @gff = $ovlp->get_EST_NEMBASE($chromosome)
    Function:   reads the GFF data for EST NEMBASE
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_EST_NEMBASE {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "BLAT_NEMBASE",
     gff_source			=> "NEMBASE_cDNAs-BLAT",
     gff_type			=> "expressed_sequence_match",
     ID_after			=> "Target\\s+\"Sequence:",
     reverse_orientation        => 1,
     homology                   => 1,
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_EST_WASHU
    Usage   :   my @gff = $ovlp->get_EST_WASHU($chromosome)
    Function:   reads the GFF data for EST WASHU
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_EST_WASHU {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "BLAT_WASHU",
     gff_source			=> "NEMATODE.NET_cDNAs-BLAT",
     gff_type			=> "expressed_sequence_match",
     ID_after			=> "Target\\s+\"Sequence:",
     reverse_orientation        => 1,
     homology                   => 1,
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_OST_BEST
    Usage   :   my @gff = $ovlp->get_OST_BEST($chromosome)
    Function:   reads the GFF data for OST BEST
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_OST_BEST {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "BLAT_OST_BEST",
     gff_source			=> "BLAT_OST_BEST",
     gff_type			=> "expressed_sequence_match",
     ID_after			=> "Target\\s+\"Sequence:",
     reverse_orientation        => 1,
     homology                   => 1,

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_mRNA_BEST
    Usage   :   my @gff = $ovlp->get_mRNA_BEST($chromosome)
    Function:   reads the GFF data for mRNA BEST
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_mRNA_BEST {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "BLAT_mRNA_BEST",
     gff_source			=> "BLAT_mRNA_BEST",
     gff_type			=> "expressed_sequence_match",
     ID_after			=> "Target\\s+\"Sequence:",
     reverse_orientation        => 1,
     homology                   => 1,

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_Trinity_BEST
    Usage   :   my @gff = $ovlp->get_Trinity_BEST($chromosome)
    Function:   reads the GFF data for Trinity BEST
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_Trinity_BEST {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "BLAT_Trinity_BEST",
     gff_source			=> "BLAT_Trinity_BEST",
     gff_type			=> "expressed_sequence_match",
     ID_after			=> "Target\\s+\"Sequence:",
     reverse_orientation        => 1,
     homology                   => 1,

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_IsoSeq_BEST
    Usage   :   my @gff = $ovlp->get_IsoSeq_BEST($chromosome)
    Function:   reads the GFF data for IsoSeq BEST
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_IsoSeq_BEST {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> 'BLAT_IsoSeq_BEST',
     gff_source			=> 'BLAT_IsoSeq_BEST',
     gff_type			=> 'expressed_sequence_match',
     ID_after			=> "Target\\s+\"Sequence:",
     reverse_orientation        => 1,
     homology                   => 1,

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);
}

=head2

    Title   :   get_Nanopore_BEST
    Usage   :   my @gff = $ovlp->get_Nanopore_BEST($chromosome)
    Function:   reads the GFF data for Nanopore BEST
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_Nanopore_BEST {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> 'BLAT_Nanopore_BEST',
     gff_source			=> 'BLAT_Nanopore_BEST',
     gff_type			=> 'expressed_sequence_match',
     ID_after			=> "Target\\s+\"Sequence:",
     reverse_orientation        => 1,
     homology                   => 1,

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);
}

=head2

    Title   :   get_RST_BEST
    Usage   :   my @gff = $ovlp->get_RST_BEST($chromosome)
    Function:   reads the GFF data for RST BEST
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_RST_BEST {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "BLAT_RST_BEST",
     gff_source			=> "BLAT_RST_BEST",
     gff_type			=> "expressed_sequence_match",
     ID_after			=> "Target\\s+\"Sequence:",
     reverse_orientation        => 1,
     homology                   => 1,

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_ncRNA_BEST
    Usage   :   my @gff = $ovlp->get_ncRNA_BEST($chromosome)
    Function:   reads the GFF data for ncRNA BEST
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_ncRNA_BEST {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "BLAT_ncRNA_BEST",
     gff_source			=> "BLAT_ncRNA_BEST",
     gff_type			=> "expressed_sequence_match",
     ID_after			=> "Target\\s+\"Sequence:",
     reverse_orientation        => 1,
     homology                   => 1,

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_curated_CDS
    Usage   :   my @gff = $ovlp->get_curated_CDS($chromosome)
    Function:   reads the GFF data for curated CDS (complete CDS from start to end)
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_curated_CDS {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "curated",
     gff_source			=> "curated",
     gff_type			=> "CDS",
     ID_after			=> "CDS\\s+",

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_curated_CDS_exons
    Usage   :   my @gff = $ovlp->get_curated_CDS_exons($chromosome)
    Function:   reads the GFF data for curated_CDS exons
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_curated_CDS_exons {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "curated",
     gff_source			=> "curated",
     gff_type			=> "exon",
     ID_after			=> "CDS\\s+",

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_curated_CDS_introns
    Usage   :   my @gff = $ovlp->get_curated_CDS_introns($chromosome)
    Function:   reads the GFF data for curated_CDS introns
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_curated_CDS_introns {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "curated",
     gff_source			=> "curated",
     gff_type			=> "intron",
     ID_after			=> "CDS\\s+",

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_Coding_transcripts
    Usage   :   my @gff = $ovlp->get_Coding_transcripts($chromosome)
    Function:   reads the GFF data for Coding_transcripts (complete transcript from start to end)
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_Coding_transcripts {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "Coding_transcript",
     gff_source			=> "Coding_transcript",
     gff_type			=> "protein_coding_primary_transcript",
     ID_after			=> "Transcript\\s+",

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_Coding_transcript_exons
    Usage   :   my @gff = $ovlp->get_Coding_transcript_exons($chromosome)
    Function:   reads the GFF data for Coding_transcript exons
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_Coding_transcript_exons {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "Coding_transcript",
     gff_source			=> "Coding_transcript",
     gff_type			=> "exon",
     ID_after			=> "Transcript\\s+",

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_Coding_transcript_introns
    Usage   :   my @gff = $ovlp->get_Coding_transcript_introns($chromosome)
    Function:   reads the GFF data for Coding_transcript introns
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_Coding_transcript_introns {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "Coding_transcript",
     gff_source			=> "Coding_transcript",
     gff_type			=> "intron",
     ID_after			=> "Transcript\\s+",

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_Non_coding_transcripts
    Usage   :   my @gff = $ovlp->get_Non_coding_transcripts($chromosome)
    Function:   reads the GFF data for Non_coding transcript (complete transcript from start to end)
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_Non_coding_transcripts {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "Non_coding_transcript",
     gff_source			=> "non_coding_transcript",
     gff_type			=> "nc_primary_transcript",
     ID_after			=> "Transcript\\s+",

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_Non_coding_transcript_exons
    Usage   :   my @gff = $ovlp->get_Non_coding_transcript_exons($chromosome)
    Function:   reads the GFF data for Non_coding_transcript exons
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_Non_coding_transcript_exons {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "Non_coding_transcript",
     gff_source			=> "non_coding_transcript",
     gff_type			=> "exon",
     ID_after			=> "Transcript\\s+",

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_Non_coding_transcript_introns
    Usage   :   my @gff = $ovlp->get_Non_coding_transcript_introns($chromosome)
    Function:   reads the GFF data for Non_coding_transcript introns
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_Non_coding_transcript_introns {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "Non_coding_transcript",
     gff_source			=> "non_coding_transcript",
     gff_type			=> "intron",
     ID_after			=> "Transcript\\s+",

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_pseudogene
    Usage   :   my @gff = $ovlp->get_pseudogene($chromosome)
    Function:   reads the GFF data for pseudogene
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_pseudogene {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "Pseudogene",
     gff_source			=> "Pseudogene",
     gff_type			=> "Pseudogene",
     ID_after			=> "Pseudogene\\s+",

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_rRNA_transcripts
    Usage   :   my @gff = $ovlp->get_rRNA_transcipts($chromosome)
    Function:   reads the GFF data for rRNA transcipts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_rRNA_transcripts {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "rRNA",
     gff_source			=> "rRNA",
     gff_type			=> "rRNA",
     ID_after			=> "Transcript\\s+",

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_rRNA_exons
    Usage   :   my @gff = $ovlp->get_rRNA_exons($chromosome)
    Function:   reads the GFF data for rRNA exons (usually the same thing as the transcript)
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_rRNA_exons {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "rRNA",
     gff_source			=> "rRNA",
     gff_type			=> "exon",
     ID_after			=> "Transcript\\s+",

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_genefinder_transcripts
    Usage   :   my @gff = $ovlp->get_genefinder_transcripts($chromosome)
    Function:   reads the GFF data for genefinder transcripts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_genefinder_transcripts {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "Genefinder",
     gff_source			=> "Genefinder",
     gff_type			=> "CDS",
     ID_after			=> "CDS\\s+",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_genefinder_exons
    Usage   :   my @gff = $ovlp->get_genefinder_exons($chromosome)
    Function:   reads the GFF data for genefinder exons
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_genefinder_exons {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "Genefinder",
     gff_source			=> "Genefinder",
     gff_type			=> "coding_exon",
     ID_after			=> "CDS\\s+",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_genelets
    Usage   :   my @gff = $ovlp->get_genelets($chromosome)
    Function:   reads the GFF data for RNASeq genelet transcripts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_genelets {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "genelets",
     gff_source			=> "RNASEQ.Hillier",
     gff_type			=> "CDS",
     ID_after			=> "CDS\\s+",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_RNASeq_splice
    Usage   :   my @gff = $ovlp->get_RNASeq_splice($chromosome)
    Function:   reads the GFF data for RNASeq spliced introns
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_RNASeq_splice {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "RNASeq_splice",
     gff_source			=> "RNASeq_splice",
     gff_type			=> "intron",
     # no ID as it is a Feature_data object
     score                      => 1, # force it to store the score
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_Aggregated_CDS_introns
    Usage   :   my @gff = $ovlp->get_Aggregated_CDS_introns($chromosome)
    Function:   reads the GFF data for the Hilllier Aggregated CDS introns
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_Aggregated_CDS_introns {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "Aggregated_CDS_introns",
     gff_source			=> "RNASEQ.Hillier.Aggregate",
     gff_type			=> "intron",
     ID_after			=> "CDS\\s+",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_Aggregated_CDS_exons
    Usage   :   my @gff = $ovlp->get_Aggregated_CDS_exons($chromosome)
    Function:   reads the GFF data for the Hilllier Aggregated CDS exons
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_Aggregated_CDS_exons {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "Aggregated_CDS_exons",
     gff_source			=> "RNASEQ.Hillier.Aggregate",
     gff_type			=> "exon",
     ID_after			=> "CDS\\s+",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_twinscan_transcripts
    Usage   :   my @gff = $ovlp->get_twinscan_transcripts($chromosome)
    Function:   reads the GFF data for twinscan transcripts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_twinscan_transcripts {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "twinscan",
     gff_source			=> "twinscan",
     gff_type			=> "CDS",
     ID_after			=> "CDS\\s+",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_twinscan_exons
    Usage   :   my @gff = $ovlp->get_twinscan_exons($chromosome)
    Function:   reads the GFF data for twinscan exons
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_twinscan_exons {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "twinscan",
     gff_source			=> "twinscan",
     gff_type			=> "coding_exon",
     ID_after			=> "CDS\\s+",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_mgene_transcripts
    Usage   :   my @gff = $ovlp->get_mgene_transcripts($chromosome)
    Function:   reads the GFF data for mgene transcripts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_mgene_transcripts {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "mGene",
     gff_source			=> "mGene",
     gff_type			=> "CDS",
     ID_after			=> "CDS\\s+",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_mgene_exons
    Usage   :   my @gff = $ovlp->get_mgene_exons($chromosome)
    Function:   reads the GFF data for mgene exons
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_mgene_exons {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "mGene",
     gff_source			=> "mGene",
     gff_type			=> "coding_exon",
     ID_after			=> "CDS\\s+",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_jigsaw_CDS
    Usage   :   my @gff = $ovlp->get_jigsaw_CDS($chromosome)
    Function:   reads the GFF data for jigsaw coding transcripts (START to STOP)
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_jigsaw_CDS {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "jigsaw",
     gff_source			=> "jigsaw",
     gff_type			=> "CDS",
     ID_after			=> "CDS\\s+",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_jigsaw_exons
    Usage   :   my @gff = $ovlp->get_jigsaw_exons($chromosome)
    Function:   reads the GFF data for jigsaw exons
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_jigsaw_exons {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "jigsaw",
     gff_source			=> "jigsaw",
     gff_type			=> "coding_exon",
     ID_after			=> "CDS\\s+",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_modencode_CDS
    Usage   :   my @gff = $ovlp->get_jigsaw_CDS($chromosome)
    Function:   reads the GFF data for jigsaw coding transcripts (START to STOP)
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_modencode_CDS {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "RNASEQ.Hillier.Aggregate",
     gff_source			=> "RNASEQ.Hillier.Aggregate",
     gff_type			=> "CDS",
     ID_after			=> "CDS\\s+",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_modencode_exons
    Usage   :   my @gff = $ovlp->get_jigsaw_exons($chromosome)
    Function:   reads the GFF data for jigsaw exons
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_modencode_exons {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "RNASEQ.Hillier.Aggregate",
     gff_source			=> "RNASEQ.Hillier.Aggregate",
     gff_type			=> "coding_exon",
     ID_after			=> "CDS\\s+",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_genblastg_exons
    Usage   :   my @gff = $ovlp->get_genblastg_exons($chromosome)
    Function:   reads the GFF data for genblastg exons (Jack Chen's homology predictions in Tier II species)
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_genblastg_exons {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "genblastg",
     gff_source			=> "genBlastG",
     gff_type			=> "coding_exon",
     ID_after			=> "CDS\\s+",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_transposons
    Usage   :   my @gff = $ovlp->get_transposons($chromosome)
    Function:   reads the GFF data for transposons
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_transposons {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "Transposon",
     gff_source			=> "Transposon",
     gff_type			=> "transposable_element",
     ID_after			=> "Transposon\\s+",

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_transposon_exons
    Usage   :   my @gff = $ovlp->get_transposon_exons($chromosome)
    Function:   reads the GFF data for transposon exons
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_transposon_exons {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "Transposon_CDS",
     gff_source			=> "Transposon_CDS",
     gff_type			=> "coding_exon",
     ID_after			=> "CDS\\s+",

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_polyA_signal
    Usage   :   my @gff = $ovlp->get_polyA_signal($chromosome)
    Function:   reads the GFF data for polyA signal sequences
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_polyA_signal {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "polyA_signal_sequence",
     gff_source			=> "polyA_signal_sequence",
     gff_type			=> "polyA_signal_sequence",
     ID_after			=> 'Feature\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_polyA_site
    Usage   :   my @gff = $ovlp->get_polyA_site($chromosome)
    Function:   reads the GFF data for polyA site features
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_polyA_site {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "polyA_site",
     gff_source			=> "polyA_site",
     gff_type			=> "polyA_site",
     ID_after			=> 'Feature\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_TSL_SL1
    Usage   :   my @gff = $ovlp->get_TSL_SL1($chromosome)
    Function:   reads the GFF data for TSL SL1
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_TSL_SL1 {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "SL1",
     gff_source			=> "SL1",
     gff_type			=> "SL1_acceptor_site",
     ID_after			=> 'Feature\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_TSL_SL2
    Usage   :   my @gff = $ovlp->get_TSL_SL2($chromosome)
    Function:   reads the GFF data for TSL SL2
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_TSL_SL2 {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "SL2",
     gff_source			=> "SL2",
     gff_type			=> "SL2_acceptor_site",
     ID_after			=> 'Feature\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_SAGE_tags
    Usage   :   my @gff = $ovlp->get_SAGE_tags($chromosome)
    Function:   reads the GFF data for SAGE tags
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_SAGE_tags {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "SAGE",
     gff_source			=> "",
     gff_type			=> "SAGE_tag",
     ID_after			=> "Sequence\\s+",
     other_data			=> 1,           # we want to store the other_data field to get the 'count' 
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);
}

=head2

    Title   :   get_blastx_homologies
    Usage   :   my @gff = $ovlp->get_blastx_homologies($chromosome)
    Function:   reads the GFF data for the protein blastx homologies
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_blastx_homologies {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "wublastx",
     gff_source			=> "(proteins-BLASTX|proteins-BLASTX)", # use '|' and double the source name to force it to do a regex comparison to get a partial match excluding the species name at the start
     gff_type			=> "protein_match",
     homology			=> "1",	# this is a GFF with homology data that we need to store
     ID_after			=> "Target\\s+\"Protein:",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_waba_coding
    Usage   :   my @gff = $ovlp->get_waba_coding($chromosome)
    Function:   reads the GFF data for the WABA coding potential
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_waba_coding {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "waba_coding",
     gff_source			=> "waba_coding",
     gff_type			=> "nucleotide_match",
     homology			=> "1",	# this is a GFF with homology data that we need to store
     ID_after			=> "Target\\s+\"Sequence:",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_repeatmasked
    Usage   :   my @gff = $ovlp->get_repeatmasked($chromosome)
    Function:   reads the GFF data for the RepeatMasked regions
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_repeatmasked {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "RepeatMasker",
     gff_source			=> "RepeatMasker",
     gff_type			=> "repeat_region",
     ID_after			=> "Target\\s+\"Motif:",
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_5_UTRs
    Usage   :   my @gff = $ovlp->get_5_UTRs($chromosome)
    Function:   reads the GFF data for the 5-prime UTRs
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_5_UTRs {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "UTR",
     gff_source			=> "Coding_transcript",
     gff_type			=> "five_prime_UTR",
     ID_after			=> 'Transcript\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_3_UTRs
    Usage   :   my @gff = $ovlp->get_3_UTRs($chromosome)
    Function:   reads the GFF data for the 3-prime UTRs
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_3_UTRs {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "UTR",
     gff_source			=> "Coding_transcript",
     gff_type			=> "three_prime_UTR",
     ID_after			=> 'Transcript\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_Genes
    Usage   :   my @gff = $ovlp->get_Genes($chromosome)
    Function:   reads the GFF data for the Genes
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_Genes {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "gene",
     gff_source			=> "gene",
     gff_type			=> "gene",
     ID_after			=> 'Gene\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_Operons
    Usage   :   my @gff = $ovlp->get_Operons($chromosome)
    Function:   reads the GFF data for the Operons
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_Operons {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "operon",
     gff_source			=> "operon",
     gff_type			=> "operon",
     ID_after			=> 'Operon\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_miRNA_mature
    Usage   :   my @gff = $ovlp->get_miRNA_mature($chromosome)
    Function:   reads the GFF data for the mature miRNA transcipts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_miRNA_mature {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "miRNA",
     gff_source			=> "miRNA_mature",
     gff_type			=> "miRNA",
     ID_after			=> 'Transcript\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_miRNA_primary
    Usage   :   my @gff = $ovlp->get_miRNA_primary($chromosome)
    Function:   reads the GFF data for the primary miRNA transcipts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_miRNA_primary {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "miRNA_primary_transcript",
     gff_source			=> "miRNA_primary_transcript",
     gff_type			=> "miRNA_primary_transcript",
     ID_after			=> 'Transcript\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_ncRNA
    Usage   :   my @gff = $ovlp->get_ncRNA($chromosome)
    Function:   reads the GFF data for the ncRNA transcripts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_ncRNA {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "ncRNA",
     gff_source			=> "ncRNA",
     gff_type			=> "ncRNA",
     ID_after			=> 'Transcript\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_scRNA
    Usage   :   my @gff = $ovlp->get_scRNA($chromosome)
    Function:   reads the GFF data for the scRNA transcripts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_scRNA {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "scRNA",
     gff_source			=> "scRNA",
     gff_type			=> "scRNA",
     ID_after			=> 'Transcript\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_asRNA
    Usage   :   my @gff = $ovlp->get_asRNA($chromosome)
    Function:   reads the GFF data for the asRNA transcripts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_asRNA {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "asRNA",
     gff_source			=> "asRNA",
     gff_type			=> "asRNA",
     ID_after			=> 'Transcript\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_lincRNA
    Usage   :   my @gff = $ovlp->get_lincRNA($chromosome)
    Function:   reads the GFF data for the lincRNA transcripts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_lincRNA {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "lincRNA",
     gff_source			=> "lincRNA",
     gff_type			=> "lincRNA",
     ID_after			=> 'Transcript\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_piRNA
    Usage   :   my @gff = $ovlp->get_piRNA($chromosome)
    Function:   reads the GFF data for the piRNA transcripts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_piRNA {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "piRNA",
     gff_source			=> "piRNA",
     gff_type			=> "piRNA",
     ID_after			=> 'Transcript\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_snRNA
    Usage   :   my @gff = $ovlp->get_snRNA($chromosome)
    Function:   reads the GFF data for the snRNA transcripts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_snRNA {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "snRNA",
     gff_source			=> "snRNA",
     gff_type			=> "snRNA",
     ID_after			=> 'Transcript\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_snoRNA
    Usage   :   my @gff = $ovlp->get_snoRNA($chromosome)
    Function:   reads the GFF data for the snoRNA transcripts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_snoRNA {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "snoRNA",
     gff_source			=> "snoRNA",
     gff_type			=> "snoRNA",
     ID_after			=> 'Transcript\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_stRNA
    Usage   :   my @gff = $ovlp->get_stRNA($chromosome)
    Function:   reads the GFF data for the stRNA stranscripts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_stRNA {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "stRNA",
     gff_source			=> "stRNA",
     gff_type			=> "stRNA",
     ID_after			=> 'Transcript\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_tRNA
    Usage   :   my @gff = $ovlp->get_tRNA($chromosome)
    Function:   reads the GFF data for the tRNA transcripts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_tRNA {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "tRNA",
     gff_source			=> "tRNA",
     gff_type			=> "tRNA",
     ID_after			=> 'Transcript\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_tRNAscan_SE_1_23RNA
    Usage   :   my @gff = $ovlp->get_tRNAscan_SE_1_23RNA($chromosome)
    Function:   reads the GFF data for the tRNAscan_SE_1_23RNA transcripts
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_tRNAscan_SE_1_23RNA {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "tRNAscan-SE-1.23",
     gff_source			=> "tRNAscan-SE-1.23",
     gff_type			=> "tRNA",
     ID_after			=> 'Transcript\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_repeatmasked_complex
    Usage   :   my @gff = $ovlp->get_(@repeat)
    Function:   filters a list of RepeatMasker GFF data to leave only the ones with no ()n of '_rich' in their names.
    Returns :   list of lists for GFF data
    Args    :   lists of lists of repeat GFF data from $self->get_repeatmasked()

=cut

 sub get_repeatmasked_complex {
   my $self = shift;
   my (@repeatmasked) = @_;
   
   # want the ones without ')n' in their ID
   my @r2 = grep { $_->[0] !~ /\)n/ } @repeatmasked; 
   # want the ones without '-rich' or in '_rich' their ID
   my @repeatmasked_complex = grep { $_->[0] !~ /[\-\_]rich/ } @r2; 

   return @repeatmasked_complex;

}


=head2

    Title   :   get_CDS_introns
    Usage   :   my @gff = $ovlp->get_CDS_introns($chromosome)
    Function:   reads the GFF data for the curated CDS introns
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_CDS_introns {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "gene_introns",
     gff_source			=> "curated",
     gff_type			=> "intron",
     ID_after			=> 'CDS\s+',
     other_data                 => 1,
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_ignored_introns
    Usage   :   my @gff = $ovlp->get_ignored_introns($chromosome)
    Function:   reads the GFF data for the confirmed EST/mRNA introns that have been marked UTR/false/inconsistent
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_ignored_introns {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
    method			=> "none",
    gff_source			=> "",   # any source, as long as it has a Confirmed_* intron
    gff_type			=> "intron",
    ID_after			=> '.', # any ID 
    other_data                  => 1,
   );

  my @introns = $self->read_GFF_file($chromosome, \%GFF_data);

  # now sift for only the ignored Confirmed_UTR/false/inconsistent introns
  my @out;

  # get the ones with confirmed introns
  my ($id, $start, $end, $sense, $hit_start, $hit_end, $score);
  foreach my $in (@introns) {
    if ($in->[4] =~ /Confirmed_UTR/ || $in->[4] =~ /Confirmed_false/ || $in->[4] =~ /Confirmed_inconsistent/) {
      push @out, $in;
    }
  }

  # return result sorted by chromosomal start,end position
  return sort {$a->[1] <=> $b->[1] or $a->[2] <=> $b->[2]} @out;
}


=head2

    Title   :   get_check_introns_EST
    Usage   :   my @gff = $ovlp->get_check_introns_EST($chromosome)
    Function:   reads the GFF data for the CHECKS/ confirmed EST introns that dont match the curated introns
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_check_introns_EST {
  my $self = shift;
  my ($chromosome) = @_;

  my $file;
  $file = $self->wormbase->autoace . "/CHECKS/${chromosome}.check_intron_*.gff";
  my @files = glob($file);
  if (!@files || ! -e $files[0]) {
    $file = "/nfs/WWWdev/SANGER_docs/htdocs/Projects/C_elegans/WORMBASE/development_release/GFF/${chromosome}.check_intron_*.gff"; # web director
  }


  my %GFF_data = 
   (
     file			=> "$file",
     gff_source			=> "",
     gff_type			=> "intron",
     ID_after			=> 'Confirmed_EST\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


=head2

    Title   :   get_check_introns_cDNA
    Usage   :   my @gff = $ovlp->get_check_introns_cDNA($chromosome)
    Function:   reads the GFF data for the the confirmed cDNA introns that dont match the curated introns
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_check_introns_cDNA {
  my $self = shift;
  my ($chromosome) = @_;

  my $file;
  $file = $self->wormbase->autoace . "/CHECKS/${chromosome}.check_intron_*.gff";
  my @files = glob($file);
  if (!@files || ! -e $files[0]) {
    $file = "/nfs/WWWdev/SANGER_docs/htdocs/Projects/C_elegans/WORMBASE/development_release/GFF/${chromosome}.check_intron_*.gff"; # web director
  }


  my %GFF_data = 
   (
    file			=> "$file",
     gff_source			=> "",
     gff_type			=> "intron",
     ID_after			=> 'Confirmed_cDNA\s+',

   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}

=head2

    Title   :   get_mass_spec_peptides
    Usage   :   my @gff = $ovlp->get_mass_spec_peptides($chromosome)
    Function:   reads the GFF data for the the mass spectromtry peptides aligned to the genome
    Returns :   list of lists for GFF data
    Args    :   chromosome number

=cut

sub get_mass_spec_peptides {
  my $self = shift;
  my ($chromosome) = @_;

  my %GFF_data = 
   (
     method			=> "mass_spec_genome",
     gff_source			=> "mass_spec_genome",
     gff_type			=> "translated_nucleotide_match",
     ID_after			=> 'Target\s+',
     homology                   => 1,
   );

  return $self->read_GFF_file($chromosome, \%GFF_data);

}


############################################################################

=head2

    Title   :   match
    Usage   :   @result = $overlap_obj->match($this_line);
    Function:   routine to do generalised matching of the region in this_line versus all of the regions in @secondary_list
    Returns :   an array of the secondary_list entries matched 
                the status hash also holds the array of the secondary_list entries matched in 'matching_data'
    Args    :   $this_line - ref to array holding ID, chrom_start, chrom_end, chrom_strand and possibly other values on the end
                $secondary_list - ref to array of arrays sorted by start position - each of which holds ID, chrom_start, chrom_end, chrom_strand and possibly other values on the end

=cut

sub match {
  my $self = shift;

  my ($main) = @_;

  my $no_of_matches = 0;		# result of the match - no of matches found

  my $main_start = $main->[1];
  my $main_end = $main->[2];
  my $main_strand = $main->[3];
  #print "THIS LINE: $main->[0] $main_start $main_end $main_strand\n";

  my $near_3;
  my $near_5;
  if ($main_strand eq '+') {
    $near_3 = $self->{state}->{near_3};
    $near_5 = $self->{state}->{near_5};
  } else {			# swap the values around in the reverse sense
    $near_3 = $self->{state}->{near_5};
    $near_5 = $self->{state}->{near_3};
  }

  @{$self->{state}->{matching_data}} = (); # no matching IDs found yet

  # if searching for same/opposite sense matches, then set last_used to be the
  # minimum of last_used_forward and last_used_reverse
  if ($self->{state}->{same_sense} || $self->{state}->{other_sense}) {
    $self->{state}->{last_used} = $self->{state}->{last_used_forward}; 
    if ($self->{state}->{last_used_reverse} < 
	$self->{state}->{last_used_forward}) {
      $self->{state}->{last_used} = $self->{state}->{last_used_reverse};
    }
  }

  #print "last_used $self->{state}->{last_used}\n";

  for (my $i = $self->{state}->{last_used}; defined $self->{secondary_list}->[$i]; $i++) {

    my $secondary = $self->{secondary_list}->[$i];

    # secondary_start = $secondary->[1]
    # secondary_end = $secondary->[2]

    # test if there is an overlap (allowing possible nearby matches)
    # and optionally test if the senses are the same
    #print "SECONDARY LINE: $secondary->[0] $secondary->[1] $secondary->[2] $secondary->[3]\n";
    if (
	(
	 !$self->{state}->{exact_match} &&
         $secondary->[1] <= $main_end + $near_5 && 
	 $secondary->[2] >= $main_start - $near_3 && 
	 ($self->{state}->{same_sense}?($main_strand eq $secondary->[3]):1) &&
	 ($self->{state}->{other_sense}?($main_strand ne $secondary->[3]):1)
	 ) ||

	# for when we require an exact match
	(
	 $self->{state}->{exact_match} &&
	 $secondary->[1] == $main_start &&
	 $secondary->[2] == $main_end &&
	 ($self->{state}->{same_sense}?($main_strand eq $secondary->[3]):1) &&
	 ($self->{state}->{other_sense}?($main_strand ne $secondary->[3]):1)
	 )
	) {

      # note that we have a match
      $no_of_matches++;

      # if we have not yet noted a match for this line, then remember
      # we got to here and found the first match
      if ($no_of_matches == 1) {
	# see if we are testing for them to be in the same/opposite sense
	if ($self->{state}->{same_sense} || $self->{state}->{other_sense}) {
	  if ($main_strand eq '+') {
	    # remember where we got up to in this sense
	    $self->{state}->{last_used_forward} = $i;
	  } else {
	    $self->{state}->{last_used_reverse} = $i;
	  }
	}
	$self->{state}->{last_used} = $i;
      }

    
      # save information about this match
      $self->{state}->{main_entry} = $main;
      push @{$self->{state}->{matching_data}}, $secondary; # the secondary_list entry 

    } else {
      #print "no match\n";
      # don't search any further for this one if no overlap and the secondary_start > this_end
      if ($secondary->[1] > $main_end + $near_5) {last;}	# we have gone past the end of all possible matches
    }

  }
  #print "out of SECONDARY loop\n";

  # return the matching secondary_list entries
  # this is a count of the matches found when in scalar mode
  return @{$self->{state}->{matching_data}}; 

}

=head2

    Title   :   matching_IDs
    Usage   :   my @IDs = $overlap_obj->matching_IDs;
    Function:   return array of the IDs of the matches of the latest comparison
    Returns :   array of the IDs
    Args    :   none

=cut

sub matching_IDs {
  my $self = shift;

  my @ids;
  foreach my $secondary (@{$self->{state}->{matching_data}}) {
    push @ids, $secondary->[0];
  }
  return @ids;
}


=head2

    Title   :   matching_data
    Usage   :   my @gff = $overlap_obj->matching_data;
    Function:   Return array of the secondary list matches. This is also returned by match(), but it may be useful to be able to access it again
    Returns :   list of lists for GFF data
    Args    :   none

=cut

sub matching_data {
  my $self = shift;
  return @{$self->{state}->{matching_data}};
}


=head2

    Title   :   matching_sense
    Usage   :   my @same_sense = $overlap_obj->matching_sense;
    Function:   return array of the logical AND of the senses of the matches, 1 if the senses were the same, 0 if they differ
    Returns :   list of logical AND of the senses of the matches
    Args    :   none

=cut

sub matching_sense {
  my $self = shift;

  my @senses;
  foreach my $secondary (@{$self->{state}->{matching_data}}) {
    # 1 if the senses are the same, 0 if not
    push @senses, ($self->{state}->{main_entry}->[3] eq $secondary->[3]?1:0);
  }

  return @senses;
}

=head2

    Title   :   matching_proportions
    Usage   :   my ($prop1, $prop2) = $overlap_obj->matching_proportions($secondary_entry);
    Function:   return the proportions of the overlap of the main entry and a specified secondary entry
    Returns :   the proportion of the main entry which is covered by the secondary entry
                the proportion of the secondary entry which is covered by the main entry
    Args    :   an array ref to a secondary_list entry

=cut

 sub matching_proportions {
  my $self = shift;
  my ($secondary) = @_;

  my ($prop1, $prop2);

  my $main = $self->{state}->{main_entry};
  my $main_start = $main->[1];
  my $main_end = $main->[2];
  my $secondary_start = $secondary->[1];
  my $secondary_end = $secondary->[2];


  # get the proportion of the main entry that the secondary entry overlaps with
  if ($main_start > $secondary_start && $main_end < $secondary_end)  {
    $prop1 = 1.0;      # 100% overlap
  } elsif ($secondary_end < $main_start || $main_end < $secondary_start) {
    $prop1 = 0.0;
  } else {
    if ($main_start > $secondary_start) { # secondary thing overlaps main at the start
      $prop1 = ($secondary_end - $main_start + 1) / ($main_end - $main_start + 1);
    } else {
      $prop1 = ($main_end - $secondary_start + 1) / ($main_end - $main_start + 1);
    }
  }

  # get the proportion of the secondary entry that the main entry overlaps with
  if ($secondary_start > $main_start && $secondary_end < $main_end)  {
    $prop2 = 1.0;      # 100% overlap
  } elsif ($main_end < $secondary_start || $secondary_end < $main_start) {
    $prop2 = 0.0;
  } else {
    if ($secondary_start > $main_start) { # main thing overlaps secondary at the start
      $prop2 = ($main_end - $secondary_start + 1) / ($secondary_end - $secondary_start + 1);
    } else {
      $prop2 = ($secondary_end - $main_start + 1) / ($secondary_end - $secondary_start + 1);
    }
  }


  return ($prop1, $prop2);

  }




1;

__END__


