#!/usr/local/bin/perl5.8.0 -w
#
# activate.pl
# 
# by Gary Williams                         
#
# This does stuff with what is in the active zone
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2014-03-13 12:40:19 $      




use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
use Sequence_extract;
use Coords_converter;
use Tk;
use Tk::Dialog; 
use Tk::FileDialog;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($species, $database, $notsl);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "database:s" => \$database, # database being curated
	    "notsl"      => \$notsl, # don't try to make TSL isoforms - for debugging purposes
	    );

# always in debug mode
$debug = $ENV{USER};


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

if (! defined $species) {$species = $wormbase->species}

my $USER;
if (! defined $database) {
  $USER = $ENV{USER};
  if ($species eq 'elegans') {
    $database = "/nfs/wormpub/camace_${USER}";
  } else {
    $database = "/nfs/wormpub/${species}_${USER}";
  }
}

##########################
# MAIN BODY OF SCRIPT
##########################

# splices are held as array of splice structures

# splice node hash structure
# int - flag to ignore this splice
# char - seq - chromosome or contig
# int - start
# int - end
# char - sense
# int - score
# int array - list of possible child introns in the structure
# int - pos (current position in list of child nodes)

my $coords = Coords_converter->invoke($database, undef, $wormbase);

mkdir "$database/CHROMOSOMES", 0777; # we may need to write a new set of chromosome files if they do not yet exist
my $seq_obj = Sequence_extract->invoke($database, undef, $wormbase);

#print "Connecting to Ace\n";
my $db = Ace->connect(-path => $database) || die "cannot connect to database at $database\n";

# load the method to display the isoformer CDS structures
&load_isoformer_method();

# get all RNASeq splices starting and ending in the region
# get all SL1, SL2 sites
my %all_splices;
my %all_TSL;
print "Reading splice data\n";
foreach my $chromosome ($wormbase->get_chromosome_names(-mito => 1, -prefix => 1)) {
  if (!exists $all_splices{$chromosome}) {
    my ($splice_data, $TSL_data) = read_splice_and_TSL_data($chromosome, $log);
    $all_splices{$chromosome} = $splice_data;
    $all_TSL{$chromosome} = $TSL_data;
  }
}

if ($notsl) {%all_TSL=()} # don't use TSL data - for debugging purposes

while (1) {

  my ($chromosome, $region_start, $region_end, $sense, $biotype);
  do {
    print "Next region +-start +-end> ";
    my $userinput =  <STDIN>;
    chomp ($userinput);
    print "user input: $userinput\n";
    if (!defined $userinput || $userinput eq '') {$userinput = 'F22E12.4a'} # for debugging

    # get the region of interest from the CDS name or clone positions
    ($chromosome, $region_start, $region_end, $sense, $biotype) = get_active_region($userinput);
  } while (! defined $chromosome);

  my @TSL = &get_TSL_in_region($chromosome, $region_start, $region_end, $sense, $all_TSL{$chromosome});

  # add the $region_start/end to the TSL list as a dummy TSL site as it is one of the positions that enclose the regions to consider
  if ($sense eq '+') {
    push @TSL, {
		seq => $chromosome,
		start => $region_start,
		end => $region_start,
		sense => $sense,
		tsl => '',
	       }
  } else {
    push @TSL, {
		seq => $chromosome,
		start => $region_end,
		end => $region_end,
		sense => $sense,
		tsl => '',
	       }
  }

  my $next_isoform = 1; # used in constructing the unique name of the object
  my @confirmed;
  my @created;

  foreach my $TSL (@TSL) {
    if ($sense eq '+') {
      $region_start = $TSL->{start};
    } else {
      $region_end = $TSL->{start};
    }

    my @splices = &get_splices_in_region($chromosome, $region_start, $region_end, $sense, $all_splices{$chromosome});
    
    # prepend dummy head splice for all other splice structures to hang from
    unshift @splices, {head => 1, seq => $chromosome, start => $region_start-1, end => $region_start-1}; # don't set a score
    
    # ignore overlapped introns with a score <10% of the other intron
    foreach my $splice1 (@splices) {
      if (exists $splice1->{ignore}) {next}
      if (exists $splice1->{head}) {next}
      foreach my $splice2 (@splices) {
	if (exists $splice2->{ignore}) {next}
	if (exists $splice2->{head}) {next}
	#    if the two splices are the same then next
	if ($splice1->{start} == $splice2->{start} && $splice1->{end} == $splice2->{end} && $splice1->{sense} eq $splice2->{sense}) {next}
	# if the two splices overlap and splice1 is < 10% of score of splice2
	if (($splice1->{start} >= $splice2->{start} && $splice1->{start} <= $splice2->{end}) ||
	    ($splice1->{end} >= $splice2->{start} && $splice1->{end} <= $splice2->{end})) {
	  if ($splice1->{score} * 10 < $splice2->{score}) {
	    $splice1->{ignore} = 1;
	  }
	}
      }
    }
    
    # get average of active splices
    my $average_score = 0;
    my $count = 0;
    foreach my $splice (@splices) {
      if (exists $splice->{ignore}) {next}
      if (exists $splice->{head}) {next}
      $average_score += $splice->{score};
      $count++;
    }
    if ($count == 0) {
      $average_score = 0;
    } else {
      $average_score /= $count;
    }
    
    # if splice1 < 10% of average of still active splices
    #   mark splice1 as 'ignore'
    # this removes low-scoring spurious introns inside exons that are not
    # overlapped by other introns
    foreach my $splice (@splices) {
      if (exists $splice->{head}) {next}
      if ($splice->{score} * 10 < $average_score) {$splice->{ignore} = 1}
    }
    
    # discard splices marked as 'ignore'
    my @splices2=();
    foreach my $splice (@splices) {
      if (!exists $splice->{ignore}) {push @splices2, $splice}
    }
    
    # sort splices by start position
    @splices = sort {$a->{start} <=> $b->{start}} @splices2;
    
    # reset all splice nodes lists and 'current position in list of child nodes' values to 0
    foreach my $splice (@splices) {
      $splice->{list} = [];
      $splice->{pos} = 0;
    }
    
    # get the list of possible child introns of each intron
    &determine_child_nodes(@splices);
    
    # iterate through the valid structures
    my ($confirmed, $created) = &iterate_through_the_valid_structures(\$next_isoform, $TSL, $chromosome, $region_start, $region_end, $sense, @splices);
    push @confirmed, @{$confirmed};
    push @created, @{$created};

  } # foreach TSL

  print "\n*** The following structures were confirmed: @confirmed\n\n";
  print "\n*** The following novel structures were created: @created\n\n";

}

# close the ACE connection
$db->close;

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);



##############################################################
#
# Subroutines
#
##############################################################

# load the method to display the isoformer CDS structures

sub load_isoformer_method {

  if (!-e "$database/tmp") {mkdir "$database/tmp", 0777};
  my $output = "$database/tmp/isoformer_method$$";

  open (OUT, "> $output") || die "Can't open $output to write the Method\n";
  print OUT "\n\n";
  print OUT "Method : \"isoformer\"\n";
  print OUT "Colour DARKGREEN\n";
  print OUT "CDS_colour DARKGREEN\n";
  print OUT "Strand_sensitive Show_up_strand\n";
  print OUT "Right_priority   2.63\n";
  print OUT "Overlap_mode Bumpable\n";
  print OUT "Show_text\n";
  print OUT "EMBL_feature CDS\n";
  print OUT "Remark \"This method is used by acedb to display curation CDS structures generated by isoformer.pl\"\n";
  print OUT "\n\n";
  print OUT "\n\n";
  print OUT "Method : \"non_coding_isoformer\"\n";
  print OUT "Colour ORANGE\n";
  print OUT "CDS_colour ORANGE\n";
  print OUT "Strand_sensitive Show_up_strand\n";
  print OUT "Right_priority   2.6\n";
  print OUT "Overlap_mode Bumpable\n";
  print OUT "Show_text\n";
  print OUT "EMBL_feature ncRNA\n";
  print OUT "Remark \"This method is used by acedb to display curation non-coding structures (ncRNA/Pseudogene/non_coding_isoform) generated by isoformer.pl\"\n";
  print OUT "\n\n";
  close(OUT);

  my $return_status = system("xremote -remote 'parse $output'");
  if ( ( $return_status >> 8 ) != 0 ) {
    print STDERR "WARNING - X11 connection appears to be lost\n";
    #      &error_warning("WARNING", "X11 connection appears to be lost");
  } else {
    print "Loaded isoformer Method\n";
  }

}

###############################################################################
# return the chromosomal location and sense of a CDS or a clone:start-end triplet

# get the region of interest from the CDS name or clone positions
# region can be specified as either a CDS name or a 'seq start end' triplet
# with either space or ':' or '-' separators

# reverse sense regions are indicated by having the end < start these
# are reversed when the positions are returned, giving start end in
# the forward sense with a sense character

sub get_active_region {
  my ($region) = @_;
  my ($chromosome, $start, $end, $sense, $biotype);
  $sense = '+';

  my @region = split /[\s]+/, $region;
  if ($region[0]) { # CDS or Pseudogene
    my $obj = $db->fetch(CDS => "$region[0]");
    if (defined $obj) {
      my $clone_start;
      my $clone_end;
      my $clone = $obj->Sequence;
      my @clone_CDSs = $clone->CDS_child;
      foreach my $CDS ( @clone_CDSs ) {
	next unless ($CDS->name eq "$region[0]");
	$clone_start = $CDS->right->name;
	$clone_end = $CDS->right->right->name;
	last;
      }
      ($chromosome, $start) = $coords->Coords_2chrom_coords($clone, $clone_start);
      ($chromosome, $end) = $coords->Coords_2chrom_coords($clone, $clone_end);
      $biotype = 'CDS';

    } else { # try Pseudogene
      my $obj = $db->fetch(Pseudogene => "$region[0]");
      if (! defined $obj) {
	print STDERR "Region: '$region' is an unknown way of specifying a region: not a CDS or a Pseudogene";
	return undef;
      }
      my $clone_start;
      my $clone_end;
      my $clone = $obj->Sequence;
      my @clone_Pseudogene = $clone->Pseudogene;
      foreach my $Pseudogene ( @clone_Pseudogene ) {
	next unless ($Pseudogene->name eq "$region[0]");
	$clone_start = $Pseudogene->right->name;
	$clone_end = $Pseudogene->right->right->name;
	last;
      }
      ($chromosome, $start) = $coords->Coords_2chrom_coords($clone, $clone_start);
      ($chromosome, $end) = $coords->Coords_2chrom_coords($clone, $clone_end);
      $biotype = 'Pseudogene';

    }

    if ($region[1]) {
      $start += $region[1];
    }

    if ($region[2]) {
      $end = $start - $region[1] + $region[2];
    }

  } else {
    print STDERR "Region: '$region' is an unknown way of specifying a region";
    return undef;
  }

  if ($start > $end) {
    ($start, $end) = ($end, $start);
    $sense = '-';
  }

  return ($chromosome, $start, $end, $sense, $biotype);
}
###############################################################################
# read the RNASeq intron data and TSL site data from the GFF file

sub read_splice_and_TSL_data {
  my ($chromosome, $log) = @_;
  my @splice_data;
  my @TSL_data;
  my %splice;
  my $splice_method = 'RNASeq_splice';
  my $SL1_method = 'SL1';
  my $SL2_method = 'SL2';
  my $file = open_GFF_file($chromosome, $log);

  while (<$file>) {
    my @line = split /\t/, $_;
    if ($line[1] eq $splice_method) {
      push @splice_data, {
			  seq => $line[0],
			  start => $line[3],
			  end => $line[4],
			  score => $line[5],
			  sense => $line[6],
			 }
    } elsif ($line[1] eq $SL1_method) {
      my ($id) = ($line[8] =~ /Feature\s+\"(\S+)\"/); # Feature "WBsf161341"
      push @TSL_data, {
		       seq => $line[0],
		       start => $line[3],
		       end => $line[4],
		       score => $line[5],
		       sense => $line[6],
		       tsl => 'SL1',
		       id => $id,
		      }

    } elsif ($line[1] eq $SL2_method) {
      my ($id) = ($line[8] =~ /Feature\s+\"(\S+)\"/); # Feature "WBsf161341"
      push @TSL_data, {
		       seq => $line[0],
		       start => $line[3],
		       end => $line[4],
		       score => $line[5],
		       sense => $line[6],
		       tsl => 'SL2',
		       id => $id,
		      }
    }
  }

  return (\@splice_data, \@TSL_data);
}

###############################################################################
# get the GFF file handle

sub open_GFF_file {
  my ($chromosome, $log) = @_;
  my $file;
  my $fh;

  if ($species eq 'elegans') {
    $file = "/nfs/wormpub/DATABASES/current_DB/CHROMOSOMES/$chromosome.gff";
  } else {
    if ($wormbase->assembly_type eq 'contig') {      
      $file = $wormbase->sequences . "/" . $wormbase->species.".gff";

    } else {
      $file =  "/nfs/wormpub/BUILD/briggsae/CHROMOSOMES/$chromosome.gff";
    }
  }

  open($fh, "cat $file | grep \"^$chromosome\\W\" |") or $log->log_and_die("Could not open seq-restricted stream to $file\n");
  return $fh;

}
###############################################################################
# filter the list of splices and return those enclose by the specied region

sub get_splices_in_region {
  my ($chromosome, $start, $end, $sense, $all_splices) = @_;
  my @splices;

  foreach my $splice (@{$all_splices}) {
    if ($splice->{start} < $start) {next}
    if ($splice->{end} > $end) {next}
    if ($splice->{sense} ne $sense) {next}
    if ($splice->{seq} ne $chromosome) {next}
    push @splices, $splice;
  }

  return @splices;
}
###############################################################################
# filter the list of TSLs and return those enclose by the specied region

sub get_TSL_in_region {
  my ($chromosome, $start, $end, $sense, $all_TSL) = @_;
  my @TSL;

  foreach my $TSL (@{$all_TSL}) {
    if ($TSL->{start} < $start) {next}
    if ($TSL->{end} > $end) {next}
    if ($TSL->{sense} ne $sense) {next}
    if ($TSL->{seq} ne $chromosome) {next}
    push @TSL, $TSL;
  }

  return @TSL;
}
###############################################################################
# get the list of possible child introns of each intron
# a child intron is one which the current intron can reasonably be
# followed by in a typical CDS structure

# define determine child nodes
# foreach splice1
#   foreach slice2
#     if splice2 start > splice1 end
#       if splice1's child node list is empty
#         add splice2 to splice1's child node list
#       else 
#         foreach child node in splice1's child node list
#           if splice2 start < child node end
#             add splice2 to splice1's child node list
# return list of splices populated with their child node lists

sub determine_child_nodes {
  my (@splices) = @_; # list of introns sorted by start pos

  my $min_exon_length=1; # we very rarely get intron less that 25, so this is a nice, safe min length
  my $overlaps_all_other_child_nodes;

  # The splices are sorted by start position and so we know that all
  # subsequence splices must be either overlapping this node or are
  # potential child-nodes of this node, and will either be overlapping
  # or after any already-identified child-nodes of this splice.

  foreach my $splice1 (@splices) {
    foreach my $splice2 (@splices) {
      if ($splice1->{start} == $splice2->{start} && $splice1->{end} == $splice2->{end}) {next} # ignore same splice
      if ($splice2->{start} > $splice1->{end} + $min_exon_length) { # splice 2 is after splice 1 plus an exon's length so can be a child node
	if ($#{$splice1->{list}} == -1) { # if splice1's child node list is empty
	  push @{$splice1->{list}}, $splice2; # add splice2 to splice1's child node list
	  #print "add score $splice2->{score} to score $splice1->{score} list\n";
	} else {
	  $overlaps_all_other_child_nodes = 1;
	  foreach my $splice3 (@{$splice1->{list}}) { # foreach child node in splice1's child node list
	    if ($splice2->{start} > $splice3->{end}) { 
	      $overlaps_all_other_child_nodes = 0; # we don't have an overlap here and so do not have alternate child nodes
	    }
	  }
	  if ($overlaps_all_other_child_nodes) {
	    push @{$splice1->{list}}, $splice2; # add splice2 to splice1's child node list
	    #print "add score $splice2->{score} to score $splice1->{score} list\n";
	  }
	}
      }
    }
  }
  print "All child nodes done\n";
}

###############################################################################

# define iterate through the valid structures
# repeat:
#   run down list of nodes chaining them together
#   note the positions where there is a list of > 1 child node
#   check for premature STOPs
#   check structure is not already curated
#   save this structure as a FMAP isoform_X if it looks good
#   setup for next iteration by taking the last node with multiple child nodes and incrementing the current position
#   if the current position is off the end of the list then reset it to 0 and increment the next level up node with multiple child nodes etc.
#   if no more levels up to increment then return

sub iterate_through_the_valid_structures {
  my ($next_isoform, $TSL, $chromosome, $region_start, $region_end, $sense, @splices) = @_;

  my $finished=0; # set when we have no more levels to increment up

  my @confirmed=();
  my @created=();

  do {

    my ($chain) = &get_next_chain_of_intron_nodes(@splices);

    my ($chrom_aug, $chrom_stop, $type, $exons) = &convert_to_exons($region_start, $region_end, $chain, $sense);
    
    if ($type eq 'CDS') {
      my ($clone, $clone_aug, $clone_stop) = get_clone_coords($chromosome, $chrom_aug, $chrom_stop, $sense);
      my $confirmed = &check_not_already_curated($type, $clone, $clone_aug, $clone_stop, $exons);
      if ($confirmed) {
	push @confirmed, $confirmed;
      } else {
	&make_isoform('isoformer', $TSL, $clone, $clone_aug, $clone_stop, $sense, $exons, $$next_isoform);
	push @created, "isoformer_${$next_isoform}";
	$$next_isoform++;
      }
    } elsif ($type eq 'INVALID') { # no START and STOP found in the first and last exons
      print "Found no good coding structure - should the region to look at be larger?\n";
      my ($clone, $clone_aug, $clone_stop) = get_clone_coords($chromosome, $chrom_aug, $chrom_stop, $sense);
      my $confirmed = &check_not_already_curated($type, $clone, $clone_aug, $clone_stop, $exons);
      if ($confirmed) {
	push @confirmed, $confirmed;
      } else {
	&make_isoform('non_coding_isoformer', $TSL, $clone, $clone_aug, $clone_stop, $sense, $exons, $$next_isoform);
	push @created, "non_coding_isoformer_${$next_isoform}";
	$$next_isoform++;
      }

    } else {
      die "Type '$type' not defined\n\n";
    }

    $finished = &increment_alternate_splicing($chain);

  } until ($finished);

  return (\@confirmed, \@created);
}

###############################################################################

# the array 'list' in the current splice object contains the alternate
# splices that can be added after the current splice

# the integer 'pos' points to the next alternate splice in the 'list'
# of alternate splices to be used

# this routine runs down the chain of splices getting the next splice
# to add to the growing chain to create the next structure

sub get_next_chain_of_intron_nodes {
  my (@splices) = @_;
  my (@chain);

  my $next_intron = $splices[0]; # start here

  my $count = 0;
  while (1) {
    my $start = $next_intron->{start};
    my $end = $next_intron->{end};
    print "$start-$end ";
    push @chain, $next_intron;
    if ($#{$next_intron->{list}} == -1) {last}
    my $pos = $next_intron->{pos};
    $next_intron = $next_intron->{list}[$pos];
    $count++;
  }
  print "\n";

  return (\@chain);
}
###############################################################################
# this takes a chain of splices and converts them into an array of exon hashes
#
# get the sequence
# find the pos and frame of STOP codons in the sequence
# find the longest ORF from a START codon to the next STOP codon in-frame
# convert the chain of splices into an array of exon start-end positions after the initial START pos
# convert the START and STOP positions to chromosome coords
# do a simple classification of the structure depending on whether the structure uses all exons or not
# and whether a START and STOP codon were found.

#
# FORWARD SENSE     
#      
# CDS start   
# exon 1
#  |   
# exon 2
#  |   
# exon 3
# CDS end
# 


#      
# REVERSE SENSE     
#      
# CDS end   
# exon 3
#  |   
# exon 2   
#  |   
# exon 1   
# CDS start
#      
#      


sub convert_to_exons {
  my ($region_start, $region_end, $chain, $sense) = @_;
  my $type;
  my $cds_seq;

  $type = 'CDS'; # assume we have a valid structure

  # get the sequence
  my $prev_end = $region_start;
  foreach my $splice (@{$chain}) {
    if (exists $splice->{head}) {next}
    my $seq = $splice->{seq};
    my $splice_start = $splice->{start};
    my $splice_end = $splice->{end};
    my $len = $splice_start-$prev_end;
    $cds_seq .= $seq_obj->Sub_sequence($seq, $prev_end-1, $len);
    $prev_end = $splice_end+1; # start of next exon
  }
  $cds_seq .= $seq_obj->Sub_sequence($chain->[0]{seq}, $prev_end-1, $region_end-$prev_end+1);
  if ($sense eq '-') {
    $cds_seq = $seq_obj->DNA_revcomp($cds_seq);
  }
  #print "$cds_seq\n";

  # get the STOP codons
  my @stop_frames;
  foreach my $frame (0..2) { # make three empty arrays to ensure they exist
    $stop_frames[$frame] = [];
  }
  foreach my $stop ('taa', 'tga', 'tag') {
    my $st_offset = 0;
    my $st_result = index($cds_seq, $stop, $st_offset);
    while ($st_result != -1) {
      my $frame = $st_result % 3;
      push @{$stop_frames[$frame]}, $st_result+2; # +2 gets the end base of the codon
      $st_offset = $st_result + 1;
      $st_result = index($cds_seq, $stop, $st_offset);
    }
  }
  # and sort each frame by pos of STOP codon
  foreach my $frame (0..2) {
    @{$stop_frames[$frame]} = sort {$a <=> $b} @{$stop_frames[$frame]};
  }

  # find longest ORF
  my $aug = 'atg';
  my $offset = 0;
  my $result = index($cds_seq, $aug, $offset);
  my $best_orf_len = -1;
  my $best_aug = 0;
  my $best_stop = 0;
  while ($result != -1) {
    my $frame = $result % 3;
    foreach my $stop (@{$stop_frames[$frame]}) {
      if ($stop > $result) {
	if ($stop - $result + 1 > $best_orf_len) {
	  $best_orf_len = $stop - $result + 1;
	  $best_aug = $result;
	  $best_stop = $stop;
	}
	last; # don't check any further STOPs for this AUG
      }
    }
    $offset = $result + 1;
    $result = index($cds_seq, $aug, $offset);
  }
  print "best_aug: $best_aug best_stop: $best_stop best_orf_len: $best_orf_len cds length: ",length $cds_seq,"\n";

  # convert splices chain to exons - this includes the UTRs still as
  # we haven't done anything yet with the AUG and STOP codon positions
  my @exons = convert_splices_to_exons($region_start, $region_end, $chain, $sense);

  # strip out the UTRs if we have START and STOP codons in the first and last exons
  # get the chromosomal positions of the START and STOP codons
  my ($chrom_aug, $chrom_stop, $valid_aug, $valid_stop, @new_exons) = add_aug_and_stop_to_structure($region_start, $region_end, $best_aug, $best_stop, \@exons, $sense);
  print "valid_aug: $valid_aug valid_stop: $valid_stop\n";

  # is this a valid ORF?
  if ($valid_aug && $valid_stop) { # START and STOP codons in first and last introns
    $type = 'CDS';
  } else {
    $type = 'INVALID'; # possibly a pseudogene or non-coding transcript - START or STOP inside structure
  }
  
  return ($chrom_aug, $chrom_stop, $type, \@new_exons);
}

###############################################################################
# convert splices chain to exons - this includes the UTRs still
sub convert_splices_to_exons {
  my ($region_start, $region_end, $chain, $sense) = @_;
  my @exons;

  print "region_start $region_start, region_end $region_end sense $sense\n";
  if ($sense eq '+') {
    my $prev_end = $region_start; # start of first exon
    foreach my $splice (@{$chain}) {
      if (exists $splice->{head}) {next}
      my $seq = $splice->{seq};
      my $splice_start = $splice->{start};
      my $splice_end = $splice->{end};
      
      push @exons, {
		    start => ($prev_end - $region_start + 1),
		    end   => ($splice_start - $region_start),
		   };
      $prev_end = $splice_end+1; # start of next exon
    }
    push @exons, {
		  start => ($prev_end - $region_start + 1),
		  end   => ($region_end - $region_start + 1),
		 };
  } else {
    my $stored_end = $region_end - $region_start + 1;
    foreach my $splice (@{$chain}) {
      if (exists $splice->{head}) {next}
      my $splice_start = $splice->{start};
      my $splice_end = $splice->{end};
      print "splice_start $splice_start splice_end $splice_end\n";
      push @exons, {
		    start => $region_end - $splice_start + 2,
		    end   => $stored_end,
		   };
      $stored_end = $region_end - $splice_end;
    }
    push @exons, {
		  start => 1,
		  end   => $stored_end,
		 };
    @exons = reverse @exons; 
  }

  return @exons;
}
###############################################################################
# strip out the UTRs if we have START and STOP codons in the first and last exons
# get the chromosomal positions of the START and STOP codons
sub add_aug_and_stop_to_structure {
  my ($region_start, $region_end, $best_aug, $best_stop, $exons, $sense) = @_;
  my ($chrom_aug, $chrom_stop, $valid_aug, $valid_stop, @new_exons);

  $valid_aug = 0; # not seen it yet
  $valid_stop = 0;

  $best_aug++; # we are working in first base == 1 coords here
  $best_stop++;

  my $start;
  my $end;
  my $cumulative_len = 0;
  my $previous_cumulative_len = 0;
  my $count = 0;
  foreach my $exon (@{$exons}) {
    $start = $exon->{start};
    $end = $exon->{end};
    my $len = $end-$start+1;
    $cumulative_len += $len;

    if ($valid_aug) {
      $start -= $best_aug - 1;
      $end -= $best_aug - 1;
    }
    
    if ($best_aug >= $start && $best_aug <= $end && $count == 0) { # is AUG in first exon?
      $valid_aug = 1;
      $start = 1;
      $end -= $best_aug - 1;
    }

    print "count: $count best_stop: $best_stop previous_cumulative_len: $previous_cumulative_len cumulative_len: $cumulative_len\n";
    if ($best_stop > $previous_cumulative_len && $best_stop <= $cumulative_len && $count == $#{$exons}) { # is STOP in last exon?
      $valid_stop = 1;
      $end -=  ($cumulative_len - $best_stop);
    }

    push @new_exons, {
		start => $start,
		end   => $end,
	       };

    $previous_cumulative_len = $cumulative_len;
    $count++; # exon number
  }

  if ($valid_aug && $valid_stop) {
    if ($sense eq '+') {
      $chrom_aug = $region_start + $best_aug - 1;
      $chrom_stop = $region_end - ($cumulative_len - $best_stop);
    } else {
      $chrom_aug = $region_end - ($cumulative_len - $best_stop);
      $chrom_stop = $region_start + $best_aug - 1;
    }
  } else {
    @new_exons = @{$exons}; # no good START/STOP codon found, so use original structure as a pseudogene/ncRNA Transcript
    if ($sense eq '+') {
      $chrom_aug = $region_start; # use the region start and end as the structure object's start and end
      $chrom_stop = $region_end;
    } else {
      $chrom_aug = $region_end;
      $chrom_stop = $region_start; # use the region start and end as the structure object's start and end
    }
  }

  return ($chrom_aug, $chrom_stop, $valid_aug, $valid_stop, @new_exons)
}

###############################################################################
# check to see if the clone start-end region and exon start-end pairs
# match those of a 'curated' object or a 'isoformer' CDS object in the database

sub check_not_already_curated {
  my ($biotype, $clone, $clone_aug, $clone_stop, $exons) = @_;
  my $name = '';

  # parent clone coords
  my $clone_obj = $db->fetch(Sequence => "$clone");
  my @clone_entry;
  if ($biotype eq 'CDS') {
    @clone_entry = $clone_obj->CDS_child;
  } elsif ($biotype eq 'INVALID') {
    @clone_entry = $clone_obj->Pseudogene;
  }
  my $start;
  my $end;
  my $obj;
  foreach my $CDS ( @clone_entry ) {
    $start = $CDS->right->name;
    $end = $CDS->right->right->name;
    if ($start != $clone_aug || $end != $clone_stop) {next}
    if ($CDS->Method->name ne 'curated' && 
	$CDS->Method->name ne 'Pseudogene' && 
	$CDS->Method->name ne 'isoformer' &&
	$CDS->Method->name ne 'non_coding_isoformer'
       ) {next} # don't want to report a match to a history object or an ab-initio prediction etc.
    $name = $CDS->name;
    if (&check_exons_match($CDS, $exons)) {return $name}
  }

  return '';
}

###############################################################################
# check to see if the exon start-end pairs match those of an object in the database

sub check_exons_match {
  my ($obj, $exons) = @_;

  my $count = 0;
  foreach ($obj->Source_exons) {
    if ($#{$exons} < $count) {
      return 0;
    }
    my ($start, $end) = $_->row(0);
    if ($start != $exons->[$count]{start} || $end != $exons->[$count]{end}) {
      return 0;
    }
    $count++;
  }
  if ($#{$exons} > $count-1) {
    return 0;
  }
  
  return 1;
}
###############################################################################
sub increment_alternate_splicing {
  my ($chain) = @_;
  my $finished=0;

  my $last_pos = $#{$chain};

  while (1) {
    my $splice = $chain->[$last_pos];
    my $number_of_child_nodes = $#{$splice->{list}};
    $splice->{pos}++; # point to the next alternate child node
    if ($splice->{pos} > $number_of_child_nodes) {
      $splice->{pos} = 0; # reset the pointer
      $last_pos--;
      if ($last_pos < 0) {
	$finished=1; # no more levels of alternate splice positions to go up so we are finished
	last;
      }
    } else {
      last; # stop incrementing
    }
  }

  return $finished;
}
###############################################################################
# get clone and coords that encapsulate the given chromosome coords

sub get_clone_coords {
  my ($chromosome, $chrom_aug, $chrom_stop, $sense) = @_;

  my @clone_coords = $coords->LocateSpan($chromosome, $chrom_aug, $chrom_stop ); 
  my ($clone, $clone_aug, $clone_stop) = @clone_coords;

  return ($clone, $clone_aug, $clone_stop);
}

###############################################################################
# make a temporary gene

sub make_isoform {
  my ($Method, $TSL, $clone, $clone_aug, $clone_stop, $sense, $exons, $next_isoform) = @_;

  my $name = "${Method}_${next_isoform}";
  my $output = "$database/tmp/isoformer_isoform$$";
  my $feature_id = $TSL->{id}; # undef if not using a TSL
  my $TSL_type = $TSL->{tsl}; # empty string if not using a TSL

  my $biotype = 'CDS';
  my $clone_tag = 'CDS_child';
  my $lab = 'HX';
  my $g_species = $wormbase->full_name();
  my $gene;
  my $pseudogene_type;
  my $transcript_type;

  if ($Method eq 'isoformer') {
    $biotype = 'CDS';
    $clone_tag = 'CDS_child';
  } elsif ($Method eq 'non_coding_isoformer') {
    $biotype = 'Transcript';
    $clone_tag = 'Transcript';
    $transcript_type = 'ncRNA';
  } else {
    die "The Method '$Method' is not recognised\n";
  }

  my $personid;
  if ($USER eq 'gw3') {
    $personid = 'WBPerson4025';
  } elsif ($USER eq 'pad') {
    $personid = 'WBPerson1983';
  } elsif ($USER eq 'es9') {
    $personid = 'WBPerson2062';
  } else {
    die "Unknown WBPerson id: $USER\n";
  }

  my ($day, $mon, $yr)  = (localtime)[3,4,5];
  my $date = sprintf("%02d%02d%02d",$yr-100, $mon+1, $day);


  if (!-e "$database/tmp") {mkdir "$database/tmp", 0777};

  #print ace format
  open (HIS,">$output") or die "cant open $output\n";
  
  # delete any existing instance of this temporary CDS
  print HIS "\n-D $biotype : \"$name\"\n";
  close HIS;
  my $return_status = system("xremote -remote 'parse $output'");
  if ( ( $return_status >> 8 ) != 0 ) {
    print STDERR "WARNING - X11 connection appears to be lost\n";
  }

  # make the object
  open (HIS,">$output") or die "cant open $output\n";

  print HIS "\nSequence : $clone\n";
  print HIS "$clone_tag \"$name\" $clone_aug $clone_stop\n";
  print "$clone_tag \"$name\" $clone_aug $clone_stop\n";
  
  print HIS "\n$biotype : \"$name\"\n";
  foreach my $exon (@{$exons}) {
    print HIS "Source_exons ", $exon->{start}, " ", $exon->{end},"\n";
    print "Source_exons ", $exon->{start}, " ", $exon->{end},"\n";
  }
  my $remark = "Remark \"[$date $USER] Created this isoform based on the RNASeq data";
  if ($TSL_type) {
    $remark .=  " and the $TSL_type site"
  }
  $remark .= ".\"";
  print HIS "$remark Curator_confirmed $personid\n";
  print HIS "$remark From_analysis RNASeq_Hillier_elegans\n";
  if ($TSL_type) {
    print HIS "$remark Feature_evidence $feature_id\n";
  }
  print HIS "Sequence $clone\n";
  print HIS "From_laboratory $lab\n";
  print HIS "Gene $gene\n" if $gene;
  print HIS "Species \"$g_species\"\n";
  print HIS "Method $Method\n";
  print HIS "CDS\n" if ($biotype eq 'CDS');
  print HIS "$pseudogene_type\n" if ($pseudogene_type);
  print HIS "Transcript $transcript_type\n" if ($transcript_type);

  close HIS;
#  sleep(1); # allow time for NFS torealize there is a file there
  $return_status = system("xremote -remote 'parse $output'");
  if ( ( $return_status >> 8 ) != 0 ) {
    print STDERR "WARNING - X11 connection appears to be lost\n";
    #      &error_warning("WARNING", "X11 connection appears to be lost");
  } else {
    print "Made $name\n";
    #      &confirm_message("Made history","History $cds:${wormpep_prefix}$version has been created");
    #      &clear;
  }

  #system("rm -f $output");

}
