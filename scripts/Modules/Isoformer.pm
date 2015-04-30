#!/software/bin/perl -w
#
# Isoformer.pm                       
# 
# by Gary Williams
#
# Methods for running the Isoformer pipeline and other useful things
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2015-04-30 14:02:39 $      

=pod

=head1 NAME

 Isoformer

=head1 SYNOPSIS

  use Isoformer;

  # initialise
  my $Isoformer = Isoformer->new($wormbase, $log, $database);


=head1 DESCRIPTION

Isoformer stuff.

=head1 CONTACT

Gary gw3@ebi.ac.uk

=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are preceded with a _

=cut




package Isoformer;


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

=head2 

    Title   :   new
    Usage   :   my $Isoformer = Isoformer->new($wormbase, $log, $database);
    Function:   initialises the object
    Returns :   Isoformer object;
    Args    :   $wormbase
                $log
                $database - path to database to get sequences, structures and genes from
                $gff - optional location of the GFF file if it is not in the normal place
                $notsl - optional - if true then don't try to make TSL isoforms - for debugging purposes

$database can be undefined

=cut

sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;

  $self->{'wormbase'} = shift;    
  $self->{'log'} = shift;
  $self->{'database'} = shift;
  $self->{gff} = shift;
  $self->{notsl} = shift;
  
  $self->{species} = $self->{wormbase}->species;
  my $species = $self->{species};

  my $USER = $ENV{USER};


  if ($species eq 'elegans') {
    $self->{curated_database} = "/nfs/wormpub/camace_${USER}";
  } else {
    $self->{curated_database} = "/nfs/wormpub/${species}_curation";
  }


  $self->{coords} = Coords_converter->invoke($self->{database}, undef, $self->{wormbase});
  $self->db_connect();

  mkdir "$self->{database}/CHROMOSOMES", 0777; # we may need to write a new set of chromosome files if they do not yet exist
  $self->{seq_obj} = Sequence_extract->invoke($self->{database}, undef, $self->{wormbase});

  # name of the method and base-name of the objects
  $self->{CDS_name} = "isoformer";
  $self->{ncRNA_name} = "non_coding_transcript_isoformer";
  $self->{not_confirmed_name} = "not_confirmed_isoformer";

  # get all RNASeq splices starting and ending in the region
  # get all SL1, SL2 sites
  my %all_splices;
  my %all_TSL;
  print "Reading splice data\n";
  my $data_dir = $self->{wormbase}->wormpub . "/CURATION_DATA";
  my $dump_splices = "$data_dir/isoformer_${species}_splices.dat";
  my $dump_TSL = "$data_dir/isoformer_${species}_TSL.dat";
  
  if (-e $dump_splices && -e $dump_TSL) {
    $self->{wormbase}->FetchData("isoformer_${species}_splices", \%all_splices, $data_dir);
    $self->{wormbase}->FetchData("isoformer_${species}_TSL", \%all_TSL, $data_dir);
    
  } else {
    foreach my $chromosome ($self->{wormbase}->get_chromosome_names(-mito => 1, -prefix => 1)) {
      if (!exists $all_splices{$chromosome}) {
	my ($splice_data, $TSL_data) = $self->read_splice_and_TSL_data($self->{gff}, $chromosome, $self->{log});
	$all_splices{$chromosome} = $splice_data;
	$all_TSL{$chromosome} = $TSL_data;
      }
    }
    # now dump data to file
    open (all_splices, ">$dump_splices") or die "Can't write to $dump_splices\n";
    print all_splices Data::Dumper->Dump([\%all_splices]);
    close all_splices;
    open (all_TSL, ">$dump_TSL") or die "Can't write to $dump_TSL\n";
    print all_TSL Data::Dumper->Dump([\%all_TSL]);
    close all_TSL;
  }

  $self->{all_splices} = \%all_splices;
  $self->{all_TSL} = \%all_TSL;

  $self->interactive(1); # default is to allow interactive actions

  $self->next_isoform(1); # used in constructing the unique name of the object

  $self->aceclear(); # initialise the output ACE string

  return $self;
}

########################################################################

# setter/getter routines for paths etc.
sub db          { my $self = shift; my $db = shift; if (defined $db) {$self->{db} = $db}; return $self->{'db'}; }
sub database    { my $self = shift; my $db = shift; if (defined $db) {$self->{database} = $db}; return $self->{'database'}; }
sub interactive { my $self = shift; my $inter = shift;  if (defined $inter) {$self->{interactive} = $inter}; return $self->{'interactive'}; }
#sub next_isoform { my $self = shift; my $next_isoform = shift;  if (defined $next_isoform) {$self->{next_isoform} = $next_isoform} else {$self->{'next_isoform'} += 1}; return $self->{'next_isoform'}; }
sub next_isoform { my $self = shift; my $next_isoform = shift;  if (defined $next_isoform) {$self->{next_isoform} = $next_isoform}; return $self->{'next_isoform'}; }
sub inc_next_isoform { my $self = shift; $self->{'next_isoform'} += 1; return $self->{'next_isoform'}; }
sub aceout      { my $self = shift; my $ace = shift;  if (defined $ace) {$self->{aceout} .= $ace}; return $self->{'aceout'}; }
sub aceclear    { my $self = shift; $self->{aceout} = '';}
sub db_connect  {my $self = shift;   $self->{db} = Ace->connect(-path => $self->{database}) || die "cannot connect to database at $self->{database}\n";}
sub db_close    {my $self = shift; $self->{db}->close}

########################################################################


=head2 

    Title   :   make_isoforms_in_region
    Usage   :   my ($confirmed, $not_confirmed, $created, $warnings) = $Iso->make_isoforms_in_region($chromosome, $region_start, $region_end, $sense, $biotype, $gene);
    Function:   makes isoforms in the given region
    Returns :   arrayrefs of confirmed, not_confirmed, created and warning messages
    Args    :   $biotype - biotype of gene defining region - if known, it is one of 'CDS', 'Pseudogene', Transposon', 'Transcript'

=cut


sub make_isoforms_in_region {
  my ($self, $sequence, $seq_start, $seq_end, $sense, $biotype, $gene) = @_;

  # convert from clone coords etc to chromosome coords
  my ($chromosome, $region_start, $region_end);
  ($chromosome, $region_start) = $self->{coords}->Coords_2chrom_coords($sequence, $seq_start);
  ($chromosome, $region_end) = $self->{coords}->Coords_2chrom_coords($sequence, $seq_end);

   
  my @TSL = &get_TSL_in_region($chromosome, $region_start, $region_end, $sense, $self->{all_TSL}{$chromosome});
  
  # add the $region_start/end to the head of the TSL list as a dummy TSL site as it is one of the positions that enclose the regions to consider
  if ($sense eq '+') {
    unshift @TSL, {
		seq => $chromosome,
		start => $region_start,
		end => $region_start,
		sense => $sense,
		tsl => '',
	       }
  } else {
    unshift @TSL, {
		seq => $chromosome,
		start => $region_end,
		end => $region_end,
		sense => $sense,
		tsl => '',
	       }
  }

  if ($self->interactive()) {
    $self->next_isoform(1); # used in constructing the unique name of the object
  }
  my %confirmed = ();
  my @created = ();
  my @warnings = ();
  my %structures_made = ();


  foreach my $TSL (@TSL) {
    if ($sense eq '+') {
      $region_start = $TSL->{start};
    } else {
      $region_end = $TSL->{start};
    }
#    print "********** in foreach TSL region_end = $region_end\n";

    my @splices = &get_splices_in_region($chromosome, $region_start, $region_end, $sense, $self->{all_splices}{$chromosome});
    
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
    my ($created, $warnings) = $self->iterate_through_the_valid_structures(\%structures_made, $gene, $biotype, \%confirmed, $TSL, $chromosome, $region_start, $region_end, $sense, @splices);
    push @created, @{$created};
    push @warnings, @{$warnings};
    
  } # foreach TSL
  
  # print a summary
  my @confirmed=();
  my @not_confirmed=();
  my $g_species = $self->{wormbase}->full_name();

  foreach my $confirmed_name (keys %confirmed) {
    if ($confirmed{$confirmed_name}) {
      push @confirmed, $confirmed_name;

      # add the Evidence to show that this looks good
      $self->aceout("\n$biotype : $confirmed_name\n");
      $self->aceout("Evidence Inferred_automatically make_isoforms.pl\n");

    } else {
      push @not_confirmed, $confirmed_name;
      
      my $not_confirmed_name = $self->{not_confirmed_name};
      my $method = $self->{not_confirmed_name};
      my $nc_name = "${confirmed_name}_${not_confirmed_name}_".$self->inc_next_isoform();
      my $last = $seq_end - $seq_start + 1;

      $self->aceout("\nCDS : $nc_name\n");
      $self->aceout("Evidence Inferred_automatically make_isoforms.pl\n");
      $self->aceout("Method $method\n");
      $self->aceout("Source_exons 1 $last\n");
      $self->aceout("Species \"$g_species\"\n");
      $self->aceout("CDS\n");


      $self->aceout("\nSequence : $sequence\n");
      if ($sense eq '+') {
	$self->aceout("CDS_child \"$nc_name\" $seq_start $seq_end\n");
      } else {
	$self->aceout("CDS_child \"$nc_name\" $seq_end $seq_start\n");
      }
    }
  }

  return (\@confirmed, \@not_confirmed, \@created, \@warnings);

}

########################################################################


##############################################################
#
# Subroutines
#
##############################################################

# load the method to display the isoformer CDS structures

sub load_isoformer_method {
  my ($self) = @_;


  $self->aceout("\n\n");
  $self->aceout("Method : \"".$self->{CDS_name}."\"\n");
  $self->aceout("Colour DARKGREEN\n");
  $self->aceout("CDS_colour DARKGREEN\n");
  $self->aceout("Strand_sensitive Show_up_strand\n");
  $self->aceout("Right_priority   2.63\n");
  $self->aceout("Overlap_mode Bumpable\n");
  $self->aceout("Show_text\n");
  $self->aceout("EMBL_feature CDS\n");
  $self->aceout("Remark \"This method is used by acedb to display curation CDS structures generated by isoformer.pl\"\n");
  $self->aceout("\n\n");
  $self->aceout("\n\n");
  $self->aceout("Method : \"".$self->{ncRNA_name}."\"\n");
  $self->aceout("Colour ORANGE\n");
  $self->aceout("CDS_colour ORANGE\n");
  $self->aceout("Strand_sensitive Show_up_strand\n");
  $self->aceout("Right_priority   2.6\n");
  $self->aceout("Overlap_mode Bumpable\n");
  $self->aceout("Show_text\n");
  $self->aceout("EMBL_feature ncRNA\n");
  $self->aceout("Remark \"This method is used by acedb to display curation non-coding structures (ncRNA/Pseudogene/non_coding_transcript_isoform) generated by isoformer.pl\"\n");
  $self->aceout("\n\n");
  $self->aceout("\n\n");
  $self->aceout("Method : \"".$self->{not_confirmed_name}."\"\n");
  $self->aceout("Colour RED\n");
  $self->aceout("CDS_colour RED\n");
  $self->aceout("Strand_sensitive Show_up_strand\n");
  $self->aceout("Right_priority   1.45\n");
  $self->aceout("Overlap_mode Bumpable\n");
  $self->aceout("Show_text\n");
  $self->aceout("EMBL_feature CDS\n");
  $self->aceout("Remark \"This method is used by acedb to display curation CDS structures that are NOT confirmed by isoformer.pl\"\n");
  $self->aceout("\n\n");

}

###############################################################################
# return the chromosomal location and sense of a CDS or a clone:start-end triplet

# get the region of interest from the CDS name or clone positions
# region can be specified as either a CDS name or a 'seq start end' triplet
# with either space or ':' or '-' separators

# reverse sense regions are indicated by having the end < start these
# are reversed when the positions are returned, giving start end in
# the forward sense with a sense character

# returns biotype as one of 'CDS', 'Pseudogene', 'Transcript'

sub get_active_region {
  my ($self, $region) = @_;
  my ($chromosome, $start, $end, $sense, $biotype, $gene);
  $gene = '';
  $sense = '+';

  my $coords = $self->{coords};

  my @region = split /[\s]+/, $region;
  if ($region[0]) { # CDS or Pseudogene
    my $obj = $self->{db}->fetch(CDS => "$region[0]");
    if (defined $obj) {
      my $gene_obj = $obj->Gene;
      if (defined $gene_obj) {
	$gene = $obj->Gene->name;
      }
      my $clone_start;
      my $clone_end;
      my $clone = $obj->Sequence;
      if (!defined $clone) {
	print STDERR "Sequence is not defined for the CDS object $region[0]\n";
	return undef;
      }
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
      my $obj = $self->{db}->fetch(Pseudogene => "$region[0]");
      if (defined $obj) {
	my $gene_obj = $obj->Gene;
	if (defined $gene_obj) {
	  $gene = $obj->Gene->name;
	}
	my $clone_start;
	my $clone_end;
	my $clone = $obj->Sequence;
	if (!defined $clone) {
	  print STDERR "Sequence is not defined for the Pseudogene object $region[0]\n";
	  return undef;
	}
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
	
      } else { # try Transcript
	my $obj = $self->{db}->fetch(Transcript => "$region[0]");
	if (! defined $obj) {
	  print STDERR "Region: '$region' is an unknown way of specifying a region: not a CDS or a Pseudogene or a Transcript\n";
	  return undef;
	}
	my $gene_obj = $obj->Gene;
	if (defined $gene_obj) {
	  $gene = $obj->Gene->name;
	}
	my $clone_start;
	my $clone_end;
	my $clone = $obj->Sequence;
	if (!defined $clone) {
	  print STDERR "Sequence is not defined for the Transcript object $region[0]\n";
	  return undef;
	}
	my @clone_Transcript = $clone->Transcript;
	foreach my $Transcript ( @clone_Transcript ) {
	  next unless ($Transcript->name eq "$region[0]");
	  $clone_start = $Transcript->right->name;
	  $clone_end = $Transcript->right->right->name;
	  last;
	}
	($chromosome, $start) = $coords->Coords_2chrom_coords($clone, $clone_start);
	($chromosome, $end) = $coords->Coords_2chrom_coords($clone, $clone_end);
	$biotype = 'Transcript';
	
      }
    }

    if ($region[1]) {
      $start += $region[1];
    }

    if ($region[2]) {
      $end = $end + $region[2];
    }

  } else {
    print STDERR "Region: '$region' is an unknown way of specifying a region";
    return undef;
  }

  if ($start > $end) {
    ($start, $end) = ($end, $start);
    $sense = '-';
  }

  return ($chromosome, $start, $end, $sense, $biotype, $gene);
}
###############################################################################
# read the RNASeq intron data and TSL site data from the GFF file

sub read_splice_and_TSL_data {
  my ($self, $gff, $chromosome, $log) = @_;
  my @splice_data;
  my @TSL_data;
  my %splice;
  my $splice_method = 'RNASeq_splice';
  my $SL1_method = 'SL1';
  my $SL2_method = 'SL2';
  my $file = $self->open_GFF_file($gff, $chromosome, $log);

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
  my ($self, $gff, $chromosome, $log) = @_;
  my $file;
  my $fh;

  if (defined $gff) {
    $file = $gff;
  } else {

    if ($self->{species} eq 'elegans') {
      $file = $self->{wormbase}->autoace."/CHROMOSOMES/$chromosome.gff";
    } else {
      if ($self->{wormbase}->assembly_type eq 'contig') {      
	$file = $self->{wormbase}->sequences . "/" . $self->{species}.".gff";
	
      } else {
	$file =  $self->{wormbase}->autoace."/CHROMOSOMES/$chromosome.gff";
      }
    }
  }
  open($fh, "cat $file | grep \"^$chromosome\\W\" |") or $log->log_and_die("Could not open seq-restricted stream to $file\n");
  return $fh;

}
###############################################################################
# filter the list of splices and return those enclosed by the specied region

sub get_splices_in_region {
  my ($chromosome, $start, $end, $sense, $all_splices) = @_;
  my @splices;

  foreach my $splice (@{$all_splices}) {
    if ($splice->{start} <= $start) {next}
    if ($splice->{end} >= $end) {next}
    if ($splice->{sense} ne $sense) {next}
    if ($splice->{seq} ne $chromosome) {next}
    push @splices, $splice;
  }

  return @splices;
}
###############################################################################
# filter the list of TSLs and return those enclosed by the specied region

sub get_TSL_in_region {
  my ($chromosome, $start, $end, $sense, $all_TSL) = @_;
  my @TSL;

  foreach my $TSL (@{$all_TSL}) {
    if ($TSL->{start} < $start) {next}
    if ($TSL->{end} > $end) {next}
    if ($TSL->{sense} ne $sense) {next}
    if ($TSL->{seq} ne $chromosome) {next}

    # now see if we have two TSLs at the same position
    my $notfound = 1;
    foreach my $found_tsl (@TSL) {
      if ($found_tsl->{start} != $TSL->{start}) {next}
      if ($found_tsl->{end} != $TSL->{end}) {next}
      if ($found_tsl->{sense} ne $TSL->{sense}) {next}
      if ($found_tsl->{seq} ne $TSL->{seq}) {next}
      $notfound = 0;
      $found_tsl->{tsl} .= (" ".$TSL->{tsl});
      $found_tsl->{id}  .= (" ".$TSL->{id});
    }
    if ($notfound) {
      push @TSL, $TSL;
    }
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
  my ($self, $structures_made, $gene, $biotype, $confirmed, $TSL, $chromosome, $region_start, $region_end, $sense, @splices) = @_;
  
  my $finished=0; # set when we have no more levels to increment up
  
  my @created = ();
  my @warnings = ();
  
  # get the clone coords of the region start/end
  #  print "********** in iterate_through_the_valid_structures region_end = $region_end\n";
  my ($region_clone, $region_clone_start, $region_clone_end) = $self->get_clone_coords($chromosome, $region_start, $region_end, $sense);
  # get any structures in this initial region
  $self->get_structures_in_region($confirmed, $sense, $region_clone, $region_clone_start, $region_clone_end);
  
  do {
    
    my ($chain) = &get_next_chain_of_intron_nodes(@splices);
    
    my ($chrom_aug, $chrom_stop, $type, $exons) = $self->convert_to_exons($region_start, $region_end, $chain, $sense);
    # avoid bug where we have one exon of length one base
    if (scalar @{$exons} == 1 && $exons->[0]->{start} == $exons->[0]->{end}) {
      print "SKIPPING STRUCTURE CONSISTING OF A OF 1 BASE EXON\n";
    } else {
      if (check_for_duplicate_structures($chrom_aug, $chrom_stop, $type, $exons, $structures_made)) {
	print "Duplicate structure ignored\n";
      } else {
	my $structure_warning = sanity_check($type, $chrom_aug, $chrom_stop, $exons);
	if ($structure_warning) {push @warnings, $structure_warning}
	
	if ($type eq 'CDS') {
	  my ($clone, $clone_aug, $clone_stop) = $self->get_clone_coords($chromosome, $chrom_aug, $chrom_stop, $sense);
	  
	  # get structures overlapping the isoformer structure - doing
	  # this again in case we are at a clone boundary and some
	  # structures are on the superlink and some on the clone
	  $self->get_structures_in_region($confirmed, $sense, $clone, $clone_aug, $clone_stop);
	  
	  my $confirmed_ok = $self->check_not_already_curated($confirmed, $sense, $type, $clone, $clone_aug, $clone_stop, $exons);
	  if (!$confirmed_ok) {	  
	    push @created, $self->make_isoform($self->{CDS_name}, $TSL, $clone, $clone_aug, $clone_stop, $sense, $exons, $gene);
	    $self->inc_next_isoform();
	  }
	} elsif ($type eq 'INVALID') { # no START and STOP found in the first and last exons
	  print "Found no good coding structure - should the region to look at be larger?\n";
	  my ($clone, $clone_aug, $clone_stop) = $self->get_clone_coords($chromosome, $chrom_aug, $chrom_stop, $sense);
	  
	  # get structures overlapping the isoformer structure - doing
	  # this again in case we are at a clone boundary and some
	  # structures are on the superlink and some on the clone
	  $self->get_structures_in_region($confirmed, $sense, $clone, $clone_aug, $clone_stop);
	  
	  my $confirmed_ok = $self->check_not_already_curated($confirmed, $sense, $type, $clone, $clone_aug, $clone_stop, $exons);
	  if (!$confirmed_ok) {
	    push @created, $self->make_isoform($self->{ncRNA_name}, $TSL, $clone, $clone_aug, $clone_stop, $sense, $exons, $gene);
	    $self->inc_next_isoform();
	  }
	  
	} else {
	  die "Type '$type' not defined\n\n";
	}
      }
    }
      
    $finished = &increment_alternate_splicing($chain);
  } until ($finished);

  return (\@created, \@warnings);
}

###############################################################################
# check against the hash of structures already made in this region and return false if not made before

sub check_for_duplicate_structures {
  my ($chrom_aug, $chrom_stop, $type, $exons, $structures_made) = @_;
  
  my $result;
  
  # make key from exon boundaries
  my $key = join ':', ($chrom_aug, $chrom_stop, $type);
  foreach my $ex (@{$exons}) {
    $key .= (":" . $ex->{start} . ":" . $ex->{end});
  }
  #print "key = $key\n";
  
  # see if this exists already in hash
  if (exists $structures_made->{$key}) {$result = 1} else {$result = 0}
  
  # add key to hash
  $structures_made->{$key} = 1;
  
  return $result;
}

###############################################################################
  # sanity check that the structure is ok

sub sanity_check {
  
  my ($type, $chrom_aug, $chrom_stop, $exons) = @_;
  
  my $warning = "";

  my $length_in_sequence = abs($chrom_aug - $chrom_stop) + 1;
  my $max_exon = 0;
  foreach my $exon (@{$exons}) {
    if ($exon->{end} > $max_exon) {$max_exon = $exon->{end}}
  }
  if ($length_in_sequence != $max_exon) {
    $warning = "WARNING - $type sequence coords $chrom_aug-$chrom_stop have length: $length_in_sequence, but exons have length: $max_exon\n";
  }
 
  return $warning;
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
  my ($self, $region_start, $region_end, $chain, $sense) = @_;
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
    $cds_seq .= $self->{seq_obj}->Sub_sequence($seq, $prev_end-1, $len);
    $prev_end = $splice_end+1; # start of next exon
  }
  $cds_seq .= $self->{seq_obj}->Sub_sequence($chain->[0]{seq}, $prev_end-1, $region_end-$prev_end+1);
  if ($sense eq '-') {
    $cds_seq = $self->{seq_obj}->DNA_revcomp($cds_seq);
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
  my $best_aug = -1;
  my $best_stop = -1;
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
      my $score = $splice->{score};
      
      push @exons, {
		    start => ($prev_end - $region_start + 1),
		    end   => ($splice_start - $region_start),
		    score => $score,
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
      my $score = $splice->{score};

      print "splice_start $splice_start splice_end $splice_end\n";
      push @exons, {
		    start => $region_end - $splice_start + 2,
		    end   => $stored_end,
		    score => $score,
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

  # $best_aug and $best_stop will now be zero if they were not seen
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
    
    if ($best_aug ne 0 && $best_aug >= $start && $best_aug <= $end && $count == 0) { # is AUG in first exon?
      $valid_aug = 1;
      $start = 1;
      $end -= $best_aug - 1;
    }

    print "exon count: $count best_stop: $best_stop previous_cumulative_len: $previous_cumulative_len cumulative_len: $cumulative_len\n";
    if ($best_stop ne 0 && $best_stop > $previous_cumulative_len && $best_stop <= $cumulative_len && $count == $#{$exons}) { # is STOP in last exon?
      $valid_stop = 1;
      $end -=  ($cumulative_len - $best_stop);# + ($best_aug - 1); # add the ($best_aug - 1) back on as it was removed earlier
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
      $chrom_aug = $region_end - $best_aug + 1;
      $chrom_stop = $region_start + ($cumulative_len - $best_stop);
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
# get all structures that are in the region so that we can confirm them later

sub get_structures_in_region {
  my ($self, $confirmed, $sense, $clone, $clone_start, $clone_end) = @_;

#  print "*** in get_structures_in_region clone_end = $clone_end\n";

  # parent clone coords
  my $clone_obj = $self->{db}->fetch(Sequence => "$clone");
  my $start;
  my $end;
  my $obj;
  foreach my $structure ( $clone_obj->CDS_child, $clone_obj->Pseudogene ) {
    if (!defined $structure->Method) {print "Method not defined for object ".$structure->name."\n"; next}
    if ($structure->Method->name ne 'curated' && 
	$structure->Method->name ne 'Pseudogene'
       ) {next} # don't want to report a match to a history object or an ab-initio prediction etc.
    $start = $structure->right->name;
    $end = $structure->right->right->name;
    if (($sense eq '+' && ($start > $end)) || ($sense eq '-' && ($end > $start))) {next} # only want same sense
    if (($sense eq '+' && (($start >= $clone_start && $start <= $clone_end) || ($end >= $clone_start && $end <= $clone_end))) ||
	($sense eq '-' && (($start >= $clone_start && $start <= $clone_end) || ($end >= $clone_start && $end <= $clone_end)))
       ) { # is it in our region
      my $name = $structure->name;
#      print "*** $name is in the region\n";
      if (!exists $confirmed->{$name}) {$confirmed->{$name} = 0;}
    }
  }

}

###############################################################################
# check to see if the clone start-end region and exon start-end pairs
# match those of a 'curated' object or a 'isoformer' CDS object in the database

sub check_not_already_curated {
  my ($self, $confirmed, $sense, $biotype, $clone, $clone_aug, $clone_stop, $exons) = @_;
  my $name = '';

  # parent clone coords
  my $clone_obj = $self->{db}->fetch(Sequence => "$clone");
  my @clone_entry;
  if ($biotype eq 'CDS') {
    @clone_entry = $clone_obj->CDS_child;
  } elsif ($biotype eq 'INVALID') {
    $biotype = 'Pseudogene';
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
	$CDS->Method->name ne $self->{CDS_name} &&
	$CDS->Method->name ne $self->{ncRNA_name}
       ) {next} # don't want to report a match to a history object or an ab-initio prediction etc.
    $name = $CDS->name;
    print "checking exons against $name\n";
    if (&check_exons_match($CDS, $exons)) {
      $confirmed->{$name} = $biotype;
      print "Found all exons match, so $name is confirmed.\n";
      return $name;
    }
  }

  return '';
}

###############################################################################
# check to see if the exon start-end pairs match those of an object in the database

sub check_exons_match {
  my ($obj, $exons) = @_;

  my $count = 0;
  my @exon_coord1 = sort by_number ($obj->follow->get('Source_exons',1));
  my @exon_coord2 = sort by_number ($obj->follow->get('Source_exons',2));
  if ($#{$exons} != $#exon_coord1) {return 0;}
  foreach my $exon (@{$exons}) {
    #print "$exon_coord1[$count] $exon_coord2[$count] ", $exon->{start}," ",$exon->{end},"\n";
    if ($exon_coord1[$count] != $exon->{start} || $exon_coord2[$count] != $exon->{end}) {
      return 0;
    }
    $count++;
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
  my ($self, $chromosome, $chrom_aug, $chrom_stop, $sense) = @_;

  my @clone_coords = $self->{coords}->LocateSpan($chromosome, $chrom_aug, $chrom_stop ); 
  my ($clone, $clone_aug, $clone_stop) = @clone_coords;

  return ($clone, $clone_aug, $clone_stop);
}

###############################################################################
# make an isofrm structure
# return the name of the isoform that was made.

sub make_isoform {
  my ($self, $Method, $TSL, $clone, $clone_aug, $clone_stop, $sense, $exons, $gene) = @_;

  my $name;
  if ($self->{interactive}) {
    $name = "${Method}_".$self->next_isoform();
  } else {
    # make the whole genome run have different names to the interactive isoforms
    $name = $self->{species}."_${Method}_${gene}_".$self->next_isoform(); 
  }

  my $feature_id = $TSL->{id}; # undef if not using a TSL
  my $TSL_type = $TSL->{tsl}; # empty string if not using a TSL

  my $biotype = 'CDS';
  my $clone_tag = 'CDS_child';
  my $lab = 'HX';
  my $g_species = $self->{wormbase}->full_name();
  my $pseudogene_type;
  my $transcript_type;
  my $total_score;

  if ($Method eq $self->{CDS_name}) {
    $biotype = 'CDS';
    $clone_tag = 'CDS_child';
  } elsif ($Method eq $self->{ncRNA_name}) {
    $biotype = 'Transcript';
    $clone_tag = 'Transcript';
    $transcript_type = 'ncRNA';
  } else {
    die "The Method '$Method' is not recognised\n";
  }

  my $personid = $self->person();

  my ($day, $mon, $yr)  = (localtime)[3,4,5];
  my $date = sprintf("%02d%02d%02d",$yr-100, $mon+1, $day);


  # force a removal of any existing object with this number
  $self->clear("delete ".$self->next_isoform());

  $self->aceout("\nSequence : $clone\n");
  $self->aceout("$clone_tag \"$name\" $clone_aug $clone_stop\n");
  print "$clone_tag \"$name\" $clone_aug $clone_stop\n";

  $self->aceout("\n$biotype : \"$name\"\n");
  if (scalar @{$exons} == 1 && $exons->[0]->{start} == 1 && $exons->[0]->{end} == 1) { #  get an error sometimes where both of these are '1'
    print "HAVE THE ERROR WHERE THE SOURCE_EXONS ARE PRINTED AS 1 1 IN EXONS:\n@{$exons}\n";
  } else {
    foreach my $exon (@{$exons}) {
      $self->aceout("Source_exons ". $exon->{start}. " ". $exon->{end}."\n");
      print "Source_exons ". $exon->{start}. " ". $exon->{end}."\n";
      if (exists $exon->{score}) {$total_score += $exon->{score}}
    }
  }

  my $remark;
  my $DB_remark;
  my $Brief_identification;
  my $USER = $ENV{USER};
  if ($clone_tag eq 'Transcript') { 
    $remark = "Remark \"[$date $USER] Created this non-coding transcript isoform based on the RNASeq splice data";
    if ($total_score) { # will be zero in a one-exon gene, so don't bother to report it
      $remark .=  " from $total_score reads spanning introns";
    }
    if ($TSL_type) {
      $remark .=  " and the $TSL_type site"
    }
    $DB_remark = "DB_remark \"C. elegans probable non-coding RNA\"";
    $Brief_identification = "Brief_identification \"non-coding Transcript Isoform\"";
  }
  else {
    $remark = "Remark \"[$date $USER] Created this isoform based on the RNASeq splice data";
    if ($total_score) { # will be zero in a one-exon gene, so don't bother to report it
      $remark .=  " from $total_score reads spanning introns";
    }
    if ($TSL_type) {
      $remark .=  " and the $TSL_type site"
    }
  }
  $remark .= ".\"";
  $self->aceout("$remark Curator_confirmed $personid\n");
  $self->aceout("$remark From_analysis RNASeq\n");
  if ($TSL_type) {
    my @ids = split ' ', $feature_id;
    foreach my $id (@ids) {
      $self->aceout("$remark Feature_evidence $id\n");
      $self->aceout("Isoform Feature_evidence $id\n");
    }
  }
  $self->aceout("Isoform Curator_confirmed $personid\n");
  $self->aceout("$DB_remark\n") if $DB_remark;
  $self->aceout("$Brief_identification\n") if $Brief_identification;
  $self->aceout("Sequence $clone\n");
  $self->aceout("From_laboratory $lab\n");
  $self->aceout("Species \"$g_species\"\n");
  $self->aceout("Method $Method\n");
  $self->aceout("Gene $gene\n") if ($gene);
  $self->aceout("CDS\n") if ($biotype eq 'CDS');
  $self->aceout("$pseudogene_type\n") if ($pseudogene_type);
  $self->aceout("Transcript $transcript_type\n") if ($transcript_type);
  $self->aceout("Evidence Curator_confirmed $personid\n");

  return $name;
}
###############################################################################
# clear selected isoformer structures
# this accepts two command words: 
# clear - checks for existing structures to be deleted
# delete - just deletes without checking a single numeric parameter
sub clear {
  my ($self, $userinput) = @_;

  if (! $self->{interactive}) {return}

  my @f = split /\s+/, $userinput;
  my %todelete;

#  if ($f[0] eq 'clear') {
#    print "Please save your database and hit RETURN>";
#    my $answer = <STDIN>;
#    $self->{db}->close;
#    $self->{db} = Ace->connect(-path => $self->{database}) || die "cannot connect to database at $self->{database}\n";
#  }

  my $CDS_name = $self->{CDS_name};
  my $ncRNA_name = $self->{ncRNA_name};

  # get all the isoformer objects
  if ($f[0] eq 'clear' || $f[0] eq 'clean') {
    $f[0] = 'clear';
    my @CDS_objects = $self->{db}->fetch(CDS => "${CDS_name}_*");
    my @ncRNA_objects = $self->{db}->fetch(Transcript => "${ncRNA_name}_*");
    
    foreach my $obj (@CDS_objects) {
      $todelete{$obj->name} = 'CDS';
    }
    foreach my $obj (@ncRNA_objects) {
      $todelete{$obj->name} = 'Transcript';
    }
  } elsif ($f[0] eq 'delete') {
    $todelete{"${CDS_name}_$f[1]"} = 'CDS';
    $todelete{"${ncRNA_name}_$f[1]"} = 'Transcript';
    $f[1] = 'all';
  }

  if ($f[0] eq 'clear' && $#f > 0 && $f[1] ne 'all') { # if nothing after 'clear' then delete everything
    my $cmd = shift @f;
    my %explicitly_delete;
    foreach my $param (@f) {
      if ($param =~ /^${CDS_name}_\d+$/) {
	if (exists $todelete{$param}) {
	  $explicitly_delete{$param} = $todelete{$param};
	}
      }

      if ($param =~ /^${ncRNA_name}_\d+$/) {
	if (exists $todelete{$param}) {
	  $explicitly_delete{$param} = $todelete{$param};
	}
      }

      if ($param =~ /^\d+$/) {
	if (exists $todelete{"${CDS_name}_$param"}) {
	  $explicitly_delete{"${CDS_name}_$param"} = $todelete{"${CDS_name}_$param"};
	}
	if (exists $todelete{"${ncRNA_name}_$param"}) {
	  $explicitly_delete{"${ncRNA_name}_$param"} = $todelete{"${ncRNA_name}_$param"};
	}
      }
    }
    %todelete = %explicitly_delete;
  }

  # delete any existing instance of this temporary CDS
  foreach my $delete (keys %todelete) {
  
    my $biotype = $todelete{$delete};
    $self->aceout("\n-D $biotype : \"$delete\"\n");
    print "Clear $biotype : \"$delete\"\n";

    my $output = $self->{database}."/tmp/isoformer_$$.ace";
    open (CLEAR,">$output") or die "cant open $output\n";
    print CLEAR $self->aceout();
    $self->aceclear();
    close CLEAR;

    my $return_status = system("xremote -remote 'parse $output'");
    if ( ( $return_status >> 8 ) != 0 ) {
      die("WARNING - X11 connection appears to be lost\n");
    }

  }



}
###############################################################################
# what - report the isoformer objects in the database

sub what {
  my ($self, $userinput) = @_;

  if (! $self->{interactive}) {return}

  my @f = split /\s+/, $userinput;

  my @cds;
  my @transcript;

  my $CDS_name = $self->{CDS_name};
  my $ncRNA_name = $self->{ncRNA_name};

  # get all the isoformer objects
  my @CDS_objects = $self->{db}->fetch(CDS => "${CDS_name}_*");
  my @ncRNA_objects = $self->{db}->fetch(Transcript => "${ncRNA_name}_*");
  
  foreach my $obj (@CDS_objects) {
    push @cds, $obj->name;
  }
  foreach my $obj (@ncRNA_objects) {
    push @transcript, $obj->name;
  }

  if (@cds || @transcript) {
    print "\n*** The following objects are in the database:\n";
    print "@cds\n";
    print "@transcript\n";
  } else {
    print "\n*** No isoformer objects in the database.\n";
  }
  print "\n";

}

###############################################################################
# fix selected isoformer structures to a specified name
# fix isoformer_1 AC3.3b - fix specified object to CDS/Transcript, creating it if necessary

# check that the isoformer object exists
# check whether or not the target object exists
# check whether an existing target object is of the same class
# rename the isoformer object if there is no existing target of the same class
# else transfer location and tags to the existing object
# if the target object name ends with a letter, set the Isoform tag
# check that the Gene tag is populated correctly
# warn the user if the Gene tag is not populated
# warn the user if History should be made


sub fix {
  my ($self, $userinput) = @_;

  if (! $self->{interactive}) {return}

  my $new_name='';
  my %tofix;

  my @f = split /\s+/, $userinput;

  my $biotype;
  my $subject = $f[1];
  my $target = $f[2];
  my $subject_obj;
  my $target_obj;
  my $target_exists = 0;
  my $gene_obj;
  my $gene;
  my $g_species = $self->{wormbase}->full_name();

  my $CDS_name = $self->{CDS_name};
  my $ncRNA_name = $self->{ncRNA_name};

  if (!defined $subject || ! defined $target) {
    print "fix what to what?\n";
    return;
  }

  # check name looks OK and get biotype
  if ($subject =~ /^${ncRNA_name}_\d+/ || $subject =~ /${ncRNA_name}_WBGene\d+/) {
    $biotype = 'Transcript';
    $subject_obj = $self->{db}->fetch(Transcript => "$subject");
  } elsif ($subject =~ /^${CDS_name}_\d+/ || $subject =~ /${CDS_name}_WBGene\d+/) {
    $biotype = 'CDS';
    $subject_obj = $self->{db}->fetch(CDS => "$subject");
  } else {
    print "fix what?\n";
    return;
  }
  
  if (! defined $subject_obj) {
    print "\n*** Can't find '$subject' - have you saved your session?\n";
    return;
  }

  # test to see if the Gene tag is set in the subject
  $gene_obj = $subject_obj->Gene;
  if (defined $gene_obj) {
    $gene = $gene_obj->name;
  }

  # check if target exists
  my $message = '';
  if ($biotype eq 'CDS') {
    $target_obj = $self->{db}->fetch(CDS => "$target");
  } else {
    $target_obj = $self->{db}->fetch(Transcript => "$target");
  }
  if (! defined $target_obj) {
    $self->aceout("\n-R $biotype $subject $target\n");
    $self->aceout("\n$biotype : $target\n");
    if ($biotype eq 'CDS') {
      $self->aceout("Method curated\n");
    } else {
      $self->aceout("Method non_coding_transcript\n");
    }

    $message .= "\n*** Fix: $biotype $target does not exist - creating it.\n";

  } else {
    $target_exists = 1;

    # delete the target location
    $self->aceout("\n$biotype : $target\n");
    $self->aceout("-D Source_exons\n\n");

    my $target_sequence_obj = $target_obj->Sequence;
    my $target_sequence_name = $target_sequence_obj->name;
    $self->aceout("\nSequence : $target_sequence_name\n");
    if ($biotype eq 'CDS') {
      $self->aceout("-D CDS_child $target\n\n");
    } else {
      $self->aceout("-D Transcript $target\n\n");
    }

    # get the Sequence
    my $subject_sequence_obj = $subject_obj->Sequence;
    my $subject_sequence_name = $subject_sequence_obj->name;

    # get the subject location
    my @clone_locations;
    my $subject_start;
    my $subject_end;
    if ($biotype eq 'CDS') {
      @clone_locations = $subject_sequence_obj->CDS_child;
    } else {
      @clone_locations = $subject_sequence_obj->Transcript;
    }    
    foreach my $subject_location ( @clone_locations ) {
      next unless ($subject_location->name eq $subject);
      $subject_start = $subject_location->right->name;
      $subject_end = $subject_location->right->right->name;
      last;
    }

    # write the new target location
    $self->aceout("\nSequence : $subject_sequence_name\n");
    if ($biotype eq 'CDS') {
      $self->aceout("CDS_child $target $subject_start $subject_end\n\n");
    } else {
      $self->aceout("Transcript $target $subject_start $subject_end\n\n");
    }

    # write the new target exons
    $self->aceout("\n$biotype : $target\n");
    foreach ($subject_obj->Source_exons) {
      my ($start,$end) = $_->row(0);
      $self->aceout("Source_exons ".$start->name." ".$end->name."\n");
    }
    $self->aceout("Species \"$g_species\"\n");
    $self->aceout("CDS\n") if ($biotype eq 'CDS');

    # transfer Remarks, etc tags
    foreach ($subject_obj->Remark) {
      my ($remark, $evidence, $evidence_value1, $evidence_value2) = $_->row(0);
      my @evidence_type = $_->at('[1]');
      my @evidence_value = $_->at('[2]');
      foreach my $evidence_type (@evidence_type) {
	my $value = shift @evidence_value;
	$self->aceout("Remark \"". $remark->name ."\"");
	$self->aceout(" \"". $evidence_type->name."\"") if ($evidence_type);
	$self->aceout(" \"". $value->name."\"") if ($value);
	$self->aceout("\n");
      }
    }
    foreach ($subject_obj->Isoform) {
      my ($isoform, $evidence, $evidence_value1, $evidence_value2) = $_->row(0);
      $self->aceout("Isoform \"". $isoform->name ."\"");
      $self->aceout(" \"". $evidence->name."\"") if ($evidence);
      $self->aceout(" \"". $evidence_value1->name."\"") if ($evidence_value1);
      $self->aceout(" \"". $evidence_value2->name."\"") if ($evidence_value2);
      $self->aceout("\n");
    }
    foreach ($subject_obj->Evidence) { # for Michael
      my ($isoform, $evidence, $evidence_value1, $evidence_value2) = $_->row(0);
      $self->aceout("Evidence \"". $isoform->name ."\"");
      $self->aceout(" \"". $evidence->name."\"") if ($evidence);
      $self->aceout(" \"". $evidence_value1->name."\"") if ($evidence_value1);
      $self->aceout(" \"". $evidence_value2->name."\"") if ($evidence_value2);
      $self->aceout("\n");
    }
 
    $self->aceout("\n\n");

    $message .= "\n*** Fix: structure updated for $target\n "
  }

  # do we need to make any existing object into an isoform?
  if ($target =~ /(\S+?)b$/) { # making the 'b' isoform
    my $orig_structure = $1;
    my $orig_structure_obj;
    my $orig_structure_class;
    $orig_structure_obj = $self->{db}->fetch(CDS => "$orig_structure");
    $orig_structure_class = 'CDS';
    if (!defined $orig_structure_obj) {
      $orig_structure_obj = $self->{db}->fetch(Pseudogene => "$orig_structure");
      $orig_structure_class = 'Pseudogene';
      if (!defined $orig_structure_obj) {
	$orig_structure_obj = $self->{db}->fetch(Transcript => "$orig_structure");
	$orig_structure_class = 'Transcript';
      }
    }
    if (defined $orig_structure_obj) {
      my $personid = $self->person();
      $self->aceout("\n-R $orig_structure_class $orig_structure ${orig_structure}a\n");
      $self->aceout("\n");
      $self->aceout("\n$orig_structure_class : ${orig_structure}a\n");
      $self->aceout("Isoform Curator_confirmed $personid\n");
      $self->aceout("\n");
      $message .= "\n*** Fix: $orig_structure renamed to ${orig_structure}a\n";
    } else {
      $message .= "\n*** Well this is odd, can't find object '$orig_structure' to rename to the 'a' isoform!\n";
    }
  }

  print $message;

  if (!defined $gene_obj && !$target_exists) {
    print "\n*** Now set the Gene tag\n";
  }
  if ($target_exists) {
    print "\n*** Now make History for $biotype $target\n";
  }

}

###############################################################################
# this checks a few basic things about the structure of the given sequence

sub check {
  my ($self, $userinput) = @_;

  if (! $self->{interactive}) {return}

  my $status = 0;
  my $message = '';

  my @f = split /\s+/, $userinput;

  my $target = $f[1];
  if (! defined $target) {
    print "check what?\n";
    return $status;
  }

  # get the object
  my $target_obj = $self->{db}->fetch(CDS => "$target");
  if (! defined $target_obj) {
    $target_obj = $self->{db}->fetch(Transcript => "$target");
  }
  if (! defined $target_obj) {
    $message .= "Can't find '$target' - did you save it?\n";
    $status = 1;
  } else {

    # get the Method
    my $method = $target_obj->Method;
    if (!defined  $method) {
      $message .= "\n*** WARNING: No Method found\n";
      $status = 1;
    } else {
      $method = $method->name;
      if ($method ne "curated" && $method ne "Pseudogene" && $method ne "isoformer") {
	$message .= "\n*** WARNING: Invalid Method: '$method'\n";
	$status = 1;
      } else {
	if ($method ne "isoformer") { # warn about a Method="isoformer", but carry on checking
	  $message .= "\n*** WARNING: Invalid Method: '$method'\n"
	}
	# get the sorted exons
	my @exon_coord1 = sort by_number ($target_obj->get('Source_exons',1));
	my @exon_coord2 = sort by_number ($target_obj->get('Source_exons',2));
	my $target_hash_key = join(':', @exon_coord1) . ',' . join(':', @exon_coord2);
	print "$target exons: $target_hash_key\n";

	# get the Sequence start-end span 
	my $sequence_obj = $target_obj->Sequence;

	my @clone_locations;
	if ($method eq 'curated' || $method eq 'isoformer') {
	  @clone_locations = $sequence_obj->CDS_child;
	} else {
	  @clone_locations = $sequence_obj->Transcript;
	}    
	my ($target_start, $target_end);
	foreach my $target_location ( @clone_locations ) {
	  next unless ($target_location->name eq $target);
	  $target_start = $target_location->right->name;
	  $target_end = $target_location->right->right->name;
	  last;
	}

	# the exons should fit the Sequence span exactly
	if ($exon_coord2[$#exon_coord2] - $exon_coord1[0] != abs($target_start - $target_end)) {
	  $message .= "\n*** WARNING: The Source_exons do not fit in the span of this object on the Sequence (bogus extra exons?)\n";
	  $status = 1;
	} else {

	  # get any other objects with the same start-end
	  my ($test_name, $test_start, $test_end);
	  foreach my $test_location ( @clone_locations ) {
	    $test_name = $test_location->name;
	    $test_start = $test_location->right->name;
	    $test_end = $test_location->right->right->name;
	    if ($test_start == $target_start && $test_end == $target_end && $test_name ne $target) {
	      # there should not be any other objects with the same start-end and exon lengths
	      my @test_exon_coord1 = sort by_number ($test_location->follow->get('Source_exons',1));
	      my @test_exon_coord2 = sort by_number ($test_location->follow->get('Source_exons',2));
	      my $test_hash_key = join(':', @test_exon_coord1) . ',' . join(':', @test_exon_coord2);
	      my $test_name = $test_location->name;
	      print "$test_name exons: $test_hash_key\n";
	      if ($test_hash_key eq  $target_hash_key) {
		$message .= "\n*** WARNING: '$target' is the same structure as '$test_name'\n";
		$status = 1;
	      }
	    }
	  }	
	}
      }
    }
  }

  if ($status == 0) {
    $message .= "'$target' looks OK\n";
  }

  print $message;
  return $status;
}
###############################################################################
# command to display information about the TSLs in the region of interest
sub check_tsls {
  my ($self, $userinput) = @_;

  if (! $self->{interactive}) {return}

  my $status = 0;
  my $message = '';

  my @f = split /\s+/, $userinput;

  my $target = $f[1];
  if (! defined $target) {
    print "tsl what?\n";
    return $status;
  }

  # get the region of interest from the CDS name or clone positions
  my ($chromosome, $region_start, $region_end, $sense, $biotype, $gene) = $self->get_active_region($target);

  my @TSL = &get_TSL_in_region($chromosome, $region_start, $region_end, $sense, $self->{all_TSL}{$chromosome});
  foreach my $TSL (@TSL) {
    my $Feature = $TSL->{id};
    my $type = $TSL->{tsl};
    print "$Feature $type ";
    foreach my $FT (split /\s+/, $Feature) {
      
      my $FT_obj = $self->{db}->fetch(Feature => $FT);
      if (!defined $FT_obj) {
	print "$Feature $type - can't find this in the database!\n";
      } else {
	my @sequences = $FT_obj->Defined_by_sequence;
	my @papers = $FT_obj->Defined_by_paper;
	my @persons = $FT_obj->Defined_by_person;
	my @analyses = $FT_obj->Defined_by_analysis;
	my @remark = $FT_obj->Remark;
	my $result=0;
	my $report='';
	my $c = scalar @sequences;
	if ($c) {$result=1; $report.=" $c sequences"}
	$c = scalar @papers;
	if ($c) {$result=1; $report.=" $c papers"}
	$c = scalar @persons;
	if ($c) {$result=1; $report.=" $c persons"}
	$c = scalar @analyses;
	my $d = scalar @remark;
	my $e=0;
	if (defined $remark[0]) {
	  ($e) = ($remark[0] !~ /\s\d+\sreads$/);
	}
	if ($d < $e) {$d = $e}
	if ($c > 1) {$result=1; $report.=" $c analyses"} 
	elsif ($d > 1) {$result=1; $report.=" $d remarks"}
	my $r = ($result ? "Looks Good" : "Too weak");
	print " $r ($report)  ";
      }
    }
    print "\n";
  }

}
###############################################################################
# force a structure to become a functional CDS by adding in small
# introns to jump around premature STOP codons

sub force {

  my ($self, $userinput) = @_;

  if (! $self->{interactive}) {return}

  my $status = 0;
  my $message = '';
  my $end_not_found;
  my $start_not_found;

  my @f = split /\s+/, $userinput;

  my $target = $f[1];
  if (! defined $target) {
    print "force what?\n";
    return $status;
  }

  # get the object
  my $biotype = 'CDS';
  my $target_obj = $self->{db}->fetch(CDS => "$target");
  if (! defined $target_obj) {
    $biotype = 'Transcript';
    $target_obj = $self->{db}->fetch(Transcript => "$target");
  }
  if (! defined $target_obj) {
    $message .= "Can't find '$target' - did you save it?\n";
    $status = 1;
  } else {

    my $sense = '+';

    # get the sorted exons
    my @exon_coord1 = sort by_number ($target_obj->get('Source_exons',1));
    my @exon_coord2 = sort by_number ($target_obj->get('Source_exons',2));
    
    # get the Sequence start-end span 
    my $clone = $target_obj->Sequence;
    if (!defined $clone) {
      print STDERR "Sequence is not defined for the object $f[0]\n";
      return undef;
    }
    my @clone_CDSs;
    if ($biotype eq 'CDS') {
      @clone_CDSs = $clone->CDS_child;
    } elsif ($biotype eq 'Transcript') {
      @clone_CDSs = $clone->Transcript;
    }
    my ($clone_start, $clone_end);
    foreach my $CDS ( @clone_CDSs ) {
      next unless ($CDS->name eq "$target");
      $clone_start = $CDS->right->name;
      $clone_end = $CDS->right->right->name;
      last;
    }

    my @chrom_exons1; # get exons in chromosome forward sense coords
    my @chrom_exons2;
    my ($chromosome, $start, $end);
    my $coords = $self->{coords};
    ($chromosome, $start) = $coords->Coords_2chrom_coords($clone, $clone_start);
    ($chromosome, $end) = $coords->Coords_2chrom_coords($clone, $clone_end);
    if ($start > $end) {
      $sense = '-';
      foreach my $e (@exon_coord1) {push @chrom_exons1, $start - $e + 1}
      foreach my $e (@exon_coord2) {push @chrom_exons2, $start - $e + 1}
    } else { # forward sense
      foreach my $e (@exon_coord1) {push @chrom_exons1, $start + $e - 1}
      foreach my $e (@exon_coord2) {push @chrom_exons2, $start + $e - 1}
    }

    @chrom_exons1 = sort by_number @chrom_exons1;
    @chrom_exons2 = sort by_number @chrom_exons2;

    # make the sequence from the exons
    my $cds_seq;
    for (my $i=0; $i < scalar @chrom_exons1; $i++) {  
      my $exon_start = $chrom_exons1[$i];
      my $exon_end = $chrom_exons2[$i];
      my $len = $exon_end - $exon_start + 1;
      $cds_seq .= $self->{seq_obj}->Sub_sequence($chromosome, $exon_start-1, $len);
    }

    if ($sense eq '-') {
      $cds_seq = $self->{seq_obj}->DNA_revcomp($cds_seq);
    }

    # find the STOP codons, including the last codon
    my @stops;
    foreach my $stop ('taa', 'tga', 'tag') {
      my $st_offset = 0;
      my $st_result = index($cds_seq, $stop, $st_offset);
      while ($st_result != -1) {
	my $frame = $st_result % 3;
	if ($frame == 0) {
	  push @stops, $st_result; # 0 is the first base
	}
	$st_offset = $st_result + 1;
	$st_result = index($cds_seq, $stop, $st_offset);
      }
    }
    
    # check the last codon is at the end of the Sequence span of the structure
    my $last_stop = $stops[$#stops];
    if (!defined $last_stop || length $cds_seq != $last_stop - 3) {
      $message .= "No final STOP codon found - setting 'End_not_found' tag\n";
      $end_not_found = 1;
    } else { # remove the final STOP codon
      pop @stops;
    }
    if ($cds_seq !~ /^atg/) {
      $message .= "No START codon found - setting 'Start_not_found' tag\n";
      $start_not_found = 1;
    }

    # find each premature STOP codon
    # if it is less than three bases away from an intron splice site, then print warning and return
    # add an intron to splice around the STOP codon
    foreach my $stop (@stops) {
      for (my $i=0; $i < scalar @exon_coord1; $i++) { 
	my $done=0;
	if ($stop > $exon_coord1[$i] && $stop < $exon_coord2[$i] - 3) { # $stop is zero-base coord, $exon_coord is one-based
	  push @exon_coord1, $stop;
	  push @exon_coord2, $stop+4;
	  $done = 1;
	  last;
	}
	if (!$done) {$message .= "STOP number $i is too close to the exon edge to splice around\n";}
      }
    }
    @exon_coord1 = sort by_number @exon_coord1;
    @exon_coord2 = sort by_number @exon_coord2;


    # print the new Source_exons
    $self->aceclear();
    $self->aceout("\n$biotype : $target\n");
    $self->aceout("End_not_found\n") if (defined $end_not_found && $biotype eq 'CDS');
    $self->aceout("Start_not_found\n") if (defined $start_not_found && $biotype eq 'CDS');
    for (my $i=0; $i < scalar @exon_coord1; $i++) { 
      my ($exon1, $exon2) = ($exon_coord1[$i], $exon_coord2[$i]);
      $self->aceout("Source_exons ".$exon1." ".$exon2."\n");
    }
  }


  if ($status == 0) {
    $message .= "'$target' looks OK\n";
  }

  print $message;

}

###############################################################################
sub by_number{ $a <=> $b;}
sub by_number_rev{ $b <=> $a;}
###############################################################################
sub person {
  my ($self) = @_;
  my $personid;
  my $USER = $ENV{USER};
  if ($USER eq 'gw3') {
    $personid = 'WBPerson4025';
  } elsif ($USER eq 'pad') {
    $personid = 'WBPerson1983';
  } elsif ($USER eq 'mh6') {
    $personid = 'WBPerson4055';
  } elsif ($USER eq 'jl16') { # Jane Lomax
    $personid = 'WBPerson28994';
  } elsif ($USER eq 'jane') { # Jane Lomax
    $personid = 'WBPerson28994';
  } else {
    die "Unknown WBPerson id: $USER\n";
  }
  return $personid;
}
###############################################################################

1;

__END__

