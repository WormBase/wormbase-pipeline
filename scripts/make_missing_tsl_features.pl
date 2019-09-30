#!/software/bin/perl -w
#
# script to extract polyA and TSL information from the aligned
# EST/mRNA/Trinity/Nanopore Sequence objects and to create Feature objects
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
use Ace;
use Sequence_extract;
use Coords_converter;
use Feature_mapper;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $database, $new_feature_id, $output);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "database:s" => \$database,
	    "next_feature_id:s" => \$new_feature_id, # specify the first Feature ID to be used when making new Features, subsequent ones will be formed by incrementing this
	    "output:s"   => \$output, # output ACE file name
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

if (! defined $database) {$database = $wormbase->autoace}

if (! defined $new_feature_id) {$new_feature_id = 'WBsf990000';}
$log->write_to("First Feature ID to use set to '$new_feature_id'\n");

if (! defined $output) {$output = 'new_TSL_features.ace'}

#################################

my %mol_types = ( 'elegans'          => [qw( EST mRNA ncRNA OST tc1 RST Trinity Nanopore)],
                  'briggsae'         => [qw( mRNA EST Trinity)],
                  'remanei'          => [qw( mRNA EST Trinity)],
                  'brenneri'         => [qw( mRNA EST Trinity)],
                  'japonica'         => [qw( mRNA EST Trinity)],
                  'brugia'           => [qw( mRNA EST Trinity IsoSeq)],
                  'pristionchus'     => [qw( mRNA EST Trinity)],
                  'ovolvulus'        => [qw( mRNA EST Trinity)],
                  'sratti'           => [qw( mRNA EST Trinity)],
                  'tmuris'           => [qw( mRNA EST Trinity IsoSeq)],
		);

# valid Feature_data methods
my %valid_methods = (
                     "SL1"               => 1, # SL1 TSL data
                     "SL2"               => 1, # SL2 TSL data
                     "polyA_site"        => 1, # polyA tail sequences
		     "polyA_signal_sequence" => 1, # polyA signal sequence
                    );



# Set up top level base directories (these are different if in test mode)

my $chromosomes_dir = $wormbase->chromosomes; # AUTOACE CHROMSOMES
my $gff_dir         = $wormbase->gff;         # AUTOACE GFF
my $gff_splits_dir  = $wormbase->gff_splits;  # AUTOACE GFF SPLIT



##########################

my $species = $wormbase->species;
my $seq_obj = Sequence_extract->invoke($database, undef, $wormbase);
my $coords = Coords_converter->invoke($database, undef, $wormbase);
my $mapper = Feature_mapper->new($database, undef, $wormbase);
my %sl1_features = &get_genome_features('SL1');
my %sl2_features = &get_genome_features('SL2');
my %polya_features = &get_genome_features('polyA_site');
my %polya_signal_features = &get_genome_features('polyA_signal_sequence');
my %seq2feature = &fetch_features();
my %estorientation = $wormbase->FetchData('estorientation'); 

my $tsl_count = 0;
my $polya_count = 0;
my $tsl_false_count = 0;
my $polya_false_count = 0;
my $used_existing_features = 0;
my $made_new_features = 0;
my $non_splice_site_tsl = 0;
my $existing_tsl_features = 0;
my $new_tsl_features = 0;
my $existing_polya_features = 0;
my $new_polya_features = 0;

open (ACE, ">$output") || die "Can't open $output\n";

foreach my $mol_type (@{$mol_types{$species}}) {
  
  print "Processing $mol_type data\n";
  
  my ($first, $last) = &process_file($mol_type);
  
  foreach my $id (keys %{$first}) {
    my @feature_data = &match_sequence($id);
    if (scalar @feature_data) {
      my @position = &validated($id, 'first', $first->{$id}, \@feature_data);
      if (@position) {check_for_existing_Feature(\@position)}
    }
  }
  
  foreach my $id (keys %{$last}) {
    my @feature_data = &match_sequence($id);
    if (scalar @feature_data) {
      my @position = &validated($id, 'last', $last->{$id}, \@feature_data);
      if (@position) {check_for_existing_Feature(\@position)}
    }
  }
  
}

close (ACE);

$log->write_to("Found $non_splice_site_tsl TSL sites that do not look like a splice site but have been observed in other projects\n");
$log->write_to("$tsl_count TSL sites tested - false positive: $tsl_false_count (".$tsl_false_count*100/$tsl_count."%)\n");
$log->write_to("$polya_count POLYA sites tested - false positive: $polya_false_count (".$polya_false_count*100/$polya_count."%)\n");
$log->write_to("Used $used_existing_features existing Feature objects\n");
$log->write_to("Used $existing_tsl_features existing TSL Feature objects.\n");
$log->write_to("Used $existing_polya_features existing polyA Feature objects.\n");
$log->write_to("Made $new_tsl_features new TSL Feature objects.\n");
$log->write_to("Made $new_polya_features new polyA Feature objects.\n");
$log->write_to("Made $made_new_features new Feature objects. (The next unused Feature ID would now be: '$new_feature_id')\n");




$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################

############################################################
# is there a Feature there already?
sub check_for_existing_Feature {
  my ($new_feature) = @_;
  
  foreach my $new (@{$new_feature}) {
    my ($id, $method, $seq, $pos, $sense) = @{$new};
    
    if ($method eq 'polyA_site') {
      if (exists $polya_features{"$seq:$pos:$sense"}) {
	&add_to_feature($id, $polya_features{"$seq:$pos:$sense"}); 
	$existing_polya_features++;
      } else {
	&make_new_feature($id, $method, $seq, $pos, $sense)
      }
      
    } elsif ($method eq 'polyA_signal_sequence') {
      if (exists $polya_signal_features{"$seq:$pos:$sense"}) {
	&add_to_feature($id, $polya_signal_features{"$seq:$pos:$sense"})
      } else {
	&make_new_feature($id, $method, $seq, $pos, $sense)
      }
      
    } elsif ($method eq 'SL1') {
      if (exists $sl1_features{"$seq:$pos:$sense"}) {
	&add_to_feature($id, $sl1_features{"$seq:$pos:$sense"}); 
	$existing_tsl_features++;
      } else {
	&make_new_feature($id, $method, $seq, $pos, $sense)
      }
      
    } elsif ($method eq 'SL2') {
      if (exists $sl2_features{"$seq:$pos:$sense"}) {
	&add_to_feature($id, $sl2_features{"$seq:$pos:$sense"}); 
	$existing_tsl_features++;
      } else {
	&make_new_feature($id, $method, $seq, $pos, $sense)
      }
    }
  }
}


############################################################
# add Sequence to existing Feature

sub add_to_feature {
  my ($sequence_id, $feature) = @_;

  $used_existing_features++;

  my ($id, $seq, $start, $end, $sense) = @{$feature};

  print ACE "\n\n";
  print ACE "Feature : $id\n";
  print ACE "Defined_by_sequence \"$sequence_id\" Inferred_automatically \"make_missing_tsl_features.pl\"\n";

}

############################################################
# create new Feature

sub make_new_feature {
#  transcript id, Method of Feature_data, genomic sequence id, pos of start of Feature (always smaller than the end pos), sense + or -
  my ($id, $method, $seq, $pos, $sense) = @_;

  # note we have used the site now
  if ($method eq 'polyA_site') {
    $polya_features{"$seq:$pos:$sense"} = [$new_feature_id, $seq, $pos, 0, $sense];
    
  } elsif ($method eq 'polyA_signal_sequence') {
    $polya_signal_features{"$seq:$pos:$sense"} = [$new_feature_id, $seq, $pos, 0, $sense];
    
  } elsif ($method eq 'SL1') {
    $sl1_features{"$seq:$pos:$sense"} = [$new_feature_id, $seq, $pos, 0, $sense];
    
  } elsif ($method eq 'SL2') {
    $sl2_features{"$seq:$pos:$sense"} = [$new_feature_id, $seq, $pos, 0, $sense];
    
  }
  

  # make the new feature
  my $feature_start = $pos;
  my $feature_end = $pos+1;
  if ($method eq 'polyA_signal_sequence') {$feature_end = $pos + 5}

  # define a padding size to help get the flanking sequences entirely within one sequence object
  my $flank_size = 40;

  # get the clone
  my $prefix = $wormbase->chromosome_prefix();
  # get the clone for this Feature
  # add the flanking sequence length to the start and end to force LocateSpan to return a
  # large enough sequence object
  my $seq_flank_start = $feature_start - $flank_size; 
  my $seq_flank_end = $feature_end + $flank_size;
  if ($sense eq '-') { # reversed orientation on reverse strand for acedb coords
    $seq_flank_start = $feature_end + $flank_size;
    $seq_flank_end = $feature_start - $flank_size;
  }
  my @clone_coords = $coords->LocateSpan($seq, $seq_flank_start, $seq_flank_end ); # end is before start if reverse strand
  my $clone_name = $clone_coords[0];
  my $clone_start_pos = $clone_coords[1];       # the start position of the flanking sequences in clone coordinates
  my $clone_end_pos = $clone_coords[2];         # the end position of the flanking sequences in clone coordinates
  # now knock off the padding we added to ensure that LocateSpan returned a big enough sequence object
  if ($sense eq '+') {
    $clone_start_pos += $flank_size;
    $clone_end_pos -= $flank_size;
  } else {
    $clone_start_pos -= $flank_size;
    $clone_end_pos += $flank_size;
  }
  print "strand $sense clone_start_pos $clone_start_pos clone_end_pos $clone_end_pos\n" if ($verbose);

  # get the unique flanking sequences for this Feature
  my $is_zero_length = (abs($clone_end_pos - $clone_start_pos) == 1) ? 1 : 0; # the only non-zero length features here are the 6-base polyA signal sequences
  my $min_flank_length = 30;
  my $no_unique_check = 0;
  my $short_end_flanks_allowed = 0;
  my ($bases_upstream, $bases_downstream) = $mapper->get_flanking_sequence_for_feature($clone_name, $clone_start_pos, $clone_end_pos, $is_zero_length, $min_flank_length, $no_unique_check, $short_end_flanks_allowed);

#($clone_name, $clone_start_pos, $clone_end_pos);
  if (! defined $bases_upstream) {
    print "Couldn't find flanking sequences in clone $clone_name, trinig in $seq instead\n" if ($verbose);
    if ($sense eq '-') {($feature_start, $feature_end) = ($feature_end, $feature_start)}
    ($bases_upstream, $bases_downstream) = $mapper->get_flanking_sequence_for_feature($seq, $feature_start, $feature_end, $is_zero_length, $min_flank_length, $no_unique_check, $short_end_flanks_allowed);
    $clone_name = $seq;
    ($clone_start_pos, $clone_end_pos) = ($feature_start, $feature_end);
    if (! defined $bases_upstream) {
      die "******* ERROR ******* FLANKING SEQUENCE NOT FOUND FOR $id\n";
    }
  }
  

  # Write ace to define the Feature mapping data in the clone -
  # this will be done again by feature_mapper.pl in the Build, but it will help in
  # debugging to see the new features.
  
  my $full_species_name = $wormbase->full_name;
  
  print ACE "\n";
  print ACE "Sequence : $clone_name\n";
  # reversed if in reverse strand
  print ACE "Feature_object $new_feature_id $clone_start_pos $clone_end_pos\n";     
  
  # write ace to define this Feature
  print ACE "\n";
  print ACE "Feature : \"$new_feature_id\"\n";
  print ACE "Flanking_sequences \"$bases_upstream\" \"$bases_downstream\"\n";
  print ACE "Sequence \"$clone_name\"\n";
  print ACE "Mapping_target \"$clone_name\"\n";
  print ACE "Species \"$full_species_name\"\n";
  print ACE "Method \"$method\"\n"; # SL1 or SL2 or polyA_site or polyA_signal (if $method = polyA)
  print ACE "Defined_by_sequence \"$id\"\n";
  if ($method =~ /SL/) {
    print ACE "Description \"$method trans-splice leader acceptor site\"\n";
    print ACE "SO_term \"SO:0000706\"\n";
    $new_tsl_features++;
  } elsif ($method = "polyA_site") {
    print ACE "Description \"$method\"\n";
    print ACE "SO_term \"SO:0000553\"\n";
    $new_polya_features++;
  } elsif ($method = "polyA_signal_sequence") {
    print ACE "Description \"polyA signal sequence\"\n";
    print ACE "SO_term \"SO:0000551\"\n";
  } else {
    die "$method is not recognised in make_new_feature()\n";
  }
  print ACE "\n";
  
  $new_feature_id++; # get next new feature id
  
  $made_new_features++; # count of features made

}


############################################################
# validate that this is a reasonable position for the Feature
sub validated {
  my ($id, $place, $exon, $feature_data_aref) = @_;
  
  my @validated = ();
  

  my ($seq, $start, $end, $sense, $qstart, $qend) = @{$exon};

  if (! exists $estorientation{$id}) {return @validated}
  my $EST5 = ($estorientation{$id} eq '5');
  my $polyA_site_is_valid = 0; # note whether we have a valid polyA site - allow a polyA signal sequence only if there is a valid polyA site
  my @possible_polyA_signal_sequence = (); # temporary store for signal sequence

  foreach my $feature_data (@{$feature_data_aref}) {

    my ($method, $begin, $stop) = @{$feature_data};

    # is it a TSL at the start (first exon)?
    if ($feature_data->[0] eq 'SL1' || $feature_data->[0] eq 'SL2') {
      $tsl_count++;

      if ($EST5) { # forward read
	if ($place eq 'first') {
	  if ($stop+1 == $qstart) { # check for base after TSL being at start of first exon
	    my $pos; # start of exon aligned to chromosome (base after splice site)
	    if ($sense eq '+') {$pos = $start-1} else {$pos = $end} # so splice site is between $pos-1 and $pos in both + and - sense
	    # does the genomic sequence look like a splice site?
	    if (&check_TSL_site($id, $seq, $pos, $sense)) {    
	      push @validated, [$id, $method, $seq, $pos, $sense];
	    } else {
	      $tsl_false_count++;
	      print "$id $place exon = @{$exon}\n" if ($verbose);
	      print "feature_data = @{$feature_data }\n" if ($verbose);
	      print $seq_obj->Sub_sequence($seq, $pos-3, 7)."\n" if ($verbose);
	    }
	    
	  } else {
	    # ignore features that are not in this exon
	  }
	}
      } else { # EST3 - rev comp read
	if ($place eq 'last') {
	  if ($begin-1 == $qend) { # check for base before TSL being at end of last exon
	    my $pos; # end of exon aligned to chromosome (base after splice site)
	    if ($sense eq '+') {$sense='-'; $pos = $end} else {$sense='+'; $pos = $start-1} # so splice site is between $pos-1 and $pos in both + and - sense
	    # does the genomic sequence look like a splice site?
	    if (&check_TSL_site($id, $seq, $pos, $sense)) {
	      push @validated, [$id, $method, $seq, $pos, $sense];
	    } else {
	      $tsl_false_count++;
	      print "$id $place exon = @{$exon}\n" if ($verbose);
	      print "feature_data = @{$feature_data }\n" if ($verbose);
	      print $seq_obj->Sub_sequence($seq, $pos-3, 7)."\n" if ($verbose);
	    }
	    
	  } else {
	    # ignore features that are not in this exon
	  }
	}
      }
      
      # is it a polyA at the end (last exon)?
    } elsif ($feature_data->[0] eq 'polyA_site') {
      $polya_count++;

      if ($EST5) {
	if ($place eq 'last') {
	  if ($begin-1 == $qend || $begin == $qend) { # check for base before polyA being at end of last exon (allowing for 1 base error)
	    my $pos; # end of exon aligned to chromosome (base before polyA cut site)
	    if ($sense eq '+') {$pos = $end} else {$pos = $start-1} # so polyA cut site is between $pos-1 and $pos in both + and - sense
	    # does the genomic sequence look like a polyA site
	    if (&check_polyA_site($seq, $pos, $sense)) {
	      push @validated, [$id, $method, $seq, $pos, $sense];
	      $polyA_site_is_valid = 1;
	    } else {
	      $polya_false_count++;
	      print "$id $place exon = @{$exon}\n" if ($verbose);
	      print "feature_data = @{$feature_data }\n" if ($verbose);
	      print $seq_obj->Sub_sequence($seq, $pos-10, 20)."\n" if ($verbose);
	    }
	  } else {
	    # ignore features that are not in this exon
	  }
	}
      } else { # EST3 - rev comp read
	if ($place eq 'first') {
	  if ($begin+1 == $qstart || $begin == $qstart) { # check for base before polyA being at end of last exon (allowing for 1 base error)
	    my $pos; # start of exon aligned to chromosome (base after polyA cut site)
	    if ($sense eq '+') {$sense='-';$pos = $start-1} else {$sense='+'; $pos = $end} # so polyA cut site is between $pos-1 and $pos in both + and - sense
	    # does the genomic sequence look like a polyA site
	    if (&check_polyA_site($seq, $pos, $sense)) {    
	      push @validated, [$id, $method, $seq, $pos, $sense];
	      $polyA_site_is_valid = 1;
	    } else {
	      $polya_false_count++;
	      print "$id $place exon = @{$exon}\n" if ($verbose);
	      print "feature_data = @{$feature_data }\n" if ($verbose);
	      print $seq_obj->Sub_sequence($seq, $pos-10, 20)."\n" if ($verbose);
	    }
	  } else {
	    # ignore features that are not in this exon
	  }
	}
      }
          
      # is it a polyA_signal_sequence at the end (last exon)?
    } elsif ($feature_data->[0] eq 'polyA_signal_sequence') {
      
#      if ($EST5) {
#	if ($place eq 'last') {
#	  my $pos; # first base of signal sequence
## there is a problem with getting the genomeic sequence when the $sense = '-'
## so ignored the few signal sequence regions and don't make them Features
## so commented out this signal sequence bit
#	  if ($sense eq '+') {$pos = $start+$begin-$qstart} else {$pos = $start+$stop-$qend} # $pos is the first base of the 6-base signal
#	  if (&check_polyA_signal($seq, $pos, $sense)) {
#	    push @possible_polyA_signal_sequence, [$id, $method, $seq, $pos, $sense];
#	  } else {
#	    print "****************************** $place $sense non-canonical signal\n"
#	  }
#	  print "$id $place $sense exon = @{$exon}\n";
#	  print "feature_data = @{$feature_data }\n";
#	}
#      } else { # EST3 - rev comp read
#	if ($place eq 'first') {
#	  my $pos; # last base of signal sequence
#	  if ($sense eq '+') {$sense='-';$pos = $start+$stop-$qend} else {$sense='+'; $pos = $start+$begin-$qstart} # $pos is the first base of the 6-base signal
#	  if (&check_polyA_signal($seq, $pos, $sense)) {
#	    push @possible_polyA_signal_sequence, [$id, $method, $seq, $pos, $sense];
#	  } else {
#	    print "****************************** $place $sense non-canonical signal\n"
#	  }
#	  print "$id $place $sense exon = @{$exon}\n";
#	  print "feature_data = @{$feature_data }\n";
#	}
#      }
  }
    
  }

#  # if there is a valid polyA, accept any polyA signal sequence as being valid
#  foreach my $val (@validated) {
#    if ($val->[1] eq 'polyA_site') {
#      push @validated, @possible_polyA_signal_sequence;
#      last;
#    }
#  }
  
  return @validated;
}


###############################################################
# see if these is a AG canonical splice site at the TSL splice site
sub check_TSL_site {

  my ($id, $seq, $pos, $sense) = @_;
  my $result = 0; # 0 is failed,    1 is OK

  my $zero_pos = $pos-1; # convert to zero-index coords for Sequence_extract

  my $sl1 = $sl1_features{"$seq:$pos:$sense"};
  my $sl2 = $sl2_features{"$seq:$pos:$sense"};
  if (!defined $sl1) {$sl1 = ''} else {$sl1 = $sl1->[0]}
  if (!defined $sl2) {$sl2 = ''} else {$sl2 = $sl2->[0]}
  
  if ($sense eq '+') {
    my $sequence = $seq_obj->Sub_sequence($seq, $zero_pos-1, 2);
    if ($sequence eq 'ag') {
      $result = 1; # OK
    } else {
      if (exists $sl1_features{"$seq:$pos:$sense"} || exists $sl2_features{"$seq:$pos:$sense"}) { # we accept this as a TSL sites if a TSL has been observed here before!
	$result = 1; # I'm fairly happy to accept these
	$non_splice_site_tsl++;
	print "*************** non-canonical TSL site with existingFeature: $id, $sl1 $sl2, $seq, $pos, $sense ".$seq_obj->Sub_sequence($seq, $zero_pos-20, 25)."\n" if ($verbose);

      } else { # non-canonical slice site and no Features made here before
	$result = 0;  # these are highly suspect - reject them
	print "*************** weird un-confirmed TSL site: $id, $sl1 $sl2, $seq, $pos, $sense ".$seq_obj->Sub_sequence($seq, $zero_pos-20, 25)."\n" if ($verbose);

      }
    }
  } else {
    my $sequence = $seq_obj->Sub_sequence($seq, $zero_pos+1, 2);
    if ($sequence eq 'ct') {
      $result = 1; # OK
    } else {
      if (exists $sl1_features{"$seq:$pos:$sense"} || exists $sl2_features{"$seq:$pos:$sense"}) { # we accept this as a TSL sites if a TSL has been observed here before!
	$result = 1; # I'm fairly happy to accept these
	$non_splice_site_tsl++;
	print "*************** non-canonical TSL site with existing Feature: $id, $sl1 $sl2, $seq, $pos, $sense ".$seq_obj->DNA_revcomp($seq_obj->Sub_sequence($seq, $zero_pos-5, 25))."\n" if ($verbose);

      } else {
	$result = 0; # these are highly suspect - reject them
	print "*************** weird un-confirmed TSL site: $id, $sl1 $sl2, $seq, $pos, $sense ".$seq_obj->DNA_revcomp($seq_obj->Sub_sequence($seq, $zero_pos-5, 25))."\n" if ($verbose);

      }
    }
  }


  return $result;
}


###############################################################
# see if there is 'aataaa'

sub check_polyA_signal {

  my ($seq, $pos, $sense) = @_;
  my $result = 0; # 0 is failed, 1 is OK

  $pos--; # convert to zero-index coords for Sequence_extract

  if ($sense eq '+') {
    my $sequence = $seq_obj->Sub_sequence($seq, $pos, 6);
    if ($sequence eq 'aataaa') {$result = 1} # OK
    print "$sequence\n";
    print $seq_obj->Sub_sequence($seq, $pos-6, 6)." $sequence ".$seq_obj->Sub_sequence($seq, $pos+6, 6)."\n";
  } else {
    my $sequence = $seq_obj->Sub_sequence($seq, $pos-6, 6);
    if ($sequence eq 'tttatt') {$result = 1} # OK
    print "$sequence\n";
    print $seq_obj->Sub_sequence($seq, $pos-2, 6)." $sequence ".$seq_obj->Sub_sequence($seq, $pos, 6)."\n";
  }


  return $result;
}


###############################################################
# see if these are 5 bases of polyA in the genome at the polyA cut
# site (there is a 1% change that this will occur even if this is a
# real cut site, but we can live with that)
sub check_polyA_site {

  my ($seq, $pos, $sense) = @_;
  my $result = 0; # 0 is failed, 1 is OK

  $pos--; # convert to zero-index coords for Sequence_extract

  if ($sense eq '+') {
    my $sequence = $seq_obj->Sub_sequence($seq, $pos+1, 5);
    if ($sequence ne 'aaaaa') {$result = 1} # OK
  } else {
    my $sequence = $seq_obj->Sub_sequence($seq, $pos-4, 5);
    if ($sequence ne 'ttttt') {$result = 1} # OK
  }


  return $result;
}


###############################################################
# get the data from a GFF file of BLAT alignments
sub process_file {
  my ($mol_type) = @_;

  print "Loading BLAT alignment data\n";

  my %first; # first exon, key $id
  my %last;  # last exon, key $id

  foreach my $file (glob "${gff_splits_dir}/*BLAT_${mol_type}_BEST.gff") {
    if (open (BLAT, "<$file")) {
      while (my $line = <BLAT>) {
	next if ($line =~ /^#/);
	my ($seq, $type, $method, $start, $end, $score, $sense, $phase, $other) = split /\t/, $line;
	next if ($type eq 'Link');
	my ($target, $id, $qstart, $qend) = split /\s+/, $other;
	$id =~ s/"//g;
	$id =~ s/Sequence://;
	if ($sense eq '+' && ! exists $first{$id}) {$first{$id} = [$seq, $start, $end, $sense, $qstart, $qend]};
	if ($sense eq '-' && ! exists $last{$id}) {$last{$id} = [$seq, $start, $end, $sense, $qstart, $qend]};
	if ($sense eq '+') {$last{$id} = [$seq, $start, $end, $sense, $qstart, $qend]};
	if ($sense eq '-') {$first{$id} = [$seq, $start, $end, $sense, $qstart, $qend]};
      }
      close (BLAT);
    } else {
      print "Can't open $file\n";
    }
  }

  return (\%first, \%last);
}


###############################################################
# get the Feature_data in all EST/OST/mRNA/RST/Trinity/Nanopore Sequences
sub fetch_features {

  print "Loading Sequence Feature_data\n";

  my %seq2feature;
  my $mismatched;
  my $tm_data = $wormbase->table_maker_query($database, $wormbase->basedir."/wquery/SCRIPT:transcriptmasker.def");
  while( <$tm_data> ) {
    s/\"//g;  #"
    next if (/acedb/ or /\/\// or /^$/);
    my ($seq, $type, $start, $end) = split;
    if ($seq and $type and $start and $end) {
      $seq2feature{$seq}->{$type} = [($start, $end)];
    }
    else {
      $log->error("something wrong with $seq:$type\n");
      $mismatched++;
    }
  }
  $log->write_to( scalar(keys %seq2feature) ." feature_data found OK, $mismatched feature_data malformed (possibly stub objects from another species?)\n");
  return %seq2feature;
}

###############################################################
# get the Feature from the GFF files
sub get_genome_features {

  my ($type) = @_;

  print "Loading $type genome Features\n";

  my %features;
  
  foreach my $file (glob "${gff_splits_dir}/*${type}.gff") {
    if (open (FT, "<$file")) {
      while (my $line = <FT>) {
	next if ($line =~ /^#/);
	my ($seq, $type, $method, $start, $end, $score, $sense, $phase, $other) = split /\t/, $line;
	next if ($type eq 'Link');
	my ($target, $id) = split /\s+/, $other;
	$id =~ s/"//g;
	$features{"$seq:$start:$sense"} = [$id, $seq, $start, $end, $sense];
      }
      close (FT);
    } else {
      print "Can't open $file\n";
    }
  }

  return %features;
}

###############################################################
# 
# get an array of the valid Feature_data for this $id

sub match_sequence {

  my $id = shift;

  my @results;

  if( $seq2feature{$id} ) {
    foreach my $method (keys %{$seq2feature{$id}}) {
      my $coords = $seq2feature{$id}->{$method};
      if (defined($method)) {
	next unless ($valid_methods{$method}); # only match valid Feature_data methods
	my ($start, $stop);
	$start = $coords->[0]; # start coord
	$stop  = $coords->[1]; # stop coord

	# don't want to use dubious polyA regions that could be part of a low-complexity region
	if ((exists $seq2feature{$id}->{'low'} || exists $seq2feature{$id}->{'low_complexity'}) && 
	    ($method eq 'polyA_site' || $method eq 'polyA_signal_sequence')) {next} 

	push @results, [$method, $start, $stop];
      }
    }
  }
  return @results;
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
