#!/software/bin/perl -w
#
# make_missing_tsl_features.pl
# 
# by Gary Williams                        
#
# This is a small script to look for EST sequences with TSL Feature_data
# and to create Features if they do not already exist.
#
# In the acedb main window get all RST* Sequences
# Then Query:   Feature_data="*:TSL"   gives 1106 out of 2944
# Then Query:   !Defines_feature   gives 1106 exported this list of names
#
# Read in the list of names from file tsl_not_defines_feature.list
# and create the Features, if necessary, or re-use existing ones.
#
# - read in all EST GFF and orientation data as usual
# - read in all SL1/SL2 GFF data
# - read in list of ESTs to process
# - for each EST from the list
#   - get type of TSL: SL1/SL2
#   - get position of TSL in the sequence
#   - map SL1/2 site to the genomic chromosome position
#   - see if we have a Feature of the right type there already
#   - if so use it, if not make a new feature
#
# N.B. This has been tested and used for RST* SL1/2 features, but not on anythng else
#
#


# Last updated by: $Author: pad $     
# Last updated on: $Date: 2015-07-03 09:51:59 $      
use lib '/nfs/WWWdev/SANGER_docs/lib/Projects/C_elegans';
use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
use Ace::Sequence;
use Sequence_extract;
use Coords_converter;
use Data::Dumper;
use Feature_mapper;
use Modules::Overlap;
use NameDB_handler;
use NameDB;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $USER,);
my ($input, $output, $database);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "input:s"    => \$input, # input list of sequences to work on
	    "output:s"    => \$output, #output ace file
	    "database:s" => \$database,	# defaults to currentdb
	    "user:s"     => \$USER,
	   );

$test = 1;

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
  print "In test mode\n";
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);



#################################
# Set up some useful paths      #
#################################


$database = $wormbase->database('current') unless $database;
my $tace            = $wormbase->tace;        # TACE PATH



##########################
# MAIN BODY OF SCRIPT
##########################

my $FEATURE = "TSL";			# type of feature we are doing TSL/polyAs
#my $FEATURE = "polyA_signal";			# type of feature we are doing TSL/polyAs
#my $FEATURE = "polyA_site";			# type of feature we are doing TSL/polyAs



my $seq_extract = Sequence_extract->invoke($database, undef, $wormbase);

print "Connecting to Ace\n";
my $db = Ace->connect (-path => $database,
                       -program => $tace) || die "cannot connect to database at $wormbase->database('current')\n";

# start using the Coords_converted module
my $coords = Coords_converter->invoke($database, undef, $wormbase);

# start using feature mapper
my $mapper = Feature_mapper->new($database, undef, $wormbase);

# get the Overlap object
my $ovlp = Overlap->new($database, $wormbase);

#connect to name server and set domain to 'Gene'
my $DB      = 'wbgene_id;mcs12a;3307';
my $DOMAIN  = 'Feature';
my $feature_db = NameDB_handler->new($DB,$USER,$USER,"/nfs/WWWdev/SANGER_docs/data");
    $db->setDomain('Feature');


# look in ~wormpub/DATABASES/current_DB/COMMON_DATA/clonesize.dat
# to get the size of the specified clone
my %clonesize;
$wormbase->FetchData("clonesize", \%clonesize, "$database/COMMON_DATA/");
# and get the centre that sequenced this clone
my %clonelab;
$wormbase->FetchData("clone2centre", \%clonelab, "$database/COMMON_DATA/");
my %Show_in_reverse_orientation;
$wormbase->FetchData("estorientation", \%Show_in_reverse_orientation, "$database/COMMON_DATA/");

print "read list of ESTs to work on\n";
my @est_list = &read_est_list($input);

open (OUT, "> $output") || die "Cant open file $output";

my @chromosomes = $wormbase->get_chromosome_names(-mito => 0, -prefix => 0);

foreach my $chromosome (@chromosomes) {

  $log->write_to("Processing chromosome $chromosome\n");

  print "get ESTs\n";
  my @est = $ovlp->get_EST_BEST($chromosome);

  print "get mRNAs\n";
  my @mRNA = $ovlp->get_mRNA_BEST($chromosome);

  print "get OSTs\n";
  my @ost = $ovlp->get_OST_BEST($chromosome);

  print "get RSTs\n";
  my @rst = $ovlp->get_RST_BEST($chromosome);

  print "get SL1\n";
  my @sl1_features = $ovlp->get_TSL_SL1($chromosome);

  print "get SL2\n";
  my @sl2_features = $ovlp->get_TSL_SL2($chromosome);

  print "get polyA_signal\n";
  my @polyA_signal_features = $ovlp->get_polyA_signal($chromosome);

  print "get polyA_site\n";
  my @polyA_site_features = $ovlp->get_polyA_site($chromosome);

  print "look for ESTs with missing TSL Features\n";
  &get_missing_features(\@est_list, \@est, \@mRNA, \@ost, \@rst, \@sl1_features, \@sl2_features, \@polyA_signal_features, \@polyA_site_features, $chromosome);

}

close (OUT);

$log->mail();

#print "Finished.\n" if ($verbose);
exit(0);


##############################################################
#
# Subroutines
#
##############################################################



##########################################
# do all the work here
sub get_missing_features {

  my ($est_list_aref, $est_aref, $mrna_aref, $ost_aref, $rst_aref, $sl1_aref, $sl2_aref, $signal_aref, $site_aref, $chromosome) = @_;

  my @est_list = @{$est_list_aref};
  my @est = @{$est_aref};
  my @mrna = @{$mrna_aref};
  my @ost = @{$ost_aref};
  my @rst = @{$rst_aref};

  my @STORE_FEATURE;		# results store for the TSL/polya mask


  # debug EST to trace
  my $trace; # = "D12984";
  if (defined $trace) {$verbose = 1;}


  # get the next EST hit
  foreach my $seq (@est_list) {
    #print "Next Sequence: $seq\n";

    my $feature_obj;
    my $feature_type;
    my $feature_start;
    my $feature_end;
    my $feature_len;

    # see if we can find this sequence in the lists of est,mrna,ost,rst
    $feature_obj = $db->fetch("Feature_data", "$seq:$FEATURE");
    $feature_type = $feature_obj->at("Feature[1]"); # SL1 or SL2[a-j]
    $feature_start = $feature_obj->at("Feature.$feature_type" . "[1]");	# polyA_signal always has type 'polyA' here for some reason
    $feature_end = $feature_obj->at("Feature.$feature_type" . "[2]");
    $feature_len = $feature_obj->at("Feature.$feature_type" . "[3]");
#    print "FEATURE data for $seq: $feature_type $feature_start $feature_end $feature_len\n";

    # slow and inefficient search for this Sequence in the lists of GFF data
    # find the HSP of this EST that contains the feature
    my $this_est_hsp = &get_HSP_with_feature($seq, $feature_start, $feature_end, \@rst, \@mrna, \@ost, \@est);

    if (!defined $this_est_hsp) {next;}	# the sequence does not align to this chromosome

    # $est_id, $chrom_start, $chrom_end, $chrom_strand, $hit_start, $hit_end, $score

    my $est_id = $this_est_hsp->[0];
    print "matches $est_id\n";

    my $est_start = $this_est_hsp->[1];	# chrom coords
    my $est_end = $this_est_hsp->[2];	# chrom coords always start < end
    my $est_strand = $this_est_hsp->[3];	# + or -
    my $est_pos_start = $this_est_hsp->[4]; # hit start pos, start < end if forwards
    my $est_pos_end = $this_est_hsp->[5]; # hit end pos, start > end if reverse

    # define direction of alignment to genome
    my $FORWARD = $est_pos_start < $est_pos_end;
    my $REVERSE = ! $FORWARD;

    # define direction of EST read
    my $orientation = $Show_in_reverse_orientation{$est_id}; 

    # get the EST DNA sequence and its length
    my $est_obj = $db->fetch("Sequence" => $est_id);
    my $est_seq = $est_obj->asDNA;
    $est_seq =~ s/^\n>.+//;
    $est_seq =~ s/\n//g;
    my $est_len = length $est_seq;

    # 
    # store the EST's TSL/polyA region details
    # 
    push @STORE_FEATURE, [($est_id, $feature_start, $feature_end, $feature_len, $est_len, $orientation, $this_est_hsp, $feature_type)];

  }				# foreach EST




  #
  # foreach EST we have found with a TSL/polyA region:
  #
  foreach my $details (@STORE_FEATURE) {

  # print the results
    my @details = @{$details};
    print "EST details @details\n"; # debug
    
    my $est_id = $details[0];	# EST ID
    my $est_feature_type = $details[7];	# SL1 or SL2 or 'polyA' (when $FEATURE eq 'polyA_signal') or 'polyA_site'

    # get the genomic location of the feature from the EST alignment
    my ($feature_start, $feature_end) = &get_feature_location($details);

    # see if there is a feature there already with the right type SL1/2
    my $est_gff = $details->[6];
    my $strand = $est_gff->[3];
    my $feature_id;
    my $found_a_match = 0;
    print "Searching for $est_feature_type\n";
    if ($est_feature_type eq "SL1") {
      foreach my $tsl (@{$sl1_aref}) {
	if ($feature_start == $tsl->[1] && $feature_end == $tsl->[2] && $strand eq $tsl->[3]) {
	  print "FOUND A MATCHING FEATURE $est_feature_type\n";
	  $feature_id = $tsl->[0];
	  print "EST strand = $strand, Feature $feature_id strand: $tsl->[3]\n";
	  $found_a_match = 1;
	  last;
	}
      }
    } elsif ($est_feature_type =~ /SL2/) {
      foreach my $tsl (@{$sl2_aref}) {
	if ($feature_start == $tsl->[1] && $feature_end == $tsl->[2] && $strand eq $tsl->[3]) {
	  print "*** FOUND A MATCHING FEATURE $est_feature_type\n";
	  $feature_id = $tsl->[0];
	  print "EST strand = $strand, Feature $feature_id strand: $tsl->[3]\n";
	  $found_a_match = 1;
	  last;
	}
      }
    } elsif ($est_feature_type eq "polyA") { # this is really "polyA_signal_sequence"
      foreach my $signal (@{$signal_aref}) {
	if ($feature_start == $signal->[1] && $feature_end == $signal->[2]) { # && $strand eq $signal->[3]) {
	  print "*** FOUND A MATCHING FEATURE $est_feature_type\n";
	  $feature_id = $signal->[0];
	  print "EST strand = $strand, Feature $feature_id strand: $signal->[3]\n";
	  $found_a_match = 1;
	  last;
	}
      }
    } elsif ($est_feature_type eq "polyA_site") {
      # want to find sites within about 3 bases
      foreach my $site (@{$site_aref}) {
#	if ($feature_start == $site->[1] && $feature_end == $site->[2]) { # && $strand eq $site->[3]) {
	if ($feature_start-3 <= $site->[1] && $feature_end+3 >= $site->[2]) {
	  print "*** FOUND A MATCHING FEATURE $est_feature_type\n";
	  $feature_id = $site->[0];
	  print "EST strand = $strand, Feature $feature_id strand: $site->[3]\n";
	  $found_a_match = 1;
	  last;
	}
      }
    } else {
      die "unknown feature type: $est_feature_type in EST $est_id\n";
    }

    # output ACE data to create the feature if it is not there
    if (!$found_a_match) {
      $feature_id = &create_new_feature($chromosome, $est_feature_type, $feature_start, $feature_end, $strand);

      # Now push the details of the newly-created Feature position
      # onto the feature GFF list so we find it on the next search at
      # this location so we don't make multiple new Features at the
      # same place.
      if ($est_feature_type eq "SL1") {
	push @{$sl1_aref}, [($feature_id, $feature_start, $feature_end, $strand)];
      } elsif ($est_feature_type =~ /SL2/) {
	push @{$sl2_aref}, [($feature_id, $feature_start, $feature_end, $strand)];
      } elsif ($est_feature_type eq "polyA") { # polyA_signal
	push @{$signal_aref}, [($feature_id, $feature_start, $feature_end, $strand)];
      } elsif ($est_feature_type eq "polyA_site") {
	push @{$site_aref}, [($feature_id, $feature_start, $feature_end, $strand)];
      } else {
	die "Unknown feature type: $est_feature_type\n";
      }
    }

    # output ACE data to link the EST to the feature
    print OUT "\n";
    print OUT "Feature : \"$feature_id\"\n";
    print OUT "Defined_by_sequence \"$est_id\" Inferred_automatically \"make_missing_tsl_features.pl\"\n"; # use this script's name

  }
}

##########################################
# output the ACE for a new feature
sub create_new_feature {
  my ($chromosome, $type, $chrom_start, $chrom_end, $strand) = @_;

  # Define a new feature ID to use
  my $feature_id = $feature_db->idCreate ;

  # define a padding size to help get the flanking sequences entirely within one sequence object
  my $flank_size = 100;

  # get the clone
  my $prefix = $wormbase->chromosome_prefix();
  my $chrom = $prefix . $chromosome;
  # get the clone for this Feature
  # add the flanking sequence length to the start and end to force LocateSpan to return a
  # large enough sequence object
  my $chrom_flank_start = $chrom_start - $flank_size; 
  my $chrom_flank_end = $chrom_end + $flank_size;
  if ($strand eq '-') { # reversed orientation on reverse strand for acedb coords
    $chrom_flank_start = $chrom_end + $flank_size;
    $chrom_flank_end = $chrom_start - $flank_size;
  }
  my @clone_coords = $coords->LocateSpan($chrom, $chrom_flank_start, $chrom_flank_end ); # end is before start if reverse strand
  my $clone_name = $clone_coords[0];
  my $clone_start_pos = $clone_coords[1]; 	# the start position of the flanking sequences in clone coordinates
  my $clone_end_pos = $clone_coords[2]; 	# the end position of the flanking sequences in clone coordinates
  # now knock off the padding we added to ensure that LocateSpan returned a big enough sequence object
  if ($strand eq '+') {
    $clone_start_pos += $flank_size;
    $clone_end_pos -= $flank_size;
  } else {
    $clone_start_pos -= $flank_size;
    $clone_end_pos += $flank_size;
  }
  print "strand $strand clone_start_pos $clone_start_pos clone_end_pos $clone_end_pos\n" if ($verbose);

  # get the unique flanking sequences for this Feature
  my ($bases_upstream, $bases_downstream) = $mapper->get_flanking_sequence($clone_name, $clone_start_pos, $clone_end_pos);
  if (! defined $bases_upstream) {print "********************** FLANKING SEQUENCE NOT FOUND\n";}

  # Write ace to define the Feature mapping data in the clone -
  # this will be done again by feature_mapper.pl in the Build, but it will help in
  # debugging to see the new features.
      
  print OUT "\n";
  print OUT "Sequence : $clone_name\n";
  # reversed if in reverse strand
  print OUT "Feature_object $feature_id $clone_start_pos $clone_end_pos\n";	

  # write ace to define this Feature
  print OUT "\n";
  print OUT "Feature : \"$feature_id\"\n";
  print OUT "Flanking_sequences \"$clone_name\" \"$bases_upstream\" \"$bases_downstream\"\n";
  print OUT "Species \"Caenorhabditis elegans\"\n";
  print OUT "Method \"$type\"\n"; # SL1 or SL2 or polyA_site or polyA_signal (if $type = polyA)
  if ($type =~ /SL/) {
    print OUT "Description \"$type trans-splice acceptor\"\n";
  } else {
    print OUT "Description \"$type\"\n";
  }
  print OUT "\n";

  return $feature_id;
}

##########################################
# # get the genomic location of the feature from the EST alignment
sub get_feature_location {
  my ($details_aref) = @_;

  my $est = $details_aref->[6];	# the EST GFF line
  my $est_feature_end = $details_aref->[2]; # end of feature on EST
  my $est_genomic_start = $est->[1];
  my $est_genomic_end   = $est->[2];
  my $est_hsp_start = $est->[4]; # start of HSP alignment on EST
  my $est_hsp_end   = $est->[5]; #   end of HSP alignment on EST
  my $orientation = $details_aref->[5];

  my ($genomic_feature_start, $genomic_feature_end);		# resulting feature genommic location
  my $diff;

  if ($FEATURE eq "TSL") {			# TSL/polyA feature
    if ($orientation == 5) {	# 5' orientation read
      if ($est_hsp_start < $est_hsp_end) {
	$diff = $est_feature_end - $est_hsp_start + 1;
	$genomic_feature_end = $est_genomic_start + $diff; # genomic location of end of feature
	$genomic_feature_start = $genomic_feature_end - 1; # TSL features are only 2 bases long
      } else {
	$diff = $est_feature_end - $est_hsp_end + 1;
	$genomic_feature_start = $est_genomic_end - $diff; # genomic location of start of feature
	$genomic_feature_end = $genomic_feature_start + 1; # TSL features are only 2 bases long
      }
    } else {			# there are no 3' RSTs with TSL features, so ignore this bit
      # +++
      
    }

  } elsif ($FEATURE eq "polyA_signal") { # polyA signal feature 
    if ($orientation == 3) {	# there are no 5' RSTs with polyA features, so ignore this bit
      # +++

    } else {			# 3' orientation read
      if ($est_hsp_start < $est_hsp_end) {
	$diff = $est_feature_end - $est_hsp_end + 1;
	$genomic_feature_start = $est_genomic_end - $diff; # genomic location of start of feature
	$genomic_feature_end = $genomic_feature_start + 5; # signal features are only 2 bases long
      } else {
	$diff = $est_feature_end - $est_hsp_start + 1;
	$genomic_feature_end = $est_genomic_start + $diff; # genomic location of end of feature
	$genomic_feature_start = $genomic_feature_end - 5; # signal features are only 6 bases long
      }
    }

  } elsif ($FEATURE eq "polyA_site") { # polyA site feature 
    if ($orientation == 3) {	# there are no 5' RSTs with polyA features, so ignore this bit
      # +++

    } else {			# 3' orientation read
      if ($est_hsp_start < $est_hsp_end) {
	$diff = $est_feature_end - $est_hsp_end + 1;
	$genomic_feature_start = $est_genomic_end - $diff; # genomic location of start of feature
	$genomic_feature_end = $genomic_feature_start + 1; 
      } else {
	$diff = $est_feature_end - $est_hsp_start + 1;
	$genomic_feature_end = $est_genomic_start + $diff; # genomic location of end of feature
	$genomic_feature_start = $genomic_feature_end - 1; 
      }
    }


  }

  return ($genomic_feature_start, $genomic_feature_end);
}

##########################################
# find the HSP of this EST that contains the feature
# slow and inefficient search for this Sequence in the lists of GFF data
#    my $est = &get_HSP_with_feature($seq, \@rst, \@mrna, \@ost, \@est);
sub get_HSP_with_feature {

  my ($seq, $est_feature_start, $est_feature_end, $rst_aref, $mrna_aref, $ost_aref, $est_aref) = @_;

  my $est;
  $est = &get_HSP($seq, $est_feature_start, $est_feature_end, $rst_aref);
  $est = &get_HSP($seq, $est_feature_start, $est_feature_end, $mrna_aref) unless defined $est;
  $est = &get_HSP($seq, $est_feature_start, $est_feature_end, $ost_aref) unless defined $est;
  $est = &get_HSP($seq, $est_feature_start, $est_feature_end, $est_aref) unless defined $est;

  return $est;
}


##########################################
# search a GFF list of HSPs for a match to a sequence_id and a region of the EST
sub get_HSP {
  my ($seq, $feature_start, $feature_end, $gff_aref) = @_;

  my $est;

  if ($FEATURE eq "TSL") {			# TSL regions
    foreach my $g (@{$gff_aref}) {
      if ($seq eq $g->[0]) {
	print "Found an HSP of our sequence: $seq feature start-end: $feature_start, $feature_end HSP start end: ", $g->[4], " ", $g->[5], "\n";
	if (# TSL region is base before this hit and the first base of the hit

	    # we want the end of the TSL to be one base before the
	    # start of the HSP ($g->[4]) but allow for a few bases
	    # missing from the start of the alignment
	    # (e.g. RST5_373512, RST5_372643)
	    (($feature_end >= $g->[4]-5 && $feature_end < $g->[4]) ||
	     ($feature_end >= $g->[5]-5 && $feature_end < $g->[5]))) {
	  print "... and a match of the region\n";
	  $est = $g;
	  last;
	}
      }
    }

  } else {			# polyA region
    foreach my $g (@{$gff_aref}) {
      if ($seq eq $g->[0]) {
	print "Found an HSP of our sequence: $seq feature start-end: $feature_start, $feature_end HSP start end: ", $g->[4], " ", $g->[5], "\n";
	if (# polyA region is within this hit
	    (($feature_start >= $g->[4] && $feature_start <= $g->[5]) ||
	     ($feature_end >= $g->[4] && $feature_end <= $g->[5]))) {
	  $est = $g;
	  last;
	}
      }
    }
  }

  return $est;
}



##########################################
# read list of ESTs to work on
#my @est_list = &read_est_list();

sub read_est_list {
  my ($input) = @_;
  my @list;

  open (LIST, "<$input") || die "Can't open input file $input\n";
  while (my $line = <LIST>) {
    if ($line =~ /Sequence : \"(\S+)\"/) {
      push @list, $1;
    }
  }
  close (LIST);
  
  return @list;

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


##########################################
# block of old code that we need but want to put into subroutines

sub old_code {


    # 
    # See if this matches the location of an existing TSL/polyA-site Feature
    # 
    my $got_a_feature_site_match = "";
    #
    # Get the chromosomal location of the EST's TSL/polyA site
    # 
    my $this_est_start = $est->[1]; # pos on chromosome
    my $this_est_end = $est->[2];	# chrom coords always start < end
    my $this_est_strand = $est->[3];	# + or -
    my $this_est_pos_start = $est->[4]; # hit start pos, start < end if forwards
    my $this_est_pos_end = $est->[5]; # hit end pos, start > end if reverse
    # define direction of alignment to genome
    my $FORWARD = $this_est_pos_start < $this_est_pos_end;
    my $REVERSE = ! $FORWARD;
    # test if the location is past the end of the 3'-most hit so not mappable
    if ($FORWARD && $ORIENTATION && $feature_site > $this_est_pos_end ||
	$FORWARD && ! $ORIENTATION && $feature_site < $this_est_pos_start ||
	$REVERSE && $ORIENTATION && $feature_site > $this_est_pos_start ||
	$REVERSE && ! $ORIENTATION && $feature_site < $this_est_pos_end) {
      print "$est_id is not mappable - after the 3'-most hit\n";
      next;
    }

      # Here if the site is in the 3'-most EST hit.  We don't want to
      # define a Feature if the 3'-most EST hit was for a small A-rich
      # region which wouldn't be aligned to part of the CDS
      if ($type_of_hit eq "est hit") {
	next;
      }
    }
    my $diff_pos;
    print "$est_id with no strand data\n" if (! defined $this_est_strand);

    my $chrom_pos;   # this is the first position of the TSL/polyA-signal on the chromosome	
    if      ($FORWARD && $ORIENTATION) { # e.g. yk1304b05.5
      $diff_pos = $feature_site - $this_est_pos_start;
    } elsif ($FORWARD && !$ORIENTATION) { # e.g. EC017365 
      $diff_pos = $feature_site - $this_est_pos_start;
    } elsif ($REVERSE && $ORIENTATION) { # e.g. EC026249
      $diff_pos = $this_est_pos_start - $feature_site;
    } elsif ($REVERSE && !$ORIENTATION) { # e.g. yk827a08.3 also check CEESB85F
      $diff_pos = $this_est_pos_start - $feature_site;
    }
    $chrom_pos = $this_est_start + $diff_pos; # this is the first position of the TSL/polyA-site on the chromosome

    
#    $diff_pos = $feature_site - $this_est_pos_start;
#    my $chrom_pos = $this_est_start + $diff_pos; # this is the first position of the TSL/polyA-site on the chromosome
    print "feature_site $feature_site this_est_pos_start $this_est_pos_start\n" if ($verbose);
    print "diff_pos = $diff_pos this_est_start $this_est_start chrom_pos $chrom_pos\n" if ($verbose);
    if ($this_est_strand eq '-') {$chrom_pos--;} # GFF features always go start->end even on reverse sense

    # look to see if there is an existing Feature we can use - just loop through all of them
    my $feature_id;
    foreach my $pasft (@site_features) {
      $feature_id = $pasft->[0];
      my $pasft_start = $pasft->[1];	# chrom coords
      my $pasft_end = $pasft->[2];	# chrom coords always start < end
      my $pasft_strand = $pasft->[3];	# + or -
      if ($chrom_pos == $pasft_start && 
	  $this_est_strand eq $pasft_strand) {
	$got_a_feature_site_match = $feature_id;
	print "***** $est_id matches the feature $feature_id at $chrom_pos\n";
	last;
      }
    }

    my ($clone_name, $bases_upstream, $bases_downstream); 
    if ($got_a_feature_site_match eq "") { # not found an existing Feature to use
      print "----- $est_id has no matching feature at $chrom_pos\n";
      # Define a new feature ID to use
#      $feature_id = $next_feature_id;
#      $next_feature_id++;
      $feature_id = $feature_db->idCreate ;


}


=pod

=head2 NAME - find_anomalies.pl

=head1 USAGE

=over 4

=item find_anomalies.pl  [-options]

=back

This script populates the worm_anomaly mysql database with data describing some types of anomalies.

It can be run periodicaly e.g. every build especially useful if there
have been corrections made to the genomic sequence

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
