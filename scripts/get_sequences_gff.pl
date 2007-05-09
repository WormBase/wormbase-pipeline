#!/nfs/disk100/wormpub/bin/perl -w
#
# get_sequences_gff.pl
# 
# by Gary Williams                         
#
# This is a script to emulate the SAGE processing script
# get_sequences_gff.pl
#
# It gets details of the exons structure and sequence of transcripts
# by using the GFF files and autoace
#
# It is run to process each chromosome individually under LSF by the
# script map_tags.pl

#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2007-05-09 15:00:42 $      

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

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $chromosome);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"      => \$store,
	    "chromosome:s" => \$chromosome, # optional single chromosome to be processed
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

# some database paths
my $currentdb = $wormbase->database('current');
my $ace_dir = $wormbase->autoace;     # AUTOACE DATABASE DIR


# other paths
my $tace = $wormbase->tace;        # TACE PATH


##########################
# MAIN BODY OF SCRIPT
##########################

#my $database_path = $currentdb;	
my $database_path = $ace_dir;

print "connecting to server...";
my $program = $wormbase->tace;  # full path to tace; change as appropriate
my $db = Ace->connect(-path => $database_path,  -program => $program) || $log->log_and_die("Connection failure: " . Ace->error);  # local database
my $coords = Coords_converter->invoke($database_path, 0, $wormbase);
my $seq_obj = Sequence_extract->invoke($database_path, 0, $wormbase);
my %clonesize = $wormbase->FetchData('clonesize');
print "done\n";



my @chroms;
if ($chromosome) {
  @chroms = ($chromosome);
} else {
  @chroms = $wormbase->get_chromosome_names(-mito =>1, -prefix => 1);
}

foreach my $chrom ( @chroms ) {
  my $output = "get_sequences_gff_${chrom}.out";
  open (OUT, ">$output") or $log->log_and_die("Cant open $output: !$\n");
  my @genes = &read_chromosome($chrom);
  &process_chromosome($chrom, @genes);
  close(OUT) or $log->error("Problem closing $output: !$\n");
}


# Close log files and exit
$log->write_to("get_sequences.pl finished @chroms\n");

# get a clean exit status by closing and undef'ing things
$db->close;
$log->mail();
$wormbase = undef;
$log = undef;
exit(0);






##############################################################
#
# Subroutines
#
##############################################################



##########################################
# read the data from a chromosomes GFF file and process it
# read_chromosome($chrom)

sub read_chromosome {
  my ($chromosome) = @_;

  my @genes;

  my $source = "Coding_transcript";
  my $type = "protein_coding_primary_transcript";

  &read_chromosome_file($source, $type, \@genes);

  $source = "Pseudogene";
  $type = "Pseudogene";

  &read_chromosome_file($source, $type, \@genes);

  $source = "Non_coding_transcript";
  $type = "nc_primary_transcript";

  &read_chromosome_file($source, $type, \@genes);


  # sort the genes by start position
  return sort {$a->[2] <=> $b->[2]} @genes

}


##########################################
# read the data from a chromosomes GFF file
# read_chromosome_file($source, $type, $genes_aref)

sub read_chromosome_file {
  my ($insource, $intype, $genes_aref) = @_;

  print "reading chromosome $chromosome $insource\n";
  my $gff = "$database_path/GFF_SPLITS/${chromosome}_${insource}.gff";

  open (GFF, "< $gff") || die "Can't open $gff\n";
# loop through the file looking for lines of type Coding Transcript,
# Pseudogene or Non-coding Transcript
# store type, start, stop, strand, ID
  while (my $line = <GFF>) {
    chomp $line;
    if ($line =~ /^\s*$/) {next;}
    if ($line =~ /^#/) {next;}
    my @f = split /\t/, $line;
    my ($chromosome, $source, $type, $start, $end, $sense) = ($f[0], $f[1], $f[2], $f[3], $f[4], $f[6]);
    if ($source ne $insource) {next;}
    if ($type ne $intype) {next;}
    my ($id) = ($f[8] =~ /\S+\s+(\S+)/);
    $id =~ s/\"//g;       # remove quotes
    my $strand = ($sense eq '+')?1:-1;
    print "$chromosome, $source, $type, $start, $end, $sense, $id\n" if ($verbose);
	
    push @{$genes_aref}, [($source, $type, $start, $end, $strand, $id)]
	
  }

  close(GFF);



}


##########################################
# foreach gene
# get the start, end position
# get the CDS exons
# extend the region to make a UTR if the CDS is the same size
# but don't extend it over another gene
# get the sequence of this region
# get the enclosing clone/superlink/chromosome
#  &process_chromosome($chrom, @genes);

sub process_chromosome {
  my ($chromosome, @genes) = @_;
  my $class;

  print "processing chromosome $chromosome\n";

  for (my $i = 0; $i < @genes; $i++) {
#    print ".";

    my ($source, $type, $start, $stop, $strand, $id) = @{$genes[$i]};

    my $seq = '';

    if ($source eq 'Coding_transcript') {
      $class="Transcript";

      my $obj = $db->fetch($class => $id);
      my $cds_ace = $obj->Corresponding_CDS;
      my ($cds_clone, $cds_start, $cds_stop) = &get_details($cds_ace);

      # set up the UTR flags
      my $real5;
      my $real3;

      # convert $cds_start, $cds_stop to chromosomal coords
      $cds_start = $coords->Coords_2chrom_coords($cds_clone, $cds_start);
      $cds_stop  = $coords->Coords_2chrom_coords($cds_clone, $cds_stop);

      if ($start == $cds_start) {            # if transcript start is the same as cds start, no 5' UTR exists - will be estimated
	$real5='';
      } else {
	$real5='yes';
      }
      if ($stop == $cds_stop) {              # if transcript stop is the same as cds stop, no 3' UTR exists - will be estimated
	$real3='';
      } else {
	$real3='yes';
      }

      # get details of the current transcript to get exons
      my @exons_start = $obj->Source_exons;
      my @exons_stop = $obj->Source_exons(2);
      # convert to chromosome coords
      for (my $i = 0; $i < @exons_start; $i++) {
	if ($strand == 1) {
	  $exons_start[$i] = $exons_start[$i] + $start - 1;
	  $exons_stop[$i] = $exons_stop[$i] + $start - 1;
	} else {
	  # want the 5' exon in $exons_start[0]
	  $exons_start[$i] = $stop - $exons_start[$i] + 1;
	  $exons_stop[$i] = $stop - $exons_stop[$i] + 1;
	}
      }

      $seq = $obj->asDNA();     # NB already reverse complemented if in reverse sense
      $seq =~ s/^\n>\S+//;	# remove the title line
      $seq =~ s/\n//g;		# remove newlines
      $seq = lc $seq;

      # get the average size of UTR regions if the transcipt has the
      # same end positions as the CDS - adjust the @exons_start and
      # @exons_stop limits and the $clone_start and $clone_end
      # positions and add the UR to $seq

      # Based on latest data, 5' and 3' end UTR length can be defined
      # as 749 and 700 as 95 percentile; these are used to produce
      # estimated UTRs
      my ($leftUTRAveSize, $rightUTRAveSize) = (749, 700); 

      my ($new_start, $new_stop, $seq5, $seq3, $absstart, $absstop);
      if (!$real5) {
	if ($strand == 1) {        # 5' side first
	  $new_start=$start-1-$leftUTRAveSize; # define a 5'UTR
	  $new_stop=$start-1;
	} else {		# reverse strand
	  $new_start=$stop+1; # define a 5'UTR
	  $new_stop=$stop+1+$leftUTRAveSize;
	}

	# get the position of the closest gene to $obj on the same strand on the clone/superlink/chromosome
	my $closest = &get_closest_gene('up', $strand, $new_start, $new_stop, $i, @genes);

	if ($closest) {
	  if ($strand == 1) {
	    if ($closest < $start) {   # when gene models overlap, just get the length...
	      $new_start=$closest+1;
	      $new_stop=$start-1;
	    }
	  }
	  else {
	    if ($closest > $stop) {
	      $new_start=$stop+1;
	      $new_stop=$closest-1;
	    }
	  }

# get the sequence of the UTR
	  my $subseq_length = $new_stop - $new_start+1;
	  if ($strand == 1) {
	    $seq5= uc $seq_obj->Sub_sequence($chromosome, $new_start-1, $subseq_length); # we want $new_start-1 
	  }
	  else {
	    $seq5= reverse uc $seq_obj->Sub_sequence($chromosome, $new_start-1, $subseq_length); # we want $new_start-1
	    $seq5=~tr/tacgTACG/atgcATGC/;
	  }
	  if ($seq5) {
	    $seq=$seq5.$seq;
	    if ($strand == 1) {
	      $start = $start - length $seq5;
	      $exons_start[0] -= length $seq5;
	    }
	    else {
	      $start = $start + length $seq5;
	      $exons_start[0] += length $seq5;
	    }
	  }	  
	}
      }

      if (!$real3) {
	if ($strand == 1) {	# 3' side now
	  $new_start=$stop+1;
	  $new_stop=$stop+1+$rightUTRAveSize;
	} else {		# reverse strand
	  $new_start=$start-1-$rightUTRAveSize;
	  $new_stop=$start-1;
	}

	# get the position of the closest gene to $obj on the same strand on the clone/superlink/chromosome
	my $closest = &get_closest_gene('down', $strand, $new_start, $new_stop, $i, @genes);

	if ($closest) {
	  if ($strand == 1) {
	    if ($closest > $stop) {   # when gene models overlap, just get the length...
	      $new_start=$stop+1;
	      $new_stop=$closest-1;
	    }
	  }
	  else {
	    if ($closest < $start) {
	      $new_start=$closest+1;
	      $new_stop=$start-1;
	    }
	  }

# get the sequence of the UTR
	  my $subseq_length = $new_stop - $new_start+1;
	  if ($strand == 1) {
	    $seq3= uc $seq_obj->Sub_sequence($chromosome, $new_start-1, $subseq_length); # we want $new_start-1 
	  }
	  else {
	    $seq3= reverse uc $seq_obj->Sub_sequence($chromosome, $new_start-1, $subseq_length); # we want $new_start-1
	    $seq3=~tr/tacgTACG/atgcATGC/;
	  }
			
	  # truncate 35 nt after polyA signal in the 3' UTR
	  my $tail=$seq3=~/[a-z].+aataaa\S{35}(\S+)$/i ? $1 : '';
	  $seq3=~s/$tail//;
	  if ($tail) {
	  }

	  if ($seq3) {
	    $seq=$seq.$seq3;
	    if ($strand == 1) {
	      $stop = $stop + length $seq3;
	      $exons_stop[$#exons_stop] += length $seq3;
	    }
	    else {
	      $stop = $stop - length $seq3;
	      $exons_stop[$#exons_stop] -= length $seq3;
	    }
	  }	  
	}
      }

      # get enclosing clone/superlink/chromosome details
      my ($clone, $clone_start, $clone_end) = $coords->LocateSpan($chromosome, $start, $stop);
      my $clonelen = &get_clone_len($clone);

      # want to convert from chromosome coords to this clone
      my ($chrom_dummy, $clone_offset) = $coords->CloneOffset($clone);

      print OUT ">$id\tCoding transcript\t$clone\t$clonelen";
      if ($real3 and $real5) {
	print OUT "\treal5UTR-real3UTR";
      }
      elsif ($real3) {
	print OUT "\testimated5UTR-real3UTR";
      }
      elsif ($real5) {
	print OUT "\treal5UTR-estimated3UTR";
      }
      else {
	print OUT "\testimated5UTR-estimated3UTR";
      }

      # the exons are in chromosomal coords - want to knock off start of
      # clone on chromosome to get the position relative to the clone
      print OUT "\t$strand\t$clone_start\t$clone_end";
      for (my $i=0; $i < @exons_start; $i++) {
	if (! defined $exons_start[$i]) {print "$id Not defined exons_start[$i]\n";}
	if (! defined $exons_stop[$i]) {print "$id Not defined exons_stop[$i]\n";}
	if (! defined $clone_offset) {print "$id Not defined clone_offset\n";}
	if ($strand == 1) {
	  print OUT "\t", $exons_start[$i]-$clone_offset+1, "\t", $exons_stop[$i]-$clone_offset+1;
	} else {
	  print OUT "\t", $exons_start[$i]-$clone_offset+1, "\t", $exons_stop[$i]-$clone_offset+1;
	}
      }
      
      print OUT "\n$seq\n";


    } elsif ($source eq 'Pseudogene') {
      $class="Pseudogene";

      # get enclosing clone/superlink/chromosome details
      my ($clone, $clone_start, $clone_end) = $coords->LocateSpan($chromosome, $start, $stop);
      my $clonelen = &get_clone_len($clone);

      print OUT ">$id\tPseudogene\t$clone\t$clonelen\tN/A\t$strand\t$clone_start\t$clone_end";

      # get details of the current transcript to get exons
      my $obj = $db->fetch($class => $id);
      my @exons_start = $obj->Source_exons;
      my @exons_stop = $obj->Source_exons(2);

      # convert to chromosome coords
      for (my $i = 0; $i < @exons_start; $i++) {
	if ($strand == 1) {
	  $exons_start[$i] = $exons_start[$i] + $start - 1;
	  $exons_stop[$i] = $exons_stop[$i] + $start - 1;
	} else {
	  # want the 5' exon in $exons_start[0]
	  $exons_start[$i] = $stop - $exons_start[$i] + 1;
	  $exons_stop[$i] = $stop - $exons_stop[$i] + 1;
	}
      }

      # want to convert from chromosome coords to this clone
      my $clone_offset = $coords->CloneOffset($clone);

      for (my $i=0; $i < @exons_start; $i++) {
	if ($strand == 1) {
	  print OUT "\t", $exons_start[$i]-$clone_offset+1, "\t", $exons_stop[$i]-$clone_offset+1;
	} else {
	  print OUT "\t", $exons_start[$i]-$clone_offset+1, "\t", $exons_stop[$i]-$clone_offset+1;
	}
      }
      
      $seq = $obj->asDNA();     # NB already reverse complemented if in reverse sense
      $seq =~ s/^\n>\S+//;	# remove the title line
      $seq =~ s/\n//g;		# remove newlines
      $seq = lc $seq;
      
      print OUT "\n$seq\n";


    } else {			# 'Non_coding_transcript'
      $class="Transcript";

      # get enclosing clone/superlink/chromosome details
      my ($clone, $clone_start, $clone_end) = $coords->LocateSpan($chromosome, $start, $stop);
      my $clonelen = &get_clone_len($clone);

      print OUT ">$id\tNon-coding transcript\t$clone\t$clonelen\tN/A\t$strand\t$clone_start\t$clone_end";

      # get details of the current transcript to get exons
      my $obj = $db->fetch($class => $id);
      my @exons_start = $obj->Source_exons;
      my @exons_stop = $obj->Source_exons(2);


      # convert to chromosome coords
      for (my $i = 0; $i < @exons_start; $i++) {
	if ($strand == 1) {
	  $exons_start[$i] = $exons_start[$i] + $start - 1;
	  $exons_stop[$i] = $exons_stop[$i] + $start - 1;
	} else {
	  # want the 5' exon in $exons_start[0]
	  $exons_start[$i] = $stop - $exons_start[$i] + 1;
	  $exons_stop[$i] = $stop - $exons_stop[$i] + 1;
	}
      }

      # want to convert from chromosome coords to this clone
      my $clone_offset = $coords->CloneOffset($clone);

      for (my $i=0; $i < @exons_start; $i++) {
	if ($strand == 1) {
	  print OUT "\t", $exons_start[$i]-$clone_offset+1, "\t", $exons_stop[$i]-$clone_offset+1;
	} else {
	  print OUT "\t", $exons_start[$i]-$clone_offset+1, "\t", $exons_stop[$i]-$clone_offset+1;
	}
      }
      
      $seq = $obj->asDNA();     # NB already reverse complemented if in reverse sense
      $seq =~ s/^\n>\S+//;	# remove the title line
      $seq =~ s/\n//g;		# remove newlines
      $seq = lc $seq;
      
      print OUT "\n$seq\n";

    }
  }
}

##########################################
# get the length of a clone/superlink/chromosome
# $len = get_clone_len($clone)

sub get_clone_len {
  my ($clone) = @_;

  my $clonelen = $clonesize{$clone};

  if (! defined $clonelen) {
    if ($clone =~ /SUPERLINK/ || $clone =~ /CHROMOSOME/) {
      # get the Superlink lengths from the Coords_converter data
      $clonelen = $coords->Superlink_length($clone);
    } else {
      die "undef returned for length of $clone\n";
    }
  }

  return $clonelen;

}

##########################################
# get details of a CDS
# my ($cds_clone, $cds_start, $cds_stop) = &get_details($cds_ace);

sub get_details {
  my ($obj) = @_;

  my $clone = $obj->at('SMap.S_parent.Sequence[1]'); # get the clone/superlink/chromosome this object is on

  my $quoted_obj = $obj->name();
  $quoted_obj =~ s/\./\\./g;		# replace . with \.

# for some reason we need to fetch the clone as a Sequence object even
# though we could use $clone to get the length of the Sequence $clone
# just now.
  my $fetched_clone = $db->fetch(Sequence => $clone->name());
  my ($detail, $start, $stop) = $fetched_clone->at("SMap.S_child.CDS_child.$quoted_obj")->row; # start and end position in the clone/superlink/chromosome

  print "start, stop of $obj = $start, $stop\n" if ($verbose);

  return ($fetched_clone, $start, $stop);

}

##########################################

# get the position of the closest gene to $obj on the same strand on
# the clone/superlink/chromosome
# my $closest = &get_closest_gene('up', $strand, $new_start, $new_stop, $i, @genes);

sub get_closest_gene {

  my ($direction, $strand, $new_start, $new_stop, $i, @genes) = @_;

  my $closest = 0;

  print "getting closest gene to $strand, $new_start, $new_stop\n" if ($verbose);

  if (($direction eq 'up' && $strand == 1) || ($direction eq 'down' && $strand == -1)) { # before the gene on chromosome
    while ($i > 0) {
      my ($source, $type, $prev_start, $prev_end, $prev_strand, $id) = @{$genes[--$i]}; # get previous gene
      # is it on the same strand?
      if ($prev_strand != $strand) {next;}
      # does the previous gene end in the start-stop region?
      if ($prev_end >= $new_start && $prev_end <= $new_stop) {
	$closest = $prev_end;
      }
      # if it ended before the region then stop looking
      if ($prev_end < $new_start) {
	last;
      }

    }
  } else {			# after the gene on the chromosome
    while ($i < $#genes) {
      my ($source, $type, $next_start, $next_end, $next_strand, $id) = @{$genes[++$i]}; # get next gene
      # is it on the same strand?
      if ($next_strand != $strand) {next;}
      # does the next gene start in the start-stop region?
      if ($next_start >= $new_start && $next_start <= $new_stop) {
	$closest = $next_start;
      }
      # if it started after the region then stop looking
      if ($next_start > $new_stop) {
	last;
      }

    }
  }

  return $closest;

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

=head2 NAME - get_sequences_gff.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

=back

This script writes out a file giving details of the coding
transcripts, non-coding transcripts and pseudogenes including the
sequences, exon positions , clone or superlink they are in and the
expected UTR region

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

=over 4
    
=item -chromosome, single chromosome to process.

=back


=head1 REQUIREMENTS

=over 4

=item transcript details in autoace and GFF files, 

=back

=head1 AUTHOR

=over 4

=item Gary Williams

=back

=cut
