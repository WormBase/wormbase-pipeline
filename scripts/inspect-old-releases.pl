#!/nfs/disk100/wormpub/bin/perl
#
# inspect-old_releases.pl                           
# 
# by Gary Williams                         
#
# This inspects the difference in the clone positions and genomic
# sequence of two relreases
#
# Usage:
# foreach r ( 147 146 145 144 143 142 141 140 139 138 137 136 135 134 133 132 131 130 129 128 127 126 125 124 123 122 121 120 119 118 117 116 115 )
# foreach? @ q = $r
# foreach? @ q--
# foreach? echo "$q $r"
# foreach? ./inspect-old-releases.pl -debug gw3 -version $r -database1 ~wormpub/gary/Archeology/WS{$q} -database2 ~wormpub/gary/Archeology/WS{$r}
# foreach? end
#
# Last updated by: $Author: pad $     
# Last updated on: $Date: 2006-05-02 14:29:01 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
#use Ace;
#use Sequence_extract;
#use Coords_converter;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($database1, $database2, $version);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store"      => \$store,
	    "database1=s" => \$database1,
	    "database2=s" => \$database2, 
	    "version=i"  => \$version,
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


##########################
# MAIN BODY OF SCRIPT
##########################

my $outfile = "/nfs/disk100/wormpub/CHROMOSOME_DIFFERENCES/sequence_differences.WS$version";

open (OUT, "> $outfile") || die "Can't open $outfile";

my @chromosomes = qw(I II III IV V X);
#my @chromosomes = qw(X);
foreach my $chromosome (@chromosomes) {
  my @differences = ();

  print "Chromosome: $chromosome\n";
  print OUT "Chromosome: $chromosome\n";

  my @chromosome_pair = &get_chromosomes($chromosome, $database1, $database2);
  @differences = &compare_chromosomes(@chromosome_pair);

  foreach my $diffs (@differences) {
    my ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped);
    # the mismatch_start value is the start of the mismatch, it is the first position which doesn't match
    # the mismatch_end value is the base past the end of the mismatch region, the first base which matches again
    ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped) = @$diffs;
    print "old= $mismatch_start1, $mismatch_end1 len=$len1\tnew= $mismatch_start2, $mismatch_end2 len=$len2 flipped=$flipped\n";
    print OUT "$mismatch_start1\t$mismatch_end1\t$len1\t$mismatch_start2\t$mismatch_end2\t$len2\t$flipped\n";
  }
}

close (OUT);

# Close log files and exit
$log->write_to("Finished.\n");

$log->mail();
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

# compare_chromosomes(@chromosome_pair)
# compares two chromosomes for differences
# returns list of differences

sub compare_chromosomes {
  my @chrom_pair = @_;

  my $start1 = 0;
  my $start2 = 0;

  my @differences = ();

  while(1) {
    ($start1, $start2) = snake($start1, $start2, @chrom_pair);
    if (!defined $start1) {
#      print "No difference found\n";
      return @differences;
    } else {
#      print "Difference found at $start1, $start2\n";
      # the start positions are the first mismatch
      # the end positions are the first match after the mismatch
      my ($mismatch_start1, $mismatch_end1, $mismatch_start2, $mismatch_end2)
	  = seeker($start1, $start2, @chrom_pair);
      $start1 = $mismatch_end1;
      $start2 = $mismatch_end2;

      # get the lengths of the mismatch regions
      my $len1 = $mismatch_end1 - $mismatch_start1;
      my $len2 = $mismatch_end2 - $mismatch_start2;

      # check for flipped clones
      my ($flipped, $flipped_start1, $flipped_start2, $flipped_len) 
	  = &check_flipped($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, @chrom_pair);

      # deal with problems of flipped regions
      if (!$flipped || ($flipped && !defined $flipped_start1)) {
	# output the simple cases of either no flipped region, or a flipped region with no indels
	push @differences, [$mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, $flipped];
      } else {
	# now have a flipped region with indels to sort out
	# is there an indel at the start?
	if ($flipped_start1 != $mismatch_start1 || $flipped_start2 != $mismatch_start2) {
	  #print "indel at start of flipped region\n";
	  $len1 = $flipped_start1 - $mismatch_start1;
	  $len2 = $flipped_start2 - $mismatch_start2;
	  push @differences, [$mismatch_start1, $flipped_start1, $len1, $mismatch_start2, $flipped_start2, $len2, 0];

	  # reset the mismatch_start sites for the flipped output
	  $mismatch_start1 = $flipped_start1;
	  $mismatch_start2 = $flipped_start2
	}

	# output the flipped region
	my $flipped_end1 = $mismatch_start1 + $flipped_len;
	my $flipped_end2 = $mismatch_start2 + $flipped_len;
	push @differences, [$mismatch_start1, $flipped_end1, $flipped_len, $mismatch_start2, $flipped_end2, $flipped_len, 1];

	# is there an indel at the end?
	if ($flipped_end1 != $mismatch_end1 || $flipped_end2 != $mismatch_end2) {
	  #print "indel at end of flipped region\n";
	  $len1 = $mismatch_end1 - $flipped_end1;
	  $len2 = $mismatch_end2 - $flipped_end2;
	  push @differences, [$flipped_end1, $mismatch_end1, $len1, $flipped_end2, $mismatch_end2, $len2, 0];
	}
      }

      #print "old= $mismatch_start1, $mismatch_end1 len=$len1\tnew= $mismatch_start2, $mismatch_end2 len=$len2 flipped=$flipped\n";
      #print "On track again at $start1, $start2\n";
    }
  }

}
##########################################

sub DNA_string_reverse {
  my $revseq = reverse shift;
  $revseq =~ tr/a/x/;
  $revseq =~ tr/t/a/;
  $revseq =~ tr/x/t/;
  $revseq =~ tr/g/x/;
  $revseq =~ tr/c/g/;
  $revseq =~ tr/x/c/;
  return ($revseq);
}

##########################################
# check_flipped($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len, @chrom_pair)
# looks for inverted regions
# return array:
#               1 for flipped, 0 if not
#               start of first sequence's flipped region
#               start of second sequence's flipped region
#               length of flipped region

sub check_flipped {

  my ($mismatch_start1, $mismatch_end1, $len1, $mismatch_start2, $mismatch_end2, $len2, @chrom_pair) = @_;
  
  my $flipped = 0;

  if ($len1 > 1000 && $len2 > 1000) {
    if ($len1 == $len2) {

      if (substr($chrom_pair[0], $mismatch_start1,  $len1) eq
	  DNA_string_reverse(substr($chrom_pair[1], $mismatch_start2,  $len2))) {
	$flipped = 1;
      }
    } else {
      # check to see if the changed region has a reversed sequence together with an indel
      my $seq1 = substr($chrom_pair[0], $mismatch_start1,  $len1);
      my $seq2 = substr($chrom_pair[1], $mismatch_start2,  $len2);
      my $revpos;
      if ($len1 > $len2) {
	$revpos = index($seq1, DNA_string_reverse($seq2));
	if ($revpos != -1) {
	  $flipped = 1;
	  return ($flipped, $mismatch_start1+$revpos, $mismatch_start2, $len2);
	}
      } else {
	$revpos = index($seq2, DNA_string_reverse($seq1));
	if ($revpos != -1) {
	  $flipped = 1;
	  return ($flipped, $mismatch_start1, $mismatch_start2+$revpos, $len1);
	}
      }
    }
  }

  return ($flipped, undef, undef, undef);

}

##########################################
# seeker(pos1, pos2, @chrom_pair)
# looks for local matches to get the alignment back as soon as possible
# returns the position of the match or undef if off the end of the sequences

sub seeker {
  my ($pos1, $pos2, @chrom_pair) = @_;

  my $match_pos1;
  my $match_pos2;
  my $chunk_size = 100000;
  my $subseq_size = 200;

  my $mismatch_start1;
  my $mismatch_start2;
  my $mismatch_end1;
  my $mismatch_end2;

  my $length1 = length $chrom_pair[0];
  my $length2 = length $chrom_pair[1];

  # get a chunk extending to twice the length of the chunk that the
  # mismatch was found in
  my $chunk1 = substr($chrom_pair[0], $pos1, $chunk_size*2);
  my $chunk2 = substr($chrom_pair[1], $pos2, $chunk_size*2);
  
  # first find the start of the mismatch region
  foreach (my $i = 0; $i < $chunk_size; $i++) {
    my $chr1 = substr($chrom_pair[0], $pos1+$i, 1);
    my $chr2 = substr($chrom_pair[1], $pos2+$i, 1);
    if ($chr1 ne $chr2) {
      $mismatch_start1 = $pos1+$i;
      $mismatch_start2 = $pos2+$i;
      last;
    }
  }

 # print "Mismatch start: $mismatch_start1, $mismatch_start2\n";

  # then find the end of the mismatch region, starting from the first mismatch base
#  print "find end of mismatch\n";
  foreach  (my $i = 0; $i < $chunk_size*2; $i++) {

    # if we are at the end of the chromosome, then simply leave now
    if ($mismatch_start2+$i+$subseq_size > $length1 || $mismatch_start1+$i+$subseq_size > $length2 ) {
      $mismatch_end1 = $length1;
      $mismatch_end2 = $length2;
      #print "Hit end of chromosome\n";
      return ($mismatch_start1, $mismatch_end1, $mismatch_start2, $mismatch_end2);
    }

    my $subseq1 = substr($chrom_pair[1], $mismatch_start2+$i, $subseq_size);
    my $subseq2 = substr($chrom_pair[0], $mismatch_start1+$i, $subseq_size);
    $match_pos1 = index($chunk1, $subseq1);
    $match_pos2 = index($chunk2, $subseq2);
    #print "Trying end: ", $pos1+$match_pos, " ", $mismatch_start2+$i, "\t";

    if ($match_pos1 != -1 && $pos1+$match_pos1 >= $mismatch_start1 && $mismatch_start2+$i >= $mismatch_start2) {	# exit loop when find match
      $mismatch_end1 = $pos1+$match_pos1;
      $mismatch_end2 = $mismatch_start2+$i;
      last;
    } elsif ($i == ($chunk_size*2) - 1) {
      #print "Off end of chunk while searching for mismatch end\n";
    }

    if ($match_pos2 != -1 && $mismatch_start1+$i >= $mismatch_start1 && $pos2+$match_pos2 >= $mismatch_start2) {	# exit loop when find match
      $mismatch_end1 = $mismatch_start1+$i;
      $mismatch_end2 = $pos2+$match_pos2;
      last;
    } elsif ($i == ($chunk_size*2) - 1) {
      #print "Off end of chunk while searching for mismatch end\n";
    }
  }

#  print "Mismatch end: $mismatch_end1, $mismatch_end2\n";

  return ($mismatch_start1, $mismatch_end1, $mismatch_start2, $mismatch_end2);
}


##########################################
# snake(pos1, pos2, @chrom_pair)
# runs down the chromosome, from the position pos1, pos2 until it hits a mismatch
# returns the position of the mismatch or undef if off the end of the sequences

sub snake {
  my ($pos1, $pos2, @chrom_pair) = @_;

  my $chunk_size = 10000;

  my $length1 = length $chrom_pair[0];
  my $length2 = length $chrom_pair[1];

  my $chunk1 = substr($chrom_pair[0], $pos1, $chunk_size);
  my $chunk2 = substr($chrom_pair[1], $pos2, $chunk_size);

  while ($chunk1 eq $chunk2) {
    if ( $pos1+$chunk_size > $length1 ||
	 $pos2+$chunk_size > $length2) {
      return (undef, undef);
    }
    $pos1 += $chunk_size;
    $pos2 += $chunk_size;
    $chunk1 = substr($chrom_pair[0], $pos1, $chunk_size);
    $chunk2 = substr($chrom_pair[1], $pos2, $chunk_size);
  }
  
  return ($pos1, $pos2);

}

##########################################

# get_chromosomes
# reads in the chromosome files from the two databases
# returns the pair of chromosome seqeunces

sub get_chromosomes {
  my ($chrom, $db1, $db2) = @_;
  my @chrom_pair;

  $chrom_pair[0] = read_chromosome($chrom, $db1);
  $chrom_pair[1] = read_chromosome($chrom, $db2);

  return @chrom_pair;
}

##########################################

# read_chromosome
# reads in a chromosome sequence from its file
# returns the sequence string

sub read_chromosome() {
  my ($chromosome, $db) = @_;
  my $seq;
                                                                                                                                                            
  my $file = "$db/CHROMOSOMES/CHROMOSOME_$chromosome.dna";

  my $old_rs = $/;              # save the current value of the record separator

  if (! open (SEQ, "<$file")) {	# try to open the sequence file
    open (SEQ, "/bin/gunzip -c $file.gz |") || die "Can't open file $file\n"; # ... or try to open the gzipped sequence file
  }
  <SEQ>;                        # skip the title line
  undef $/;                     # don't use record separator when reading file
  $seq = <SEQ>;                 # slurp up the whole file
  $seq =~ s/\n//g;              # remove newline characters
  close (SEQ);
  $/ = $old_rs;                 # restore the old value of the record separator
  return $seq;
}


##########################################

# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - inspect_old_releases.pl

=head1 USAGE

=over 4

=item inspect_old_releases.pl  [-options]

=back

This script reads in the chromosomes of two different acedb releases
and creates a file of differences between them, read for mapping GFF
files between releases using the script 'map_gff_between_releases.pl'.

script_template.pl MANDATORY arguments:

=over 4

=item -database1 The first (earlier) database to use

=back

=item -database2 The second (later) database to use

=back

=item -version The version number of the second database, e.g. 150. The first database MUST be the previous release e.g. 149

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

=item Gary Williams

=back

=cut
