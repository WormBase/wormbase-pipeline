#!/software/bin/perl -w
#
# get_TSL_RNASeq_reads.pl
#
#
# Reads a fastq file looking for SL1/SL2 trans-splice sequences and
# writes out a file of reads with SL1 and SL2 subsequence removed ready
# for aligning to the genome.
#
#
#
#
# by Gary Williams



use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;


######################################
# variables and command-line options #
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $infile);
my ($outfile, $species);

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "test"        => \$test,
	    "verbose"     => \$verbose,
	    "store:s"     => \$store,
	    "infile=s"    => \$infile,
	    "outfile=s"   => \$outfile, 
	    "species=s"   => \$species,
	    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species,
			     );
}

my $log = Log_files->make_build_log( $wormbase );

# Display help if required
&usage("Help") if ($help);




################################################################


# the TSL names and sequences
my %SL = $wormbase->TSL;



# prepare the TSL sequences
my %minimal;
my %reverse_minimal;
my %reverse_SL;
my $MINIMUM_TSL = 8; # minimum length of TSL we will accept as evidence of a TSL
foreach my $sl_name (keys %SL) {
  my $seq = $SL{$sl_name};

  # get the reverse sense SL hash
  $reverse_SL{$sl_name} = $wormbase->DNA_string_reverse($seq);

  # get the minimal matching sequences
  $seq = substr($seq, -${MINIMUM_TSL});
  my $reverse_seq = $wormbase->DNA_string_reverse($seq);  
  push @{$minimal{$seq}}, $sl_name;
  push @{$reverse_minimal{$reverse_seq}}, $sl_name;

}


# read all of the hit read IDs
print "Reading hits\n";

open (OUT_SL, ">$outfile") || $log->log_and_die("can't open $outfile\n");

my $id;
my $seq;
my $line3;
my $line4;
my $tslname;
my $tsllen;

# stats
my $total_reads=0;
my $negative_reads=0;
my $sl1_reads=0;
my $sl2_reads=0;

if ($infile =~ /z$/) {
  open (READ, "gunzip -c $infile|") || $log->log_and_die("can't open read file: $infile\n");
} else {
  open (READ, "<$infile") || $log->log_and_die("can't open read file: $infile\n");
}

while (my $id_line = <READ>) {

  my ($id) = ($id_line =~ /(\S+)/);

#@SRR006514.1 length=36
#AAAGCTATGCGGATTATGTACTGAACTAGGATCTGG
#+SRR006514.1 length=36
#I?8:=9I).'&%&,/-+%+)#+#&&"%'#""""%%$

# read the sequence
    $seq = <READ>;
# and the other two lines
    $line3 = <READ>;
    $line4 = <READ>;

    $total_reads++;

    # check for sequence in forward sense
    ($tslname, $tsllen) = match_tsl($seq, '+');
    if (defined $tslname) {
      #print "$seq has $tslname, len = $tsllen\n\n";
      # add the TSL type to the ID
      chomp $id;
      $id =~ s/\//_/g; # change / to _ because tophat removed anything after the slash when it writes out the bam data
      $id  =~ s/$id/${id}\.\+\.${tslname}/;
      chomp $line3;
      $line3 =~ s/$id/${id}\.\+\.${tslname}/;	
      # remove the TSL from the sequence
      # and put back the 'AG' on the front of the sequence
      substr($seq, 0, $tsllen) = 'AG';
      # remove the TSL quality
      # stick the quality code 'II' on the front
      substr($line4, 0, $tsllen) = 'II';
      # print it out
      print OUT_SL "${id}\n${seq}${line3}\n${line4}";
      if ($tslname =~ /SL1/) {
	$sl1_reads++;
      } else {
	$sl2_reads++;
      }
# there are only 1% of negative reads with a TSL sequence and these
# are probably errors anyway because the negative reads are
# sequencing the 3' end of the fragments of mRNA and so are not likely
# to sequence the very 5' end.

    } else {
      # else check for sequence in reverse sense
      chomp $seq;
      ($tslname, $tsllen) = match_tsl($seq, '-');
      if (defined $tslname) {
	#print "$seq has $tslname, len = $tsllen\n\n";
	# add the TSL type to the ID
	chomp $id;
	$id =~ s/\//_/g; # change / to _ because tophat removed anything after the slash when it writes uot the bam data
	$id  =~ s/$id/${id}\.\-\.${tslname}/;
	chomp $line3;
	$line3 =~ s/$id/${id}\.\-\.${tslname}/;	
	# remove the TSL from the sequence
	# stick 'AG' on the end of the sequence
	substr($seq, - $tsllen) = 'CT';
	# remove the TSL quality
	# stick the quality code 'II' on the end
	chomp $line4;
	substr($line4, - $tsllen) = 'II';
	# print it out
	print OUT_SL "${id}\n${seq}\n${line3}\n${line4}\n";
	if ($tslname =~ /SL1/) {
	  $sl1_reads++;
	} else {
	  $sl2_reads++;
	}
	$negative_reads++;
      }
    }
  }

close(READ);
$log->write_to("\nFinished reading $infile\n");
close(OUT_SL);

$log->write_to("Total reads processed: $total_reads\n");
my $percent = $sl1_reads*100/$total_reads;
$log->write_to("SL1 reads: $sl1_reads, ($percent% of all reads)\n");
$percent = $sl2_reads*100/$total_reads;
$log->write_to("SL2 reads: $sl2_reads, ($percent% of all reads)\n");
$percent = $negative_reads*100/($sl1_reads+$sl2_reads+$negative_reads);
$log->write_to("Reverse sense reads: $negative_reads, ($percent% of TSL)\n");


$log->write_to("\nFinished all\n");

$log->mail();
exit(0);




###############################
# Prints help and disappears  #
###############################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}


######################################
# finds matches to the TSL sequences #
######################################

# Returns name of type of TSL and length of TSL at start of sequence

sub match_tsl {

  my ($seq, $sense)=@_;

  # searches for n-mers of the TLS 
  # min 8 so that we are fairly confident that we have a TSL
  if ($sense eq '+') {

    # look for a match to the minimal allowed 
    foreach my $minimum_seq (keys %minimal) {
      # see if it matches anywhere in the read sequence
      #    index($seq, $minimum_seq) != -1
      # is slightly faster than
      #    $seq =~ /${minimum_seq}/
      # (33 seconds vs 37 seconds to do 100 million test searches)
      if (index($seq, $minimum_seq) != -1) {    
	foreach my $slname (@{$minimal{$minimum_seq}}) {
	  my $len = length $SL{$slname};
	  for (my $i=$MINIMUM_TSL; $i <= $len; $i++) {
	    my $nmer=substr($SL{$slname}, -$i);
	    if ($seq =~ /^${nmer}/) {
	      return ($slname, $i);
	    }
	  }
	}
      }
    }
  } else {

    # look for a match to the minimal allowed 
    foreach my $minimum_seq (keys %reverse_minimal) {
      if (index($seq, $minimum_seq) != -1) {    
	foreach my $slname (@{$reverse_minimal{$minimum_seq}}) {
	  my $len = length $reverse_SL{$slname};
	  for (my $i=$MINIMUM_TSL; $i <= $len; $i++) {
	    my $nmer=substr($reverse_SL{$slname}, 0, $i);
	    if ($seq =~ /${nmer}$/) {
	      return ($slname, $i);
	    }
	  }
	}
      }
    }
  }
  return undef;
}



##########################################

__END__

=pod

=head2 NAME - make_unmatched_RNASeq.pl

=head1 USAGE

=over 4

=item make_unmatched_RNASeq.pl [-options]

Reads a BAM file to get the reads that matched, then goes through
the original read files writing out the files that didn't match and
which look like the have TSL or polyA sequences.

=back


find_intergenic.pl mandatory arguments:


=over 4

=item -output, Specifies the output file name

=back

find_intergenic.pl optional arguments:

=over 4

=back

=item -debug, Verbose/Debug mode

=item -test, Test mode

=item -help, Help pages

=back

=cut
