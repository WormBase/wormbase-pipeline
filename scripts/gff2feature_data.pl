#!/software/bin/perl -w
#
# Small script to convert GFF regions to Feature_data objects
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2012-06-22 08:56:52 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
#use Ace;
#use Sequence_extract;
use Coords_converter;
use Modules::Remap_Sequence_Change;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($input, $output, $species);
my %features;

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "input:s"    => \$input, # the GFF file
	    "output:s"   => \$output, 
	    "features:s"  => \%features, # feature(s) to select column 3 of the GFF, takes value of method  
	    "species:s"  => \$species,
	    );
# use multiple features definitions e.g. -features Poly-A=RNASeq_polyA -features SL1=RNASeq_SL1 -features SL2=RNASeq_SL2

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


#################################
# check input arguments
#################################


$species = $wormbase->full_name;


#################################

my $database = $wormbase->autoace;
my $coords = Coords_converter->invoke($database, 0, $wormbase);
my $virtual;

# suck the data in
open (IN, "<$input") || die "Can't open $input\n";
open (ACE, ">$output") || die "Can't open $output\n";

my (@tiles, @whole_chromosome);

my $sequence = '';
my $sequence_len = 0;
my %scores; # a list of the scores for each feature, so that we can find the median scores

while (my $line = <IN>) {
  
  if ($line =~ /^#/) {next;}

  my @cols = split /\s+/, $line;
  my $chrom = $cols[0];
  my $ft    = $cols[2];
  my $start = $cols[3];
  my $end   = $cols[4];
  my $reads = $cols[5];
  my $sense = $cols[6];
  
  # ignore lines that are not from the feature we want
  if (!exists $features{$ft}) {next;}
  my $method = $features{$ft};
      
  # assume we are dealing with strand-sensitive data, in which case we
  # need to be able to distinguish strands by making the smallest
  # region 2 bases long so that reverse strand is different to forward
  # strand.
  if ($start == $end) {
    if ($sense eq '+') {$end++}
    if ($sense eq '-') {$start--}
  }
  if ($sense eq '-') {
    ($end, $start) = ($start, $end);
  }
    
  if ($chrom ne $sequence) { # new sequence
 
    write_tiles(\@tiles, \@whole_chromosome, $sequence, $sequence_len); #  write the old data

    $sequence = $chrom;
    $sequence_len = initialise_tiles($chrom, \@tiles, $coords);

  }

  store_feature_in_tile(\@tiles, \@whole_chromosome, $method, $start, $end, $reads, $ft);


}

# write the last sequence
write_tiles(\@tiles, \@whole_chromosome, $sequence, $sequence_len); #  write the old data

# work out the median score and write the Method object
write_method();


close(ACE);
close(IN);



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
# find the tile to store the Feature in
sub store_feature_in_tile {
  my ($tiles_aref, $whole_chromosome_aref, $method, $start, $end, $reads, $text) = @_;
  # $text is the $feature from the GFF column 3

  # store the score so that we can get the median score
  if ($reads =~ /^\d+$/) {
    push @{$scores{$method}}, $reads;
  }

  my $found = 0;
  for( my $tile_idx = 1; $tile_idx <= @{$tiles_aref}; $tile_idx++) {
    my $tile = $tiles_aref->[$tile_idx-1];
    if ($start < $end) {
      if ($start > $tile->{start} && $end <= $tile->{end}) { # find the tile containing this forward Feature
        push @{$tile->{segs}}, [$method, $start - $tile->{start} + 1, $end - $tile->{start} + 1, $reads, $text];
        $found = 1;
      }
    } else {
      if ($end > $tile->{start} && $start <= $tile->{end}) { # find the tile containing this reverse Feature
        push @{$tile->{segs}}, [$method, $start - $tile->{start} + 1, $end - $tile->{start} + 1, $reads, $text];
        $found = 1;
      }
    }
  }
  if (!$found) { # it falls between two tiles, so place it on the top-level Sequence
    push @{$whole_chromosome_aref}, [$method, $start, $end, $reads, $text];
  }
}

##########################################
sub write_tiles {
  my ($tiles_aref, $whole_chromosome_aref, $sequence, $sequence_len) = @_;


  # output the new Sequence lines

  foreach my $method (values %features) {
    if (!defined $sequence || $sequence eq '') {last}
    my $virtual = "${sequence}:${method}";
    
    my @sequence_out;
    my @feature_out;

    if (scalar @{$tiles_aref}) {
      push @sequence_out, "\nSequence : \"${sequence}\"\n";

      for(my $tile_idx = 1; $tile_idx <= @{$tiles_aref}; $tile_idx++) {
	my $tile = $tiles_aref->[$tile_idx-1];
	
	my $vseq = "${virtual}:$tile_idx";
	
	if (@{$tile->{segs}}) {
	  push @sequence_out, "S_Child Feature_data ". $vseq ." ". $tile->{start} ." ". $tile->{end} ."\n";
	  
	  push @feature_out, "\nFeature_data : \"$vseq\"\n";
	  foreach my $seg (@{$tile->{segs}}) {
	    if ($seg->[0] eq $method) {
	      push @feature_out, "Feature @$seg\n";
	    }
	  }
	}
      }
    }
    
    
    if (scalar @{$whole_chromosome_aref}) {
      push @sequence_out, "\nSequence : \"${sequence}\"\n";
      push @sequence_out, "S_Child Feature_data ${virtual} 1 $sequence_len\n";
      push @feature_out, "\nFeature_data : ${virtual}\n";
      foreach my $seg (@{$whole_chromosome_aref}) {
	if ($seg->[0] eq $method) {
	  push @feature_out,  "Feature @$seg\n";
	}
      }
    }
    
    print ACE @sequence_out;
    print ACE "\n"; # acezip.pl concatenates another line to the last line if this is not blank
    
    print ACE @feature_out;
    print ACE "\n";

  }

  @{$tiles_aref} = ();
  @{$whole_chromosome_aref} = ();
}
##########################################
sub initialise_tiles {
  my ($sequence, $tiles_aref, $coords) = @_;

  my $chr_len = $coords->Superlink_length($sequence);
  if (!defined $chr_len) {$log->log_and_die("Can't find the length of the Sequence $sequence\n")}

  for(my $i=0; $i < $chr_len; $i += 300000) {
    my $chr_start = $i + 1;
    my $chr_end = $chr_start + 300000 - 1;
    $chr_end = $chr_len if $chr_end > $chr_len;
    push @{$tiles_aref}, {
                  start => $chr_start, 
                  end   => $chr_end,
                  segs  => [],
    }
  }
  return $chr_len;
}
##########################################
# work out the median score and write the Method object

sub write_method {

  foreach my $feature (keys %features) {
    my $method = $features{$feature};
    my $median_score = median(@{$scores{$method}});
    my $score_max = 10 * $median_score;
    print "$method\t median score: $median_score\n";
    my $colour = get_colour($method); # convert the method name into a colour value

    print ACE "\n";
    print ACE "Method : $method\n";
    print ACE "Remark \"This data was produced by $ENV{USER} with the script gff2feature_data.gff from the data file $input using the feature column '$feature'.\"\n";
    print ACE "Show_up_strand\n";
    print ACE "Score_by_width\n";
    print ACE "Score_bounds 1 $score_max\n";
    print ACE "Overlap\n";
    print ACE "Right_priority 1.5\n";
    print ACE "Colour $colour\n";
    print ACE "\n";

  }


}
##########################################
# return the median value of a list of values
sub median {

    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}
##########################################
# return colour from a string that is hased to a value from 1 to 16
sub get_colour {
  my ($text) = @_;

  my @colours = qw(
		    WHITE
		    BLACK
		    LIGHTGRAY
		    DARKGRAY
		    RED
		    GREEN
		    BLUE
		    YELLOW
		    CYNA
		    MAGENTA
		    LIGHTRED
		    LIGHTGREEN
		    LIGHTBLUE
		    DARKRED
		    DARKGREEN
		    DARKBLUE
		    PALERED
		    PALEGREEN
		    PALEBLUE
		    PALEYELLOW
		    PALECYAN
		    PALEMAGENTA
		    BROWN
		    ORANGE
		    PALEORANGE
		    PURPLE
		    VIOLET
		    PALEVIOLET
		    GRAY
		    PALEGRAY
		    CERISE
		    MIDBLUE
		);

  my $value = 0;
  foreach my $letter (split //, $text) {
    $value += ord $letter
  }
  $value %= 32;
  return $colours[$value];
}
##########################################

# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - gff3ace.pl

=head1 USAGE

=over 4

=item gff3ace.pl  [-options]

=back

This script reads in a GFF file of locations and writes out an Feature_data ACE file.




script_template.pl MANDATORY arguments:

=over 4

=item -input input file of gene predictions in GFF3 format

=back

=over 4

=item -output output ACE file of CDS objects

=back

=over 4

=item -feature feature GFF3 field name to fnd and use and the method to use. specify as many times as you wish e.g. -features Poly-A=RNASeq_polyA -features SL1=RNASeq_SL1 -features SL2=RNASeq_SL2

=back




script_template.pl  OPTIONAL arguments:

=over 4

=item -species species_name. By default, this script will write ACE data specifying the species as 'elegans'. This specifies a different species.

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

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Gary Williams

=back

=cut
