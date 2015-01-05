#!/software/bin/perl -w
#
# remap_clone_feature_data.pl
# 
# by Gary Williams                         
#
# This takes a BUILD_DATA/MISC_DYNAMIC/*.ace file mapped to clone
# Feature_data and converts any coordinates that have changed between
# releases
#
# Last updated by: $Author: pad $     
# Last updated on: $Date: 2015-01-05 15:56:23 $      

use strict;                                     
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

use Modules::Remap_Sequence_Change;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($input, $output);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
            "input:s"    => \$input,
	    "output:s"   => \$output,
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

# Set up top level base directories (these are different if in test mode)
my $ace_dir         = $wormbase->autoace;     # AUTOACE DATABASE DIR

# some database paths
my $currentdb = $wormbase->database('current');

##########################
# read in the mapping data
##########################
my $version = $wormbase->get_wormbase_version;
print "Getting mapping data for WS$version\n";
my $assembly_mapper = Remap_Sequence_Change->new($version - 1, $version, $wormbase->species, $wormbase->genome_diffs);

##########################
# MAIN BODY OF SCRIPT
##########################


#
# Read the locations from the previous ace file
# in order to be able to remap them
#

#Sequence : "Y44F5A"
#S_Child Feature_data Y44F5A:TranscriptionallyActiveRegion 1 6146
#
#Sequence : "ZC178"
#S_Child Feature_data ZC178:TranscriptionallyActiveRegion 1 11560
#
#Feature_data : "H21P03:TranscriptionallyActiveRegion"
#Feature TAR_embryo_coelomocytes 14430 16077 1.99999600001E-6 modENCODE_458
#Feature TAR_embryo_AVA_neurons 4996 5204 0.0603798792402 modENCODE_459
#Feature TAR_embryo_GABA_motor_neurons 4996 5178 0.0702978594043 modENCODE_468
#Feature TAR_L3_L4_rectal_epithelial_cells 4996 5230 0.0664938670123 modENCODE_3440
#Feature TAR_L3_L4_PVD___OLL_neurons 4996 5255 0.107671784656 modENCODE_460
#Feature TAR_L3_L4_PVD___OLL_neurons 4460 4688 0.00553198893602 modENCODE_460
#Feature TAR_embryo_panneural 11415 11562 0.0791678416643 modENCODE_455

# or

#Sequence : "AC3"
#S_Child Feature_data AC3:RNASeq 1 19686
#S_Child Feature_data AC3:RNASeq_forward_reads 1 19686
#S_Child Feature_data AC3:RNASeq_reverse_reads 1 19686
#
#Feature_data : "AC3:RNASeq"
#Feature RNASeq 3507 3538 1
#Feature RNASeq 3539 3553 2
#Feature RNASeq 3554 3808 4
#Feature RNASeq 4165 4224 5
#Feature RNASeq 4276 4440 6
#Feature RNASeq 4485 4613 6
#Feature RNASeq 4688 4795 7
#Feature RNASeq 4847 4993 6
#Feature RNASeq 5454 5820 7
#Feature RNASeq 5821 5879 4
#Feature RNASeq 5880 5881 2

# or

#Feature_data : "Y44F5A:Polysome"
#Feature Polysome 202 237 1
#Feature Polysome 254 286 1
#Feature Polysome 287 289 3
#Feature Polysome 290 322 2
#Feature Polysome 336 338 2
#Feature Polysome 339 344 3
#
#Sequence : "Y44F5A"
#S_Child Feature_data Y44F5A:Polysome 1 6146
#
#Sequence : "ZC178"
#S_Child Feature_data ZC178:Polysome 1 11560



# start the coords converters
my $current_converter = Coords_converter->invoke($currentdb, 0, $wormbase);
my $autoace_converter = Coords_converter->invoke($ace_dir, 0, $wormbase);

# get the details
my $clone;
my %clonesize       = $wormbase->FetchData('clonesize'); 
my $clone_length;
my $feature_data_type; # e.g. <Polysome> in Feature_data : "Y44F5A:Polysome"
my @line;

open (IN, "< $input") || die "can't open input file $input\n";
open (OUT, "> $output") || die "can't open output file $output\n";
while (my $line = <IN>) {
  chomp $line;
  
  # Sequence : "Y44F5A"
  if ($line =~ /^Sequence\s+\:\s+(\S+)/) {
    $clone = $1; # use this $clone in subsequent "S_Child Feature_data" lines
    $clone =~ s/\"//g; # remove quotes
    print OUT "$line\n";

  # S_child
  } elsif ($line =~ /^S_Child Feature_data/) {
    @line = split /\s+/, $line;
    # get the superlink or clone length
    if ($clone =~ /CHROMOSOME/) {
      $clone_length = &get_chrom_length($clone);
    } else {
      $clone_length = $clonesize{$clone};
    }    
    $line[4] = $clone_length;
    $line = join " ", @line;
    print OUT "$line\n";
    

  # Feature_data : "Y44F5A:Polysome"
  } elsif ($line =~ /^Feature_data\s+\:\s+(\S+)/) {
    ($clone, $feature_data_type) = (split /\:/,$1); # use this $clone in subsequent "Feature" lines
    $clone =~ s/\"//g; # remove quotes
    $feature_data_type =~ s/\"//g; # remove quotes
    print OUT "$line\n";

  # Feature RNASeq 5880 5881 2
  } elsif ($line =~ /^Feature\s+\S+\s+\d+\s+\d+/) {
    my $old_clone = $clone;
    @line = split /\s+/, $line;
    # $line[2] and $line[3] are the clone coords
    my $start = $line[2];
    my $end = $line[3];
    my ($new_clone_id, $indel, $change);
    ($new_clone_id, $start, $end, $indel, $change) = 
      $assembly_mapper->remap_clone($clone, $start, $end, $current_converter, $autoace_converter);

    if ($indel) {
      $log->write_to("There is an indel in the sequence in clone $new_clone_id, $start, $end\n");
    } elsif ($change) {
      $log->write_to("There is a change in the sequence in clone $new_clone_id, $start, $end\n");
    }
    
    # see if we have remapped to a different clone
    # if so, make a new line for "Feature_data : "Y44F5A:Polysome""
    if ($new_clone_id ne $clone) {
      print "Line <$line> for clone $clone has been remapped to $new_clone_id\n";
      print OUT "\n";
      print OUT "// this has been remapped from clone $clone to $new_clone_id\n";
      $clone = $new_clone_id;
      print OUT "Feature_data : \"${clone}:${feature_data_type}\"\n";

    }

    $line[2] = $start;
    $line[3] = $end;
    $line = join " ", @line;
    print OUT "$line\n";

    if ($old_clone ne $clone) {
      # continue with the original clone ID
      $clone = $old_clone;
      print OUT "\n";
      print OUT "Feature_data : \"${clone}:${feature_data_type}\"\n";
    }


  # blank line or comment
  } elsif ($line =~ /^\s*$/) {
    print OUT "$line\n";

  # die
  } else {
    $log->log_and_die("Unknown type of line found: $line\n");
    print "Unknown type of line found: $line\n";

  }


}

close (IN);
close (OUT);


# Close log files and exit
$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################

#
#$clone_length = &get_chrom_length($new_clone_id);
#
# Gets the length of a chrmoosome
#

sub get_chrom_length {
  my ($chrom) = @_;
  my $len;

  my $in = $wormbase->autoace . "/GFF_SPLITS/${chrom}_curated.gff";

  open (CHR, "< $in") || $log->log_and_die( "Can't open file $in\n");
  while (my $line = <CHR>) {
    if ($line =~ /##sequence-region\s+\S+\s+\S+\s+(\d+)/) {
	$len = $1;
    }
  }
  close (CHR);

  #print "$chrom length = $len\n";
  return $len;
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

=head2 NAME - remap_clone_homol_data.pl

=head1 USAGE

=over 4

=item remap_clone_homol_data.pl  [-options]

=back

This script reads in an ace file of Homol_data mapping something to clones and remaps the coordinates from the previous releases genome to the current genome.

script_template.pl MANDATORY arguments:

=over 4

=item -in file, input ace file from the previous release

=back

=over 4

=item -out file, output ace file for the current release

=back

=over 4

=item -data_type type, one of '21_urna', 'expression_pattern', 'mass_spec'

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
