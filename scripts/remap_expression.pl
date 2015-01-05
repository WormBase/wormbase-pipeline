#!/usr/local/bin/perl5.8.0 -w
#
# remap_expression_between_releases.pl                     
# 
# by Gary Williams                         
#
# This takes a BUILD_DATA/MISC_DYNAMIC/*expression*.ace file and converts any coordinates that have changed between releases
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
#use Ace;
#use Sequence_extract;
#use Coords_converter;

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
# Read the Expression locations from the previous ace file
# in order to be able to remap them
#


#Homol_data : "AC3:Expr"
#Expr_homol        "Expr6346"  "Expr_pattern"  100.000000  1099  2072  1  973
#Expr_homol        "Expr2206"  "Expr_pattern"  100.000000  20039  18740  1  1300
#
#Homol_data : "AC7:Expr"
#Expr_homol        "Expr5001"  "Expr_pattern"  100.000000  13685  10884  1  2801
#Expr_homol        "Expr5002"  "Expr_pattern"  100.000000  25559  22653  1  2906
#Expr_homol        "Expr5003"  "Expr_pattern"  100.000000  25559  22653  1  2906
#
#Homol_data : "AH6:Expr"
#Expr_homol        "Expr5006"  "Expr_pattern"  100.000000  29227  32102  1  2875
#Expr_homol        "Expr5007"  "Expr_pattern"  100.000000  13370  16228  1  2858
#Expr_homol        "Expr5008"  "Expr_pattern"  100.000000  33258  30331  1  2927
#
#Sequence : "ZK1127"
#Homol_data       "ZK1127:Expr" 1 35962
#
#Sequence : "ZK1128"
#Homol_data       "ZK1128:Expr" 1 27832
#
#Sequence : "ZK1193"
#Homol_data       "ZK1193:Expr" 1 33075
                                                                                                                                      

# start the coords converters
my $current_converter = Coords_converter->invoke($currentdb, 0, $wormbase);
my $autoace_converter = Coords_converter->invoke($ace_dir, 0, $wormbase);

# get the EXPRESSION details
my $prev_clone_id = "";
my $clone_id;
my $new_clone_id;
my ($indel, $change);
my %clonesize       = $wormbase->FetchData('clonesize'); 
my $clone_length;

my %clones_seen;		# for writing out the Homol_data lines at the end
my $blank_line = 0;

open (IN, "< $input") || die "can't open input file $input\n";
open (OUT, "> $output") || die "can't open output file $output\n";
while (my $line = <IN>) {
  chomp $line;

  # Homol_data : "R02D3:Expr"
  if ($line =~ /Homol_data\s+:\s+\"(\S+):/) { 
    $clone_id = $1;

  # Expr_homol "Expr5008" "Expr_pattern" 100.0 12386 12345 1 2858
  } elsif ($line =~ /Expr_homol\s+\S+\s+\S+\s+\S+\s+\d+\s+\d+\s+\d+\s+\d+/) {
    my ($homol, $expression_id, $expression_type, $expression_score, $start, $end, $cb_start, $cb_end) = split(/\s+/, $line);
    $expression_id =~ s/"//g; #"
    print "$line\n" if ($verbose);
    # if $start > $end, then sense is -ve (i.e. normal ace convention)
    ($new_clone_id, $start, $end, $indel, $change) = 
	$assembly_mapper->remap_clone($clone_id, $start, $end, $current_converter, $autoace_converter);

    if ($indel) {
      $log->write_to("There is an indel in the sequence in EXPRESSION $expression_id, clone $new_clone_id, $start, $end\n");
    } elsif ($change) {
      $log->write_to("There is a change in the sequence in EXPRESSION $expression_id, clone $new_clone_id, $start, $end\n");
    }

    # note this clone for later
    $clones_seen{$new_clone_id} = 1;

    if ($clone_id ne $prev_clone_id || $clone_id ne $new_clone_id) { # if the clone has changed, then write a new object header
      if ($blank_line == 0) {print OUT "\n";}
      print OUT "Homol_data : \"$new_clone_id:Expr\" \n";
      $prev_clone_id = $new_clone_id;
    }
    print OUT "Expr_homol \"$expression_id\" $expression_type $expression_score $start $end $cb_start $cb_end \n";
    $blank_line = 0;

  # Homol_data R02D3:Expr 1 32218
  } elsif ($line =~ /Homol_data\s+\"(\S+):Expr\"\s+1/) {	# ignore these - write them out again later
    $prev_clone_id = "";

  } elsif ($line =~ /Sequence\s+:\s+/) {	# ignore these - write them out again later
    $prev_clone_id = "";

  } else {
    $prev_clone_id = "";
    if ($line =~ /^\s*$/) {
      if ($blank_line) {next;}	# collapse multiple blank lines into one.
      $blank_line=1;
    } else {
      $blank_line=0;
    }
    print OUT "$line\n";
  }

}

# now write out the Homol_data lines for all of the clones seen
foreach my $clone_id (keys %clones_seen) {
  # get the superlink or clone length
  if ($clone_id =~ /CHROMOSOME/) {
    $clone_length = &get_chrom_length($clone_id);
  } else {
    $clone_length = $clonesize{$clone_id};
  }
  print OUT "\nSequence : \"$clone_id\" \n";
  print OUT "Homol_data \"$clone_id:Expr\" 1 $clone_length \n";
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

  open (CHR, "< $in") || die "Can't open file $in\n";
  while (my $line = <CHR>) {
    if ($line =~ /##sequence-region\s+\S+\s+\S+\s+(\d+)/) {
	$len = $1;
    }
  }
  close (CHR);

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

=item Keith Bradnam (krb@sanger.ac.uk)

=back

=cut
