#!/usr/local/bin/perl5.8.0 -w
#
# remap_genefinder_between_releases.pl                     
# 
# by Gary Williams                         
#
# This takes the BUILD_DATA/MISC_DYNAMIC/waba.ace file and converts any coordinates that have changed between releases
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2006-05-16 15:51:12 $      

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
my @mapping_data = Remap_Sequence_Change::read_mapping_data($version - 1, $version);
 

##########################
# MAIN BODY OF SCRIPT
##########################


#
# Read the WABA locations from the previous ace file
# in order to be able to remap them
#


#Sequence : "2L52"
#Homol_data "2L52:waba" 1 4592
#                                                                                                               
#Homol_data : "2L52:waba"
#DNA_homol       "cb25.fpc0305"  WABA_coding     32.0000 1342    1361    264066  264085
#DNA_homol       "cb25.fpc0305"  WABA_weak       2.0000  1362    1382    264086  264106
#DNA_homol       "cb25.fpc0305"  WABA_strong     34.0000 1383    1399    264122  264138
#DNA_homol       "cb25.fpc0305"  WABA_weak       -8.0000 1400    1466    264139  264205
#DNA_homol       "cb25.fpc0305"  WABA_weak       -58.0000        1467    1661    264236  264428
#DNA_homol       "cb25.fpc0305"  WABA_weak       -72.0000        1728    1918    264429  264616
#DNA_homol       "cb25.fpc0305"  WABA_weak       -64.3200        2016    2357    264617  264951
#DNA_homol       "cb25.fpc0305"  WABA_weak       10.0000 2358    2400    265095  265137
#DNA_homol       "cb25.fpc0305"  WABA_coding     88.0000 2401    2508    265138  265245
#DNA_homol       "cb25.fpc0305"  WABA_weak       10.0000 2509    2539    265246  265276
#DNA_homol       "cb25.fpc0305"  WABA_weak       32.4500 2540    2902    265283  265652
#DNA_homol       "cb25.fpc0305"  WABA_coding     80.0000 3240    3311    265653  265724
#DNA_homol       "cb25.fpc0305"  WABA_weak       30.0000 3327    3393    265725  265791
#DNA_homol       "cb25.fpc0305"  WABA_weak       -98.0000        3597    3924    265792  266120
#DNA_homol       "cb25.fpc0305"  WABA_weak       8.0000  3925    4068    266148  266290
#DNA_homol       "cb25.fpc0305"  WABA_weak       -121.5500       4069    4436    266579  266937
#DNA_homol       "cb25.fpc0305"  WABA_coding     82.0000 4437    4499    266981  267043
                                                                                                               


# start the coords converters
my $current_converter = Coords_converter->invoke($currentdb, 0, $wormbase);
my $autoace_converter = Coords_converter->invoke($ace_dir, 0, $wormbase);

# get the WABA details
my $prev_clone_id = "";
my $clone_id;
my $new_clone_id;
my ($indel, $change);
my %clonesize       = $wormbase->FetchData('clonesize', "$ace_dir/COMMON_DATA"); 
my $clone_length;

open (IN, "< $input") || die "can't open input file $input\n";
open (OUT, "> $output") || die "can't open output file $output\n";
while (my $line = <IN>) {
  chomp $line;

  if ($line =~ /Sequence\s+:\s+\"(\S+)\"/) { 
    $clone_id = $1;

  } elsif ($line =~ /DNA_homol\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
    my ($waba_id, $waba_type, $waba_score, $start, $end, $cb_start, $cb_end) = ($1, $2, $3, $4, $5, $6, $7);
    print "$line\n" if ($verbose);

    # if $start > $end, then sense is -ve (i.e. normal ace convention)
    ($new_clone_id, $start, $end, $indel, $change) = 
	Remap_Sequence_Change::remap_clone($wormbase, $clone_id, $start, $end, $version, $current_converter, $autoace_converter, @mapping_data);

    if ($indel) {
      $log->write_to("There is an indel in the sequence in WABA $waba_id, clone $new_clone_id, $start, $end\n");
    } elsif ($change) {
      $log->write_to("There is a change in the sequence in WABA $waba_id, clone $new_clone_id, $start, $end\n");
    }

    # if the new clone is the same as the previous clone, continue to write out the next DNA_homol line
    # else we need to write new Sequence and Homol_data lines.
    if ($new_clone_id ne $prev_clone_id) {

      print OUT "\nSequence : \"$new_clone_id\"\n";

      # get the superlink or clone length
      if ($new_clone_id =~ /CHROMOSOME/) {
	$clone_length = &get_chrom_length($new_clone_id);
      } elsif ($new_clone_id =~ /SUPERLINK/) {
	$clone_length = $autoace_converter->Superlink_length($new_clone_id);
      } else {
	$clone_length = $clonesize{$new_clone_id};
      }

      print OUT "Homol_data \"$new_clone_id:waba\" 1 $clone_length\n\n";

      print OUT "Homol_data : \"$new_clone_id:waba\"\n";

      $prev_clone_id = $new_clone_id;
    }
    print OUT "DNA_homol $waba_id $waba_type $waba_score $start $end $cb_start $cb_end\n";

  } elsif ($line =~ /Homol_data/) {	# ignore the Homol_data lines

  } else {
    print OUT "$line\n";
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

  open (IN, "< $in") || die "Can't open file $in\n";
  while (my $line = <IN>) {
    if ($line =~ /##sequence-region\s+\S+\s+\S+\s+(\d+)/) {
	$len = $1;
    }
  }
  close ($in);

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
