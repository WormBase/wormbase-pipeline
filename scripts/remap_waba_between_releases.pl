#!/usr/local/bin/perl5.8.0 -w
#
# remap_genefinder_between_releases.pl                     
# 
# by Gary Williams                         
#
# This takes the BUILD_DATA/MISC_DYNAMIC/waba.ace file and converts any coordinates that have changed between releases
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2008-11-10 11:16:03 $      

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
#my $currentdb = $wormbase->database('current');

##########################
# read in the mapping data
##########################

my $version = $wormbase->get_wormbase_version;
print "Getting mapping data for WS$version\n";
my @mapping_data = Remap_Sequence_Change::read_mapping_data($version - 1, $version, $wormbase->species);
 

##########################
# MAIN BODY OF SCRIPT
##########################


#
# Read the WABA locations from the previous ace file
# in order to be able to remap them
#

#Sequence : "CHROMOSOME_III"
#Homol_data "CHROMOSOME_III:waba" 1 13783681
#
#Homol_data : "CHROMOSOME_III:waba"
#DNA_homol       "chrUn" WABA_strong     32.0000 8301    8322    4630632 4630653
#DNA_homol       "chrUn" WABA_weak       18.0000 8323    8419    4630654 4630750
#DNA_homol       "chrUn" WABA_coding     28.0000 8328    8353    1830139 1830164
#DNA_homol       "chrUn" WABA_weak       -6.0000 8354    8384    1830200 1830230
#DNA_homol       "chrUn" WABA_coding     36.0000 8385    8419    1830231 1830264
#DNA_homol       "chrUn" WABA_weak       28.0000 8420    8481    1830438 1830497
#DNA_homol       "chrUn" WABA_weak       -40.0000        8420    8503    4630758 4630842

# and

#Sequence : "chrUn"
#Homol_data "chrUn:waba" 1 7311690
#
#Homol_data : "chrUn:waba"
#DNA_homol       "CHROMOSOME_II" WABA_weak       -80.0000        6586954 6586674 15408   15682
#DNA_homol       "CHROMOSOME_II" WABA_weak       -10.0000        6586673 6586516 15846   16007
#DNA_homol       "CHROMOSOME_II" WABA_weak       -64.7700        6586515 6586300 16021   16236
#DNA_homol       "CHROMOSOME_II" WABA_weak       -132.3200       6586257 6585788 16237   16699
#DNA_homol       "CHROMOSOME_II" WABA_weak       32.0000 6585766 6585589 16700   16877
#DNA_homol       "CHROMOSOME_II" WABA_strong     50.0000 6585588 6585556 16878   16910

                                                                                                               


# start the coords converters
my $autoace_converter = Coords_converter->invoke($ace_dir, 0, $wormbase);

# get the WABA details
my $prev_clone_id = "";
my $chromosome;
my $chrom_length;
my ($indel, $change);

open (IN, "< $input") || die "can't open input file $input\n";
open (OUT, "> $output") || die "can't open output file $output\n";
while (my $line = <IN>) {
  print "$line\n" if ($verbose);
  chomp $line;

  if ($line =~ /Sequence\s+:\s+\"(\S+)\"/) { 
    $chromosome = $1;

    if ($chromosome =~ /CHROMOSOME/) { # this is an elegans chromosome
      $chrom_length = &get_chrom_length($chromosome);
    }
    print OUT "$line\n";

  } elsif ($line =~ /Homol_data\s+\"/) {
    if ($chromosome =~ /CHROMOSOME/) { # this is an elegans chromosome
      print OUT "Homol_data \"$chromosome:waba\" 1 $chrom_length\n";
    } else {
      print OUT "$line\n";
    }

  } elsif ($line =~ /DNA_homol\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
    my ($waba_id, $waba_type, $waba_score, $start, $end, $cb_start, $cb_end) = ($1, $2, $3, $4, $5, $6, $7);

    if ($chromosome =~ /CHROMOSOME/) { # this is an elegans chromosome

      # if $start > $end, then sense is -ve (i.e. normal ace convention)
      ($start, $end, $indel, $change) = 
	  Remap_Sequence_Change::remap_ace($chromosome, $start, $end, $version-1, $version, @mapping_data);

      if ($indel) {
	$log->write_to("There is an indel in the sequence in WABA $waba_id, chromosome $chromosome, $start, $end\n");
      } elsif ($change) {
	$log->write_to("There is a change in the sequence in WABA $waba_id, chromosome $chromosome, $start, $end\n");
      }

      print OUT "DNA_homol $waba_id $waba_type $waba_score $start $end $cb_start $cb_end\n";
    } else {			# briggsae chromosome
      my $elegans_chrom = $waba_id;
      $elegans_chrom =~ s/\"//g;
      ($cb_start, $cb_end, $indel, $change) = 
	  Remap_Sequence_Change::remap_ace($elegans_chrom, $cb_start, $cb_end, $version-1, $version, @mapping_data);

      if ($indel) {
	$log->write_to("There is an indel in the sequence in WABA $waba_id, chromosome $chromosome, $start, $end\n");
      } elsif ($change) {
	$log->write_to("There is a change in the sequence in WABA $waba_id, chromosome $chromosome, $start, $end\n");
      }

      print OUT "DNA_homol $waba_id $waba_type $waba_score $start $end $cb_start $cb_end\n";
      
    }
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

  open (CHROM, "< $in") || die "Can't open file $in\n";
  while (my $line = <CHROM>) {
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
