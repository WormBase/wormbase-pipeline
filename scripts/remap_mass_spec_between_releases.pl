#!/usr/local/bin/perl5.8.0 -w
#
# remap_mass_spec_between_releases.pl                     
# 
# by Gary Williams                         
#
# This takes a BUILD_DATA/MISC_DYNAMIC/*mass_spec*.ace file and converts any coordinates that have changed between releases
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2006-12-14 14:23:45 $      

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
# Read the Mass Spec locations from the previous ace file
# in order to be able to remap them
#


#Homol_data : "R02D3:Mass-spec"
#MSPeptide_homol "MSP:FVENCNETFDVTIR" mass_spec_genome 100.0 12386 12345 1 14
# 
#Mass_spec_peptide : "MSP:FVENCNETFDVTIR"
#Mass_spec_experiments "MH_exp_006" Protein_probability 1.00
#Mass_spec_experiments "MH_exp_006" Peptide_probability 1.00
#Mass_spec_experiments "MH_exp_006" Protein "WP:CE18098"
# 
#Protein : "MSP:FVENCNETFDVTIR"
#Peptide "MSP:FVENCNETFDVTIR"
# 
#Protein : "WP:CE18098"
#Pep_homol "MSP:FVENCNETFDVTIR" mass-spec 1 18 31 1 14
#
#Mass_spec_experiment : "MH_exp_006"
#Remark "The minimum peptide sequence is not directly given by the length rather by the mass of the peptide. The precursor ion scan window is 400-2000 m/z. This means a plus one charged peptide can be minimally 399 Da. A Trp, Arg plus a third aa will do it in principle."
#Person "WBPerson5305"
#Species "Caenorhabditis elegans"
#Strain "N2"
#Life_stage "L1"
#Sub_cellular_localization "membrane proteins, cystein-containing peptides"
#Digestion "Trypsin"
#Ionisation_source "ESI"
#Database "Wormpep proteins from release WS120"
#Program "SEQUEST"
#
#Sequence : R02D3
#Homol_data R02D3:Mass-spec 1 32218
                                                                                                                                              

# start the coords converters
my $current_converter = Coords_converter->invoke($currentdb, 0, $wormbase);
my $autoace_converter = Coords_converter->invoke($ace_dir, 0, $wormbase);

# get the MASS SPEC details
my $prev_clone_id = "";
my $clone_id;
my $new_clone_id;
my ($indel, $change);
my %clonesize       = $wormbase->FetchData('clonesize', "$ace_dir/COMMON_DATA"); 
my $clone_length;

my %clones_seen;		# for writing out the Homol_data lines at the end
my $blank_line = 0;

open (IN, "< $input") || die "can't open input file $input\n";
open (OUT, "> $output") || die "can't open output file $output\n";
while (my $line = <IN>) {
  chomp $line;

  # Homol_data : "R02D3:Mass-spec"
  if ($line =~ /Homol_data\s+:\s+\"(\S+):/) { 
    $clone_id = $1;

  # MSPeptide_homol "MSP:FVENCNETFDVTIR" mass_spec_genome 100.0 12386 12345 1 14
  } elsif ($line =~ /MSPeptide_homol\s+\"(\S+)\"\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
    my ($mass_spec_id, $mass_spec_type, $mass_spec_score, $start, $end, $cb_start, $cb_end) = ($1, $2, $3, $4, $5, $6, $7);
    print "$line\n" if ($verbose);

    # if $start > $end, then sense is -ve (i.e. normal ace convention)
    ($new_clone_id, $start, $end, $indel, $change) = 
	Remap_Sequence_Change::remap_clone($wormbase, $clone_id, $start, $end, $version, $current_converter, $autoace_converter, @mapping_data);

    if ($indel) {
      $log->write_to("There is an indel in the sequence in MASS SPEC $mass_spec_id, clone $new_clone_id, $start, $end\n");
    } elsif ($change) {
      $log->write_to("There is a change in the sequence in MASS SPEC $mass_spec_id, clone $new_clone_id, $start, $end\n");
    }

    # note this clone for later
    $clones_seen{$new_clone_id} = 1;

    if ($clone_id ne $prev_clone_id || $clone_id ne $new_clone_id) { # if the clone has changed, then write a new object header
      if ($blank_line == 0) {print OUT "\n";}
      print OUT "Homol_data : \"$new_clone_id:Mass-spec\"\n";
      $prev_clone_id = $new_clone_id;
    }
    print OUT "MSPeptide_homol \"$mass_spec_id\" $mass_spec_type $mass_spec_score $start $end $cb_start $cb_end\n";
    $blank_line = 0;

  # Homol_data R02D3:Mass-spec 1 32218
  } elsif ($line =~ /Homol_data\s+(\S+):Mass-spec\s+1\s+\d+/) {	# ignore these - write them out again later
    $prev_clone_id = "";

  } elsif ($line =~ /Sequence\s+:\s+/) {	# ignore these - write them out again later
    $prev_clone_id = "";

  } else {
    print OUT "$line\n";
    $prev_clone_id = "";
    if ($line =~ /^\s*$/) {
      $blank_line=1;
    } else {
      $blank_line=0;
    }
  }

}

# now write out the Homol_data lines for all of the clones seen
foreach my $clone_id (keys %clones_seen) {
  # get the superlink or clone length
  if ($clone_id =~ /CHROMOSOME/) {
    $clone_length = &get_chrom_length($clone_id);
  } elsif ($clone_id =~ /SUPERLINK/) {
    $clone_length = $autoace_converter->Superlink_length($clone_id);
  } else {
    $clone_length = $clonesize{$clone_id};
  }
  print OUT "\nSequence : $clone_id\n";
  print OUT "Homol_data $clone_id:Mass-spec 1 $clone_length\n";
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
