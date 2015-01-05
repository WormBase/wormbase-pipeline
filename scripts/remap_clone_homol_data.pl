#!/software/bin/perl -w
#
# remap_21_urna.pl                     
# 
# by Gary Williams                         
#
# This takes a BUILD_DATA/MISC_DYNAMIC/*.ace file mapped to clone
# Homol_data and converts any coordinates that have changed between
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
my ($input, $output, $data_type);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
            "input:s"    => \$input,
	    "output:s"   => \$output,
	    "data_type:s"  => \$data_type,
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


#Homol_data : Y57G11B:BLAT_ncRNA
#DNA_homol       EF044580        BLAT_ncRNA_BEST 100     42213   42233   1 21
#
#Sequence : Y57G11B
#Homol_data      Y57G11B:BLAT_ncRNA 1 48889
#
#Homol_data : Y51H4A:BLAT_ncRNA
#DNA_homol       EF044581        BLAT_ncRNA_BEST 100     163195  163215  1 21
#
#Sequence : Y51H4A
#Homol_data      Y51H4A:BLAT_ncRNA 1 269619
#
#Homol_data : Y41E3:BLAT_ncRNA
#DNA_homol       EF044582        BLAT_ncRNA_BEST 100     31396   31376   1 21


my $method; # as in: Homol_data : F56C11:<Expr>
my $homol_type; # as in: <Expr_homol> Expr1389 "Expr_pattern" 100.000000 14018 19205 1 5188


if ($data_type eq '21_urna') {
  $method = 'BLAT_ncRNA';
  $homol_type = 'DNA_homol';
} elsif ($data_type eq 'expression_pattern') {
  $method = 'Expr';
  $homol_type = 'Expr_homol';
} elsif ($data_type eq 'Tijsterman_G4') {
  $method = 'G4_DNA';
  $homol_type = 'Motif_homol';

} else {
  $log->log_and_die("Unknown data type: $data_type\n");
}


# start the coords converters
my $current_converter = Coords_converter->invoke($currentdb, 0, $wormbase);
my $autoace_converter = Coords_converter->invoke($ace_dir, 0, $wormbase);

# get the 21_URNA details
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

  # Homol_data : R02D3:BLAT_ncRNA
  if ($line =~ /^Homol_data\s+:\s+(\S+):/) { 
    print STDERR "$line\n" if ($verbose);
    $clone_id = $1;
    $clone_id =~ s/"//g; #"

  # DNA_homol       EF044580        BLAT_ncRNA_BEST 100     42213   42233   1 21
  } elsif ($line =~ /^${homol_type}\s+\S+\s+\S+\s+\S+\s+\d+\s+\d+\s+\d+\s+\d+/) {
    my ($homol, $id, $type, $score, $start, $end, $cb_start, $cb_end) = split(/\s+/, $line);
    $id =~ s/"//g; #"
    print STDERR "$line\n" if ($verbose);
    print STDERR " Remapping $clone_id, $start, $end\n" if $verbose;
    # if $start > $end, then sense is -ve (i.e. normal ace convention)
    ($new_clone_id, $start, $end, $indel, $change) = 
	$assembly_mapper->remap_clone($clone_id, $start, $end, $current_converter, $autoace_converter);

    if ($indel) {
      $log->write_to("There is an indel in the sequence in 21_URNA $id, clone $new_clone_id, $start, $end\n");
    } elsif ($change) {
      $log->write_to("There is a change in the sequence in 21_URNA $id, clone $new_clone_id, $start, $end\n");
    }

    # note this clone for later
    $clones_seen{$new_clone_id} = 1;

    if ($clone_id ne $prev_clone_id || $clone_id ne $new_clone_id) { # if the clone has changed, then write a new object header
      if ($blank_line == 0) {print OUT "\n";}
      print OUT "Homol_data : $new_clone_id:${method}\n";
      $prev_clone_id = $new_clone_id;
    }
    print OUT "${homol_type} $id $type $score $start $end $cb_start $cb_end\n";
    $blank_line = 0;

  # Homol_data R02D3:BLAT_ncRNA 1 32218
  } elsif ($line =~ /^Homol_data\s+\S+:${method}(\"|)\s+1/) {	# ignore these - write them out again later
    $prev_clone_id = "";

  } elsif ($line =~ /^Sequence\s+:\s+/) {	# ignore these - write them out again later
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
  print OUT "\nSequence : $clone_id\n";
  print OUT "Homol_data $clone_id:${method} 1 $clone_length\n";
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
