#!/software/bin/perl -w
#
# group_gff_entities.pl
# 
# by Gary Williams
#
# For use when a set of GFF lies should be grouped in order to display
# then better. e.g. EST expressed_sequence_match lines should be
# grouped with cDNA_match lines to produce the full-length alignment
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2015-03-25 13:50:18 $      

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
my ($input, $output);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "input:s"    => \$input,
            "output:s"   => \$output,
	    );

$debug = "gw3";

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

my $prev_seq = '';
my $prev_target_end = 0;
my $prev_start = 0;
my $prev_end = 0;
my $prev_score = 0;
my $prev_sense = '';
my $prev_id = '';
my $parent = '';

open (IN, "<$input") || $log->log_and_die("cant open input file: $input\n");
open (OUT, ">$output") || $log->log_and_die("cant open output file: $output\n");
while (my $line = <IN>) {
  my @a = split /\t/, $line;
  if (!defined $a[8]) {next}
  my $seq   = $a[0];
  my $start = $a[3];
  my $end   = $a[4];
  my $score = $a[5];
  my $sense = $a[6];
  my ($id) = ($a[8] =~ /^ID=(\S+?);/);
  my ($target_start, $target_end) = ($a[8] =~ /Target=\S+\s+(\d+)\s+(\d+)\s++/);

# want to make a fresh parent each time there is a discontinuity in
# the target_start/target_end
  if ($id ne $prev_id) {
    if ($prev_id ne '') {
      print OUT "$a[0]\t$a[1]\tcDNA_match\t$prev_start\t$prev_end\t$prev_score\t$prev_sense\t.\tID=$parent\n";
    }
    $prev_id = $id;
    $prev_start = $start;
    $prev_end = $end;
    $prev_score = $score;
    $prev_sense = $sense;
    $parent = "EST:" . $id;    
  }
  $prev_end = $end;

  $a[8] = 'Parent=' . $parent . ';' . $a[8];

  my $new_line = join "\t", @a;
  print OUT $new_line;
}
close(IN);
close(OUT);


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

=item xxx (xxx@sanger.ac.uk)

=back

=cut
