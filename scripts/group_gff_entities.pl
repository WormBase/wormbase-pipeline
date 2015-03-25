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
# Last updated on: $Date: 2015-03-25 14:52:47 $      

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
my ($input, $output, $est, $blastx);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "input:s"    => \$input,
            "output:s"   => \$output,
	    "est"        => \$est,
	    "blastx"     => \$blastx
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

my $parenttype;
if ($est) {
  $parenttype = 'cDNA_match';
} else {
  $parenttype = 'match';
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

my $prev_seq = '';
my $prev_target_start = 0;
my $prev_start = 0;
my $prev_end = 0;
my $prev_score = 0;
my $prev_sense = '';
my $prev_target = '';
my $parent = '';

my $seq;
my $source;
my $type;
my $start;
my $end;
my $score;
my $sense;
my $id;
my $target; 
my $target_start=0;
my $target_end;

my %seen;

my @input;
my @a;
open (IN, "<$input") || $log->log_and_die("cant open input file: $input\n");
while (my $line = <IN>) {
  @a = split /\t/, $line;
  if (!defined $a[8]) {next}
  $seq   = $a[0];
  $source = $a[1];
  $type  = $a[2];
  $start = $a[3];
  $end   = $a[4];
  $score = $a[5];
  $sense = $a[6];
  ($id) = ($a[8] =~ /^ID=(\S+?);/);
  ($target, $target_start, $target_end) = ($a[8] =~ /Target=(\S+)\s+(\d+)\s+(\d+)\s+/);
  push @input, [$seq, $source, $type, $start, $end, $score, $sense, $id, $target, $target_start, $target_end];
}
close(IN);

# sort by chromosome, sense, Target then position
my @sorted = sort {
  ($a->[0] cmp $b->[0]) || # seq
    ($a->[6] cmp $b->[6]) || # sense
      ($a->[8] cmp $b->[8]) || # Target
	($a->[3] <=> $b->[3]) # chromosome start
}  @input;



open (OUT, ">$output") || $log->log_and_die("cant open output file: $output\n");

# want to make a fresh parent each time there is a discontinuity in
# the target_start/target_end
foreach my $a (@sorted) {
  $seq   = $a->[0];
  $source = $a->[1];
  $type  = $a->[2];
  $start = $a->[3];

  $prev_end = $end;

  $end   = $a->[4];
  $score = $a->[5];
  $sense = $a->[6];
  $id = $a->[7];

  $prev_target_start = $target_start;

  ($target, $target_start, $target_end) = ($a->[8], $a->[9], $a->[10]);

  if (($target ne $prev_target) || 
      ($sense eq '+' && $target_start < $prev_target_start) || # discontinuity in the target alignment
      ($sense eq '-' && $target_start > $prev_target_start) || 
      ($start > $prev_end + 20000) # too large an intron
     ) {

    if ($prev_target ne '') {
      print OUT "$seq\t$source\t$parenttype\t$prev_start\t$prev_end\t$prev_score\t$prev_sense\t.\tID=$parent\n";
    }
    $prev_target = $target;
    $prev_start = $start;
    $prev_end = $end;
    $prev_score = $score;
    $prev_sense = $sense;

    # make a unique name for the parent
    my $parent_number='';
    if (exists $seen{$target}) {
      $parent_number = ".$seen{$target}";
    }

    $seen{$target}++;

    if ($est) {
      $parent = "EST:" . $target . $parent_number;
    } elsif ($blastx) {
      $parent = "BLASTX:" . $target . $parent_number;
    } else {
      $log->log_and_die("-est or -blastx should be specified\n");
    }
  }
  $prev_end = $end;

  $a[8] = 'Parent=' . $parent . ';' . $a[8];
  
  my $idstr='';
  if (defined $id) {$idstr = "ID=${id}"}
  my $new_line = "$seq\t$source\t$type\t$start\t$end\t$score\t$sense\t.\tParent=${parent};${idstr}Target=$target $target_start $target_end +\n";
  print OUT $new_line;
}

# last line
if ($prev_target ne '') {
  print OUT "$seq\t$source\t$parenttype\t$prev_start\t$prev_end\t$prev_score\t$prev_sense\t.\tID=$parent\n";
}

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
