#!/usr/local/bin/perl5.8.0 -w
#
# fix_interpolation.pl                         
# 
# by Keith Bradnam                         
#
# This fixes the gene interpolation paoitions.
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2007-04-12 09:16:59 $      

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


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"      => \$store,
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
my $tace            = $wormbase->tace;        # TACE PATH




##########################
# MAIN BODY OF SCRIPT
##########################

my @genes;
my $input = "/nfs/disk100/wormpub/BUILD/autoace/logs/rev_phys.log";
#my $input = "rev_phys.log";
open(IN, "< $input") || die "Can't open $input\n";
while (my $line = <IN>) {
  if ($line =~ /\-\-\-\-/) {next;}
  my @columns = split /\s+/, $line;
  #print "@columns[0..3]\n";
  push @genes, [(@columns[0..3], $columns[2])]; # add a duplicate of the position on the end so we know which genes we changed
}
close(IN);

my $changed_in_this_iteration;
do {
  $changed_in_this_iteration = 0;
  my $prev_chrom = "";
  my $prev_pos;
  my $next_pos;
  for (my $i = 0; $i < @genes; $i++) {
    my $chrom = $genes[$i]->[0];
    my $pos = $genes[$i]->[2];
    if ($prev_chrom ne $chrom) {
      $prev_pos = $pos;
      $prev_chrom = $chrom;
      next;			# always skip the first gene in the chromosome
    }
#    print "\t$chrom $prev_pos $pos\n";
    
    # get the next position
    # are we at the end of the array or end of the chromosome?
    if ($i + 1 < @genes && $genes[$i+1]->[0] eq $chrom) {
      $next_pos = $genes[$i+1]->[2];
    } else {
      $next_pos = $pos + 0.5;
    }

    # should this position be changed? Test for this position less than previous.
    if ($prev_pos > $pos) {
#      print "Found a position to be changed: $genes[$i][3]\n";

      # get the difference between the previous and next positions
      my $diff = $next_pos - $prev_pos;
      if ($diff > 0.0005) {
	$genes[$i]->[2] = $prev_pos + ($diff/2);
      } else {
	$genes[$i]->[2] = $prev_pos + 0.0005;
      }
      $pos = $genes[$i]->[2];
      $changed_in_this_iteration = 1;
      #print "\t new position: $pos\n";
    }

    
    # should this position be changed? Test for this position greater than next
    if ($pos > $next_pos && $next_pos > $prev_pos) {
      #print "Found a position to be changed: $genes[$i][3]\n";

      # get the difference between the previous and next positions
      my $diff = $next_pos - $prev_pos;
      if ($diff > 0.0005) {
	$genes[$i]->[2] = $prev_pos + ($diff/2);
      } else {
	$genes[$i]->[2] = $prev_pos + 0.0005;
      }
      $pos = $genes[$i]->[2];
      $changed_in_this_iteration = 1;
      #print "\t new position: $pos\n";
    }

    
    
    
    $prev_pos = $pos;
  }
} while ($changed_in_this_iteration);


# find the genes we changed and write out the ace file
my $out = "/tmp/interp.ace";
open (OUT, "> $out") || die "Can't open $out";
foreach my $columns (@genes) {
  if ($columns->[2] == $columns->[4]) {next;} # not changed
  my $chrom = $columns->[0];
  my $pos = $columns->[2];
  my $gene = $columns->[3];
  $chrom =~ s/://;		# remove the colon

  print OUT "\n";
  print OUT "Gene : $gene\n";
  print OUT "Map $chrom Position $pos\n";
  print  "\n";
  print  "Gene : $gene\n";
  print  "Map $chrom Position $pos\n";

}
close (OUT);

# load the ace file
print "\nLoad the file $out\n";

# Close log files and exit
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

=item Keith Bradnam (krb@sanger.ac.uk)

=back

=cut
