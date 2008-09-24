#!/software/bin/perl -w
#
# intersect_gff.pl
# 
# by gw3               
#
# This is a script to fidn the strand-insensitive intersection of two GFF files.
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2008-09-24 12:50:23 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Modules::Overlap;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($species, $gff1, $gff2, $output);
my ($not_matching, $near_5, $near_3, $same_sense, $other_sense, $exact_match);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "gff1:s"     => \$gff1,
	    "gff2:s"     => \$gff2,
	    "output:s"   => \$output,
	    "not_matching" => \$not_matching,
	    "near_5:i"   => \$near_5,
	    "near_3:i"   => \$near_3,
	    "same_sense" => \$same_sense,
	    "other_sense"=> \$other_sense,
	    "exact_match"=> \$exact_match,
	    );

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

if (!defined $near_5) {$near_5 = 0}
if (!defined $near_3) {$near_3 = 0}
# sanity check
if (defined $same_sense && defined $other_sense && $same_sense == 1 && $other_sense == 1) {
  croak ("You can't choose for a match to only things on the same sense and only things on the opposite sense!\n");
}

##########################
# MAIN BODY OF SCRIPT
##########################

open (OUT, ">$output") || $log->log_and_die("Can't open $output: $!\n");

# loop through the chromosomes
foreach my $chromosome ($wormbase->get_chromosome_names(-mito => 1, -prefix => 1)) {

# load GFF files
  my @gff1_list = &read_gff_file($gff1, $chromosome);
  my @gff2_list = &read_gff_file($gff2, $chromosome);

# get overlaps and output
  my %state = (
	       last_used  => 0, # reset last secondary list's line
               last_used_forward => 0, # reset last secondary list's line for the forward sense
               last_used_reverse => 0, # reset last secondary list's line for the reverse sense
               near_5     => $near_5, # 0= we don't allow near 5' matches to count as a match
               near_3     => $near_3, # ditto for 3'
               same_sense => $same_sense, # 0= we want to allow a match to an object in either sense, set to true if want the same sense only
               other_sense => $other_sense, # 0= we want to allow a match to an object in either sense, set to true if want the opposite sense only
               exact_match => $exact_match, # 0= we do not want an exact match, but some sort of sloppy overlap
	       );

  foreach my $gff1_line (@gff1_list) {
    my $match_result = match($gff1_line, \@gff2_list, \%state);
    if ($match_result && !$not_matching) { # we want GFF1 lines that do have a match to GFF2
      my $line = join "\t", @{$gff1_line};
      print OUT "$line\n";
    } elsif (!$match_result && $not_matching) {	# we want GFF1 lines that don't have a match to GFF2
      my $line = join "\t", @{$gff1_line};
      print OUT "$line\n";
    }
  }
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
#  @gff1_list = &read_gff_file($gff1);

sub read_gff_file {

  my ($file, $chromosome) = @_;
  my @result;
  open (GFF, "grep \"^$chromosome\\W\" $file|") || $log->log_and_die("Can't open $file: $!\n");
  while (my $line = <GFF>) {
    chomp $line;
    if ($line =~ /^\s*$/) {next;}
    if ($line =~ /^#/) {next;}
    my @f = split /\t/, $line;
    push @result, [@f];
  }
  close (GFF);
  @result = sort {$a->[3] <=> $b->[3] or $a->[4] <=> $b->[4]} @result; # sort by start position
  return @result;
}

##########################################
#    if (match($gff1_line, \@gff2_list, \%state)) {

sub match {
  my ($main, $gff2_list_aref, $state_href) = @_;

  my $no_of_matches = 0;                # result of the match - no of matches found

  my $main_start = $main->[3];
  my $main_end = $main->[4];
  my $main_strand = $main->[6];
  #print "THIS LINE: @{$main}\n";

  my $near_3;
  my $near_5;
  if ($main_strand eq '+') {
    $near_3 = $state_href->{near_3};
    $near_5 = $state_href->{near_5};
  } else {                      # swap the values around in the reverse sense
    $near_3 = $state_href->{near_5};
    $near_5 = $state_href->{near_3};
  }

  @{$state_href->{matching_data}} = (); # no matching IDs found yet

  # if searching for same/opposite sense matches, then set last_used to be the
  # minimum of last_used_forward and last_used_reverse
  if ($state_href->{same_sense} || $state_href->{other_sense}) {
    $state_href->{last_used} = $state_href->{last_used_forward};
    if ($state_href->{last_used_reverse} <
        $state_href->{last_used_forward}) {
      $state_href->{last_used} = $state_href->{last_used_reverse};
    }
  }

  #print "last_used $state_href->{last_used}\n";

  for (my $i = $state_href->{last_used}; defined $gff2_list_aref->[$i]; $i++) {

    my $secondary = $gff2_list_aref->[$i];

    # secondary_start = $secondary->[3]
    # secondary_end = $secondary->[4]

    # test if there is an overlap (allowing possible nearby matches)
    # and optionally test if the senses are the same
    #print "SECONDARY LINE: @{$secondary}\n";

   if (
        (
         !$state_href->{exact_match} &&
         $secondary->[3] <= $main_end + $near_5 &&
         $secondary->[4] >= $main_start - $near_3 &&
         ($state_href->{same_sense}?($main_strand eq $secondary->[6]):1) &&
         ($state_href->{other_sense}?($main_strand ne $secondary->[6]):1)
         ) ||

        # for when we require an exact match
        (
         $state_href->{exact_match} &&
         $secondary->[3] == $main_start &&
         $secondary->[4] == $main_end &&
         ($state_href->{same_sense}?($main_strand eq $secondary->[6]):1) &&
         ($state_href->{other_sense}?($main_strand ne $secondary->[6]):1)
         )
        ) {

      # note that we have a match
      $no_of_matches++;

      # if we have not yet noted a match for this line, then remember
      # we got to here and found the first match
      if ($no_of_matches == 1) {
        # see if we are testing for them to be in the same/opposite sense
        if ($state_href->{same_sense} || $state_href->{other_sense}) {
          if ($main_strand eq '+') {
            # remember where we got up to in this sense
            $state_href->{last_used_forward} = $i;
          } else {
            $state_href->{last_used_reverse} = $i;
          }
        }
        $state_href->{last_used} = $i;
      }


      # save information about this match
      $state_href->{main_entry} = $main;
      push @{$state_href->{matching_data}}, $secondary; # the secondary_list entry

    } else {
      #print "no match\n";
      # don't search any further for this one if no overlap and the secondary_start > this_end
      if ($secondary->[3] > $main_end + $near_5) {last;}        # we have gone past the end of all possible matches
    }

  }
  #print "out of SECONDARY loop\n";
  # return the matching secondary_list entries
  # this is a count of the matches found when in scalar mode
  return @{$state_href->{matching_data}};

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

=head2 NAME - intersect_gff.pl

=head1 USAGE

=over 4

=item intersect_gff -gff1 file.gff -gff2 file2.gff -output out.gff  [-options]

=back

This script does a set intersection of the features in two GFF files looking for overlaps between any of the features in the first file as compared to the second and outputs any matching lines from the first file.

script_template.pl MANDATORY arguments:

=over 4

=item -gff1, the main GFF file. The lines in this file that match the other one will be output.

=back

=over 4

=item -gff2, the secondary GFF file. The first fiel will be compared to this to look for matches.

=back

=over 4

=item -output, the output file to hold the matching ilnes from the -gff1 file.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -species, if a species other than elegans is being used then specify its name. This is used to get a list of the chromosomes in the species.

=back

=over 4

=item -not_matching, if this is set that the lines in the -gff1 file that do not match any in the -gff2 file will be output instead of the ones that do match.

=back

=over 4

=item -near_5, this takes the number of bases to add on to a -gff2 feature to allow features in -gff1 that are near the gff2 feature's 5' end to be considered a match. This allows for near misses. If it is given a negative number then the gff1 feature must overlap the gff2 5' end by at least that number of bases. This ensures only good matches are allowed.

=back

=over 4

=item -near_3, this takes the number of bases to add on to a -gff2 feature to allow features in -gff1 that are near the gff2 feature's 3' end to be considered a match. This allows for near misses. If it is given a negative number then the gff1 feature must overlap the gff2 3' end by at least that number of bases. This ensures only good matches are allowed.

=back

=over 4

=item -other_sense, if set then only matches to a gff2 feature in the other sense to the gff1 feature will be considered to be a match. The default is to allow overlaps of features in either sense to be a match.

=back

=over 4

=item -same_sense, if set then only matches to a gff2 feature in the same sense to the gff1 feature will be considered to be a match. The default is to allow overlaps of features in either sense to be a match.

=back

=over 4

=item -exact_match, if set then the features must have the same start and stop positions to be a match. The default is for any overlap to be a match.

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


=head1 AUTHOR

=over 4

=item Gary Williams (gw3@sanger.ac.uk)

=back

=cut
