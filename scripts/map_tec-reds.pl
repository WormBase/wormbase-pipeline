#!/software/bin/perl -w

###########################################################################################################################
#
#
# This writes the MISC_DYNAMIC files misc_TEC_RED-homol.ace and misc_TEC_RED_homol_data.ace
#
###########################################################################################################################

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Log_files;
use Storable;
use Ace;
use Coords_converter;
use Getopt::Long;

my ( $help, $debug, $test, $store );
my $verbose;    # for toggling extra output


##############################
# command-line options       #
##############################

GetOptions(
    "help"      => \$help,
    "debug=s"   => \$debug,
    "test"      => \$test,
    "store:s"   => \$store,
);

$test = 1;
$debug = "gw3";


# Display help if required
&usage("Help") if ($help);

my $wormbase;

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             );
}


my $log = Log_files->make_build_log($wormbase);

my $currentdb = $wormbase->database('current');	# use autoace when debugged 
my $ace_dir = $wormbase->autoace;     # AUTOACE DATABASE DIR
#my $database_path = $currentdb;     # full path to local AceDB database; change as appropriate
my $database_path = $ace_dir;     # full path to local AceDB database; change as appropriate
my $program = $wormbase->tace;  # full path to tace; change as appropriate

print "connecting to server...";
my $db = Ace->connect(-path => $database_path,  -program => $program) || die print "Connection failure: ", Ace->error;  # local database
print "done\n";
my $coords = Coords_converter->invoke($database_path, 0, $wormbase);
my %clonesize = $wormbase->FetchData('clonesize');


my $query="find Feature\nwhere Method = \"SL*\"\nshow -a\nquit\n";

my $tace = $wormbase->tace;        # TACE PATH

my $clone;
my $tec_red;
my $prev_tec_red = "";		# to stop processing useless duplicates
my $feature;
my %homol_data;

my $output = $wormbase->misc_dynamic."/misc_TEC_RED_homol.ace";
open (OUT, "> $output") || die "Can't open $output\n";
open (TACE, "echo '$query' | $tace $database_path |");
while (my $line = <TACE>) {
  next if ($line =~ /acedb\>/);
  next if ($line =~ /\/\//);
  chomp $line;

  if ($line =~ /Feature\s+:\s+\"(\S+)\"/) {
    $feature = $1;
    $prev_tec_red = "";

  } elsif ($line =~ /Sequence\s+\"(\S+)\"/) {
    $clone = $1;
    #print "Next clone: $clone\n";
  } elsif ($line =~ /Defined_by_sequence\s+\"(\S+)\"/) {
    $tec_red = $1;
    if ($prev_tec_red ne $tec_red) {
      &map_tec_red($clone, $tec_red, $feature, \%homol_data);
    }
    $prev_tec_red = $tec_red;
  }
}
close TACE;
close(OUT);

# write out the Sequence object's Homol_data
my $output2 = $wormbase->misc_dynamic."/misc_TEC_RED_homol_data.ace";
open (OUT2, "> $output2") || die "Can't open $output\n";
foreach my $clone (keys %homol_data) {
  my $clonelen = &get_clone_len($clone);
  print OUT2 "\nSequence : \"$clone\"\n";
  print OUT2 "Homol_data \"$clone:TEC_RED\" 1 $clonelen\n";
}
close(OUT2);

# get a clean exit status by closing and undef'ing things
$db->close;
$log->mail();
$wormbase = undef;
$log = undef;
exit(0);

##########################################
# map the tec-red to the position of its feature
# &map_tec_red($clone, $tec_red, $feature);

sub map_tec_red {
  my ($clone, $tec_red, $feature, $homol_data_href) = @_;

  # check that the sequence is a tec-red
  my $tec_red_obj = $db->fetch(Sequence => $tec_red);
  my $tag = $tec_red_obj->fetch('TSL_tag');
  if (! defined $tag) {return;}

  # get the length of the tec-red tag
  my $tec_red_len = $tec_red_obj->at('DNA[2]');

  # get start/end position of feature in clone
  my $clone_obj = $db->fetch(Sequence => $clone);
  my $start = $clone_obj->at("SMap.S_child.Feature_object.$feature"."[1]");
  my $end = $clone_obj->at("SMap.S_child.Feature_object.$feature"."[2]");
  
  # get the clone location of the tec_red 
  my $tec_red_start;
  my $tec_red_end;
  if ($start < $end) {
    $tec_red_start = $start - 1;
    $tec_red_end = $tec_red_start + $tec_red_len - 1;
  } else {			# reverse sense
    $tec_red_end = $start + 1;
    $tec_red_start = $tec_red_end - $tec_red_len + 1;
  }

  # get the chromosomal location of the tec_red
  my $start_chrom_coord;
  my $end_chrom_coord;
  my $chrom;
  ($chrom, $start_chrom_coord) = $coords->Coords_2chrom_coords($clone,  $tec_red_start);
  ($chrom, $end_chrom_coord) = $coords->Coords_2chrom_coords($clone,  $tec_red_end);

  # get the clone that the tec_red is on
  my @clone_coords = $coords->LocateSpan($chrom, $start_chrom_coord, $end_chrom_coord);
  $clone = $clone_coords[0];
  $tec_red_start = $clone_coords[1];
  $tec_red_end = $clone_coords[2];

  # write the Homol_data object's DNA_homol
  print OUT "\nHomol_data : \"$clone:TEC_RED\"\n";
  print OUT "DNA_homol \"$tec_red\" TEC_RED 100 $tec_red_start $tec_red_end ";
  if ($start < $end) {
    print OUT "1 $tec_red_len\n";
  } else {
    print OUT "$tec_red_len 1\n";
  }

  # store the Sequence object's Homol_data
  $homol_data_href->{$clone} = 1;

}

##########################################
# get the length of a clone/superlink/chromosome
# $len = get_clone_len($clone)

sub get_clone_len {
  my ($clone) = @_;

  my $clonelen = $clonesize{$clone};

  if (! defined $clonelen) {
    if ($clone =~ /SUPERLINK/ || $clone =~ /CHROMOSOME/) {
      # get the Superlink lengths from the Coords_converter data
      $clonelen = $coords->Superlink_length($clone);
    } else {
      die "undef returned for length of $clone\n";
    }
  }

  return $clonelen;

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

=pod

=head2 NAME - map_tec-reds.pl

=head1 USAGE

=over 4

=item map_tec-reds.pl  [-options]

=back

This script maps the TEC-RED sequences back to the genome after a
genomic sequence change by using the position of the Features that
were originally defined using these TEC-REDs.

This is using the fact that the Features will have been mapped first
by using their flanking sequences and to will be resistent to errors
caused by genomic changes.

script_template.pl MANDATORY arguments:

=over 4

=item 

=back

map_tec-reds.pl  OPTIONAL arguments:

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

=item This requires the autoace database and GFF files to read the transcript data from and to read the SAGE data from.

=back

=head1 AUTHOR

=over 4

=item Gary Williams

=back

=cut
