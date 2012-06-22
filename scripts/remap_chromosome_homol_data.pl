#!/software/bin/perl -w
#
# remap_chromosome_homol_data.pl
# 
# by Gary Williams                         
#
# This takes a BUILD_DATA/MISC_DYNAMIC/*.ace file mapped to chromosome
# Homol_data in ~50Kb virtual blocks and converts any coordinates that have changed between
# releases
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2012-06-22 08:56:53 $      

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

my $assembly_mapper = Remap_Sequence_Change->new($version -1, $version, $wormbase->species, $wormbase->genome_diffs);

##########################
# MAIN BODY OF SCRIPT
##########################


#
# Read the locations from the previous ace file
# in order to be able to remap them
#


# Homol_data : "RNASeq_Hillier.Adult_spe-9:CHROMOSOME_III_1"
# DNA_homol "Adult_spe-9_23dC_8days_post-L4_molt_bundle_of_reads_supporting_intron_III_10073_10364_+_wb170" RNASeq_Hillier.Adult_spe-9 1.0 10047 10072 1 26
# DNA_homol "Adult_spe-9_23dC_8days_post-L4_molt_bundle_of_reads_supporting_intron_III_10073_10364_+_wb170" RNASeq_Hillier.Adult_spe-9 1.0 10365 10374 27 36
# DNA_homol "Adult_spe-9_23dC_8days_post-L4_molt_bundle_of_reads_supporting_intron_III_12402_12475_-_wb170" RNASeq_Hillier.Adult_spe-9 8.5 12368 12401 68 35
# 
# Sequence : "CHROMOSOME_I"
# Confirmed_intron 10015132 10015184 EST Adult_spe-9_23dC_8days_post-L4_molt_bundle_of_reads_supporting_intron_I_10015130_10015182_+_wb170
# Confirmed_intron 10016971 10017028 EST Adult_spe-9_23dC_8days_post-L4_molt_bundle_of_reads_supporting_intron_I_10016969_10017026_+_wb170
# Confirmed_intron 10017131 10017179 EST Adult_spe-9_23dC_8days_post-L4_molt_bundle_of_reads_supporting_intron_I_10017129_10017177_+_wb170
# S_child Homol_data RNASeq_Hillier.Adult_spe-9:CHROMOSOME_I_1 1 150000
# S_child Feature_data Confirmed_intron_RNASeq:CHROMOSOME_I_1 1 150000
# S_child Homol_data RNASeq_Hillier.Adult_spe-9:CHROMOSOME_I_2 100001 250000

my $blockname;
my $chromosome;
my ($indel, $change);
my ($newstart, $newend);
my ($remapped_start, $remapped_end);
my $line;
my ($start, $end);

# pass 1
# read the ace file to get the remapped virtual blocks
my %virtual_blocks;
open (IN, "< $input") || die "can't open input file $input\n";
while ($line = <IN>) {
  chomp $line;
  my $tag;
  if ($line =~ /Homol_data\s+(\S+)\s+(\d+)\s+(\d+)/) {
    ($blockname, $start, $end) = ($1, $2, $3);
    $tag = "Homol_data";

    if (! exists $virtual_blocks{$blockname}) {
      ($chromosome) = ($blockname =~ /_([IVX]+)_\d+$/);
      ($newstart, $newend, $indel, $change) = $assembly_mapper->remap_ace($chromosome, $start, $end);
      push @{$virtual_blocks{$blockname}}, ($tag, $start, $end, $chromosome, $newstart, $newend);
    }
  } elsif ($line =~ /Feature_data\s+(\S+)\s+(\d+)\s+(\d+)/) {
    ($blockname, $start, $end) = ($1, $2, $3);
    $tag = "Feature_data";

    if (! exists $virtual_blocks{$blockname}) {
      ($chromosome) = ($blockname =~ /_([IVX]+)_\d+$/);
      ($newstart, $newend, $indel, $change) = $assembly_mapper->remap_ace($chromosome, $start, $end);
      push @{$virtual_blocks{$blockname}}, ($tag, $start, $end, $chromosome, $newstart, $newend);
    }
  }

}
close (IN);


# pass 2

# get the details, look up the chromosomal start position in
# %virtual_blocks, remap the position then subtract newstart position

open (IN, "< $input") || die "can't open input file $input\n";
open (OUT, "> $output") || die "can't open output file $output\n";
while ($line = <IN>) {
  chomp $line;

  # Homol_data : "RNASeq_Hillier.Adult_spe-9:CHROMOSOME_III_1"
  if ($line =~ /Homol_data\s+:\s+(\S+)/) { 
    $blockname = $1;
    $blockname =~ s/"//g; #"
    ($chromosome) = ($blockname =~ /_(\S+)_\d+$/);
    print OUT "$line\n";

  # DNA_homol "Adult_spe-9_23dC_8days_post-L4_molt_bundle_of_reads_supporting_intron_III_10073_10364_+_wb170" RNASeq_Hillier.Adult_spe-9 1.0 10047 10072 1 26
  } elsif ($line =~ /DNA_homol\s+\S+\s+\S+\s+\S+\s+\d+\s+\d+\s+\d+\s+\d+/) {
    my ($homol, $id, $type, $score, $hstart, $hend, $cb_start, $cb_end) = split(/\s+/, $line);
    $id =~ s/"//g; #"
    if (! defined $blockname) {die "blockname not defined\n";}
    if (! exists $virtual_blocks{$blockname}) {die "can't find virtual_blocks of blockname=$blockname\n";}
    my ($tag, $vstart, $vend, $vchromosome, $newstart, $newend) = @{$virtual_blocks{$blockname}};
    my $oldstart = $vstart + $hstart-1;
    my $oldend = $vstart + $hend-1;
    ($remapped_start, $remapped_end, $indel, $change) = $assembly_mapper->remap_ace($vchromosome, $oldstart, $oldend);
    $remapped_start -= $newstart-1;
    $remapped_end -= $newstart-1;
    # see if the length of this thing changed and if so change the cb_start/cb_end
    my $len1 = $remapped_end-$remapped_start+1; # remapped length
    my $len2 = $oldend-$oldstart+1; # original length
    my $diff = $len2 - $len1; # change in length
    if ($cb_start < $cb_end) {
      $cb_end += $diff;
    } else {
      $cb_start += $diff;
    }
    print OUT "$homol \"$id\" $type $score $remapped_start $remapped_end $cb_start $cb_end\n";

  # Sequence : "CHROMOSOME_I"
  } elsif ($line =~ /Sequence : "CHROMOSOME_(\S+)"/) {
    $chromosome = $1;
    print OUT "$line\n";

  # S_child Homol_data RNASeq_Hillier.Adult_spe-9:CHROMOSOME_II_1 1 150000
  } elsif ($line =~ /Homol_data\s+(\S+)\s+(\d+)\s+(\d+)/) {	
    ($blockname, $start, $end) = ($1, $2, $3);
    my ($tag, $vstart, $vend, $vchromosome, $newstart, $newend) = @{$virtual_blocks{$blockname}};
    print OUT "S_child Homol_data $blockname $newstart $newend\n";
 
  # S_child Feature_data Confirmed_intron_RNASeq:CHROMOSOME_II_1 1 150000  
  } elsif ($line =~ /Feature_data\s+(\S+)\s+(\d+)\s+(\d+)/) {	
    ($blockname, $start, $end) = ($1, $2, $3);
    my ($tag, $vstart, $vend, $vchromosome, $newstart, $newend) = @{$virtual_blocks{$blockname}};
    print OUT "S_child Feature_data $blockname $newstart $newend\n";

  # Confirmed_intron 10017131 10017179 EST Adult_spe-9_23dC_8days_post-L4_molt_bundle_of_reads_supporting_intron_I_10017129_10017177_+_wb170
  } elsif ($line =~ /Confirmed_intron\s+(\d+)\s+(\d+)\s+(.+)/) {
    my $rest;
    ($start, $end, $rest) = ($1, $2, $3);
    ($remapped_start, $remapped_end, $indel, $change) = $assembly_mapper->remap_ace($chromosome, $start, $end);
    print OUT "Confirmed_intron $start $end $rest\n";


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

=head2 NAME - remap_chromosome_homol_data.pl

=head1 USAGE

=over 4

=item remap_chromosome_homol_data.pl  [-options]

=back

This script reads in an ace file of Homol_data mapping something to chromosomes and remaps the coordinates from the previous releases genome to the current genome.

script_template.pl MANDATORY arguments:

=over 4

=item -in file, input ace file from the previous release

=back

=over 4

=item -out file, output ace file for the current release

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
