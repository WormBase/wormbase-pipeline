#!/software/bin/perl -w
#
# Small script to parse out onlt the high quality RNASeq introns and store them in a new file
# with a new method.
#
# Example: perl RNASeq_intron_parser.pl -input /nfs/wormpub/BUILD_DATA/MISC_DYNAMIC/RNASeq_splice_elegans.ace -output /nfs/wormpub/BUILD_DATA/MISC_DYNAMIC/RNASeq_splice_elegans_high_qual.ace
#

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Coords_converter;
use Modules::Remap_Sequence_Change;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($input, $output, $species);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "input:s"    => \$input,
	    "output:s"   => \$output,
	    "species:s"  => \$species,
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

#################################
# check input arguments
#################################

$species = $wormbase->full_name;

#################################
my $database = $wormbase->autoace;
my @data;
my @feature_data;
# suck the data in
open (IN, "<$input") || die "Can't open $input\n";
open (OUT, ">$output") || die "Can't open $output\n";

# Print the method to the file

print OUT "\nMethod : \"RNASeq_splice_high_qual\"\n";
print OUT "Remark   \"RNASeq splice region. These are the regions spanned by the SRA RNASeq short reads aligned to the genome by tophat\/bowtie. This is the set of high quality introns receiving 1000+ reads\"\n";
print OUT "Remark   \"RNASeq splice region. These are the regions spanned by the SRA RNASeq short reads aligned to the genome by tophat\/bowtie. The value stored of the splice is the sum of the number of reads observed spanning this region from all aligned RNASeq libraries.\"\n";
print OUT "Colour   DARKRED\n";
print OUT "Show_up_strand\n";
print OUT "Score_by_width\n";
print OUT "Score_bounds     1.000000 100.000000\n";
print OUT "Bumpable\n";
print OUT "Right_priority   2.725\n";
print OUT "GFF_source       \"RNASeq_splice_high_qual\"\n";
print OUT "GFF_feature      \"intron\"\n";
print OUT "GFF_SO   \"SO:0000188\"\n";
print OUT "\n\n";

# Process the .ace input.
my $line;
my $flag = 0;
while ($line = <IN>) {
  unless ($line =~ /^\S+/) {next;}

  #Sequence : "2L52" - print OUT these but next lines are massively duplicated in the original file so only print 1
  if ($line =~ /^Sequence : (\S+)/) {
    print OUT "\n$line";
    $flag = 1;
    next;
  }

  #S_Child Feature_data 2L52:Confirmed_intron_RNASeq 1 4592
  if ($line =~ /^S_Child Feature_data (\S+):/) {
    my $sequence = $1;
    if ($flag =~ 1) {
      @feature_data = split" ",$line;
      print OUT "S_Child Feature_data $sequence:Confirmed_intron_RNASeq_high_qual $feature_data[3] $feature_data[4]\n";
      $flag = 0;
    }
    else {
      next;
    }
  }


  #Feature_data : "2L52:Confirmed_intron_RNASeq" - Print these
  if ($line =~ /^Feature_data \: \"(\S+)\:/) {
    @data=split/\s+/,$line;
    print OUT "\nFeature_data : \"$1:Confirmed_intron_RNASeq_high_qual\"\n";
    next;
  }

  #Feature RNASeq_splice 3292 3344 6 "SRX208780 1" - filter these
  if($line =~ /^Feature RNASeq_splice/){
    @data=split/\s+/,$line;
    if ($data[4] > 1000) {
      my @slice = @data[5 .. $#data];
      print OUT "Feature RNASeq_splice_high_qual $data[2] $data[3] $data[4] @slice\n";
    }
    next;
  }
}

close(OUT);
close(IN);

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

