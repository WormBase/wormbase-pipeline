#!/software/bin/perl -w
#
# get_splice_flanks.pl                          
# 
# by ar2            
#
# Last updated by: $Author: ar2 $     
# Last updated on: $Date: 2010-04-28 16:00:45 $      
#
#
# This to get the sequences either side of exon splice sites.
# cat in a GFF file of the exons (or any other feature) and it will return some bases (default 20)
# either side of that feature eg
#
# grep curated $CURRENT/CHROMOSOMES/CHROMOSOME_III.gff | grep exon | perl get_splice_flanks.pl -store ~wormpub/BUILD/autoace/Elegans.store -database $CURRENT > ~/chrIII_splice_seqs.txt
#
#

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Sequence_extract;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($database, $flank);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "database:s" => \$database,
	    "flank:i"    => \$flank
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

$database = $wormbase->orgdb unless $database;
$flank = 20 unless $flank;

my $extract = Sequence_extract->invoke($database, 0, $wormbase);

while(<>) {
    my @data = split;
   #CHROMOSOME_III  curated exon    89068   89151   .       +       .       CDS "F54C4.1"
    my ($chrom, $exStart, $exEnd, $strand, $cds) = ($data[0],$data[3],$data[4],$data[6],$data[9]);
    $cds =~ s/\"//g;

    my $p5seq = $extract->Sub_sequence($chrom,$exStart-($flank+2),(2*$flank +2));
    my $p3seq = $extract->Sub_sequence($chrom,$exEnd-($flank+1),(2*$flank +2));

    next unless ($p5seq and $p3seq);
    if($strand eq '-'){
	$p5seq = $extract->DNA_revcomp($p3seq);
	$p3seq = $extract->DNA_revcomp($p5seq);
    }

    print "$cds:$exStart-$exEnd($strand)\t$p5seq\t$p3seq\n";
}

exit();






##############################################################
#
# Subroutines
#
##############################################################

