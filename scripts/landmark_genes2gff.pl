#!/usr/local/bin/perl5.8.0 -w
#
# landmark_gene2gff.pl
#
# by Keith Bradnam
#
# script for creating extra GFF lines to indicate those genes that are landmark genes
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2006-03-03 09:52:51 $
use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

my $database;     # choose another database (defaults to autoace)
my $debug;        # debug mode, output only emailed to one person
my $help;         # help mode, show perldoc
my $test;         # use test mode in ~wormpub
my %landmarks;    # hash containing gene IDs as keys and public_name field as value
my $store;                                       # to specify storable file
my $verbose;

GetOptions(
	   "help"       => \$help,
	   "debug=s"    => \$debug,
	   "test"       => \$test,
	   "verbose"    => \$verbose,
	   "database=s" => \$database,
	   'store=s'    => \$store
	   );

############################
# recreate configuration   #
############################
my $wormbase;
if ($store) { $wormbase = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wormbase = Wormbase->new( -debug => $debug, -test => $test, ) }

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


###########################################
# Variables Part II (depending on $wormbase)    #
###########################################

# database/file paths and locations
my $basedir = $wormbase->basedir;

# set default database if -database not used
$database = $wormbase->autoace if !$database;

###############################
#
#  main part of script
#
###############################

# get list of landmark genes from database, and add to hash
&get_landmark_genes;

# now loop through GFF files to look for existing gene spans
my @chromosomes = qw( I II III IV V X );

foreach (@chromosomes) {

    $log->write_to("Processing chromosome $_\n");

    # open input/output streams
    open( OUT, ">".$wormbase->gff_splits."/CHROMOSOME_${_}_landmarks.gff" ) || die "Cannot open output file\n";
    open( GFF, "<".$wormbase->chromosomes."/CHROMOSOME_${_}.gff" ) || die "Can't read CHROMOSOME_${_}.gff file\n";

    while (<GFF>) {

        # only want to match the following lines
        # CHROMOSOME_II   gene    gene    23347   24428   .       +       .       Gene "WBGene00005017"
        next unless ( /gene/ && /WBGene/ );

        # more exact check by splitting line and checking fields
        my @data = split(/\t/);
        next unless ( $data[1] eq "gene" && $data[2] eq "gene" );

        # modify 9th GFF column to look up in hash to get CGC name (or Public name)
        my $gene = $data[8];
        $gene =~ s/Gene //;
        $gene =~ s/\"//g;
        chomp($gene);

        # check gene from GFF file with genes in hash to see if it is a landmark gene, write to output if so
        if ( $landmarks{$gene} ) {
            print OUT
              "$data[0]\tlandmark\tgene\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t$data[7]\tLocus $landmarks{$gene}\n";
        }
    }
    close GFF;
    close(OUT) || die "Couldn't close output file\n";
}

# tidy up and exit
$log->mail();
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

##############################
# get list of landmark genes
##############################

sub get_landmark_genes {

    $log->write_to("Getting list of landmark genes from $database\n");

    my $tace = $wormbase->tace;
    my $db   = Ace->connect(
        -path    => $database,
        -program => $tace
      )
      or $log->log_and_die("Connection failure: ". Ace->error);

    # only want landmark genes which will be in GFF files, i.e. those with Sequence name
    # this is virtually all of them
    my @landmarks = $db->fetch( -query => "Find Gene WHERE Landmark_gene AND Sequence_name" );

    # build hash
    foreach my $gene (@landmarks) {
        my $public_name = $gene->Public_name->name;
        $landmarks{ $gene->name } = $public_name;
    }
    $db->close;

}

__END__
                                                                                                   
=pod
                                                                                                   
=head2 NAME - landmark_genes2gff.pl
                                                                                                   
=head1 USAGE
                                                                                                   
=over 4
                                                                                                   
=item landmark_genes2gff.pl  [-debug -help -database=s -test]
                                                                                                   
=back

This script takes the list of 100 or so 'landmark' genes (those gene objects with a Landmark_gene tag)
and if that gene appears in the GFF files, then the script will write extra GFF output, i.e.

From this:

CHROMOSOME_II   gene    gene    23347   24428   .       +       .       Gene "WBGene00005017"

To this:

CHROMOSOME_II framework gene    23347   24428   .       -       .       Locus bli-2


The final field of the GFF file will contain the 'Public_name' field for the gene object.  This
should always be a CGC style name.

Extra GFF output files will be written to $database/GFF_SPLITS/WSXXX and then need to be appended
to the main GFF files at the end of the build.  This process is for CSHL who need this information.

                                                                                                   
landmark_genes2gff.pl MANDATORY arguments:
                                                                                                   
=over 4

=item none
                                                                                                   
=back

map_Alleles.pl  OPTIONAL arguments:
                                                                                                   
=over 4
                                                                                                   
=item -help
                                                                                                   
this stuff
                                                                                                   
=back
                                                                                                   
=over 4
                                                                                                   
=item -debug <user>

log output goes to <user>
                                                                                                   
=back
                                                                                                   
=over 4
                                                                                                   
=item -database
                                                                                                   
specify a database path (defaults to /wormsrv2/autoace)
                                                                                                   
=back

=over 4
                                                                                                   
=item -test
                                                                                                   
use test build environment in ~wormpub
                                                                                                   
=back
                                                                                                   
=head1 AUTHOR
 
=over 4
 
=item Keith Bradnam (krb@sanger.ac.uk)
 
=back
