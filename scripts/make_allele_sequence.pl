#!/nfs/team71/worm/mh6/bin/perl
# map_Alleles.pl
#
# by Anthony Rogers
#
# This maps alleles to the genome based on their flanking sequences
#
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2006-07-17 12:31:52 $
# SubVersion :     $Revision: 1.1 $

use strict;
use lib $ENV{'CVS_DIR'};
use Feature_mapper;
use Wormbase;
use Ace;
use Getopt::Long;
use Data::Dumper;

my @opts = @ARGV;

#######################################
# command-line options                #
#######################################

my $debug;       # debug mode, output goes to /tmp
my $limit;       # limit number of alleles to map
my $database;    # specify database to map alleles in, default is autoace
my $release;     # specify release number, defaults to current release number
my $help;        # help mode, show perldoc
my $restart = 'go';    # specify an allele name from which to start mapping
my $no_parse;          # turn off loading of data to $database
my $list;              # read in alleles to map from file rather than from database
my $verbose;           # verbose mode, extra output to screen
my $store;             # defines stored configuration~/src/perl/mapper/map_Alleles.pl
my $test;
my $outdir;            # acefile output directory

GetOptions(
    'debug=s'    => \$debug,
    'limit=i'    => \$limit,
    'database=s' => \$database,
    'release=i'  => \$release,
    'help'       => \$help,
    'restart=s'  => \$restart,
    'no_parse'   => \$no_parse,
    'list=s'     => \$list,
    'verbose'    => \$verbose,
    'store=s'    => \$store,
    'test'       => \$test,
    'outdir=s'   => \$outdir,
);

# check command line options
if ($help) { print `perldoc $0`; exit }
die "No need for -release option in debug mode" if ( $release && $debug );

###########################
# recreate configuration   #
############################
my $wb;

if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug, -test => $test, ) }

###########################################
# Variables Part II (depending on $wb)    #
###########################################
$test  = $wb->test  if $wb->test;     # Test mode
$debug = $wb->debug if $wb->debug;    # Debug mode, output only goes to one user

my $maintainers = "All";
my $tace        = $wb->tace;
my $error_count = 0;                  # for tracking alleles that don't map

$release = $wb->get_wormbase_version unless ( defined $release );
my $ace_file = ( $outdir ? $outdir : $wb->acefiles ) . "/allele_seque_mapping.WS$release.ace";

$database = $wb->autoace unless $database;    # ACeDB

# set variables depending on whether in debug (test) mode or not
if ($debug) {
    $maintainers = "$debug\@sanger.ac.uk";
    print "DEBUG = \"$debug\"\n\n";
}

# check if it needs to submit stuff into lsf

################################################
# Set hashes
################################################
my %to_map;         # will be set if -list option is specified, stores list of alleles to map
my %allele_data;    # allele data from mapping data, hash of arrays:
my $log = Log_files->make_build_log($wb);    # create log file, open output file handles

###############################################################################################
#
#                                   MAIN part of script
#
###############################################################################################

&process_list if ($list);                    # read in list of alleles to map if -list is specified
&map_alleles($database);                     # map alleles, the main routine
&load_alleles_to_database unless ($no_parse);    # load acefile to autoace (unless -no_parse specified)

#####################################
# tidy up, email log, and exit
#####################################
if ( $error_count > 0 ) {
    print "\nERROR: problems with $error_count alleles\n\n";
    $log->write_to("\nProblems with $error_count alleles\n");
    $log->mail( "$maintainers", "BUILD REPORT: map_Alleles.pl $error_count ERRORS!" );
}
else {
    $log->mail( "$maintainers", "BUILD REPORT: map_Alleles.pl" );
}
exit(0);

##################################################
#                                                #
#            S U B R O U T I N E S               #
#                                                #
##################################################

#########################################################
#             Main allele mapping loop                  #
#########################################################

sub map_alleles {
    my ($database) = @_;

    # First get allele info from database
    my $db = Ace->connect( -path => $database ) || do { print "$database Connection failure: ", Ace->error; die(); };
    my $mapper = Feature_mapper->new( $database, undef, $wb );
    my @alleles = $db->fetch( -query => 'Find Variation WHERE Flanking_sequences' );    #very slow :-(

    open( OUT, ">$ace_file" ) or die "cant open $ace_file\n";

  ALLELE:
    foreach my $allele (@alleles) {
        my $name = $allele->name;
        my ( $sequence, $clone, $left, $right );
        my $go          = 0;
        my $orientation = '.';

        print "\nMapping $name\n" if $verbose;

        # debug facility - this bit is so that it can restart from given allele name
        unless ( "$restart" eq "go" ) {
            if ( "$restart" eq "$name" ) { $restart = 'go' }
            else { print "skipping $name\n" if $verbose; next }
        }

        # debug facility - after the restart means you can specify where to start and how many to do
        if ($limit) { last if $error_count++ >= $limit }
        if ($list) { next unless defined $to_map{$name} }

        # check that allele is attached to a valid sequence object
        $sequence = $allele->Sequence;
        if ( !defined $sequence ) {
            $log->write_to("ERROR: $name has missing Sequence tag\n");
            $error_count++;
            next ALLELE;
        }
        elsif ( !defined $sequence->Source ) {
            $log->write_to("ERROR: $name connects to Sequence $sequence which has no Source tag\n");
            $error_count++;
            next ALLELE;
        }
        unless ( $allele->Flanking_sequences and $allele->Flanking_sequences->right ) {
            $log->write_to( "ERROR: ", $allele->name, " is missing a flanking sequence!\n" );
            next ALLELE;
        }

        #map position on genome
        #(0)type,
        #(1)5'flank_seq ,
        #(2)3'flank_seq
        #(3)CDS
        #(4)end of 5'match
        #(5)start of 3'match
        #(6)clone
        #(7)chromosome
        #(8)strains containing this allele

        # get allele mapping
        # grab both flanking sequences from $database
        $left  = lc $allele->Flanking_sequences->name;
        $right = lc $allele->Flanking_sequences->right->name;

        # warn if flanking sequence is missing
        unless ( $left and $right ) {
            $log->write_to("ERROR: $name does not have two flanking sequences\n");
            $error_count++;
            next ALLELE;
        }
        $allele_data{$name}[0] = $allele->Method;
        $allele_data{$name}[1] = $left;
        $allele_data{$name}[2] = $right;

        # map allele using Feature_mapper.pm, store results of map in @map, warn if mapping failed
        my @map = $mapper->map_feature( $sequence->name, $left, $right );

        if ( "$map[0]" eq "0" ) {
            $log->write_to("ERROR: Couldn't map $name with $left and $right to seq $sequence\n");
            $error_count++;
            next ALLELE;
        }

        unless ( ( $allele->Type_of_mutation and $allele->Type_of_mutation eq "Insertion" )
            or ( $allele->Variation_type and $allele->Variation_type eq "Transposon_insertion" ) )
        {

            # get coords of allele
            if ( $map[2] > $map[1] ) {

                # maps to fwd strand
                $map[1]++;
                $map[2]--;
                $orientation = '+';
            }
            else {
                $map[1]--;
                $map[2]++;
                $orientation = '-';
            }
        }

        $allele_data{$name}[4] = $map[1];
        $allele_data{$name}[5] = $map[2];
        $allele_data{$name}[6] = $map[0];
        print OUT "\nSequence : \"$allele_data{$allele}[6]\"\nAllele $name $allele_data{$name}[4] $allele_data{$name}[5]\n\n";

        # retrieve list of overlapping CDSs/Transcripts

        # close database and file handles
    }
    $db->close;
    close OUT;
}
#####################################################################################################
####################################################################################
# small subroutine to load data into %to_map hash from file if -list is specified
####################################################################################

sub process_list {
    open( LIST, "<$list" ) || die "Cannot open file specified by $list\n";
    while (<LIST>) { chomp; $to_map{$_} = 1 }
    close(LIST);
}

################################
# read acefiles into autoace   #
#                              #
# use autoace_minder.pl -load  #
################################

sub load_alleles_to_database {

    $log->write_to("\nStart parsing $ace_file in to $database\n\n");
    my $command = "autoace_builder.pl -load $ace_file -tsuser map_Alleles.pl";
    my $status  = system($command);
    if ( ( $status >> 8 ) != 0 ) { die "ERROR: Loading $ace_file file failed \$\? = $status\n" }
    $log->write_to("\nFinished parsing $ace_file in to $database\n\n");
}

__END__

=pod

=head2 NAME - map_Alleles.pl

=head1 USAGE

=over 4

=item map_alleles.pl  [-debug -limit -database=s -release=s -verbose -restart=s -list=s -no_parse -noupdate -test]

=back

This script:

Gets alleles with flanking sequence from designated database and maps them to a clone or superlink using SCAN, then checks which if any CDSs they overlap with.

Also requires that the allele has as a "seed" sequence either a Sequence or CDS (a locus with an associated
genomic_sequence will also work).

Also writes two files for updating allele->Sequence and Allele->CDS in geneace, one to remove the current connections and one to enter the new ones.

Outputs acefiles which are loaded in to the same database.

map_Alleles.pl MANDATORY arguments:

=over 4

=item none

=back

map_Alleles.pl  OPTIONAL arguments:

=over 4

=item -help

this stuff

=back

=over 4

=item -debug 
   
output goes to [user] and uses current_DB

=back

=over 4

=item -noupdate

don't update/touch the mysql database as it might screw up the build

=back

=over 4

=item -limit 

limit the number of alleles mapped (debug tool)

=back

=over 4

=item -list 

list the alleles you want mapped (debug tool) - as filename ie -list "$mapping_dir/to_map"

=back

=over 4

=item -database 

specify which database to read info from and load mapping results in to

=back

=over 4

=item -restart 

choose which allele to start with. all preceding (alphabetically) alleles will be skipped

=back

=over 4

=item -verbose 

greater indication of what procedured are being used to map the allele

=back

=over 4

=item -release

select a version of the database other than that being built, do not use in debug mode

=back




=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk (if run as part
of build.

=item

=item Must be run AFTER gff splitting has produced CHROMOSOME_*.genes.gff

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers ( ar2@sanger.ac.uk)

=back

=cut
