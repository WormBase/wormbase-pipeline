#!/usr/local/bin/perl
#
# map_Oligo_set
#
# Usage: map_Oligo_set.pl [-options]
#
# 010927 :  dl : modified output print line to include speech marks. this prevents any acefile
#                parsing problems
# 010717 :  kj : modified to (hopefully) run faster and produce the final output files in one go
#
# 051114 :  mh6 : rewrote it for better readability
#
#############################################################################################

use lib $ENV{'CVS_DIR'};
use strict;
use Wormbase;
use Getopt::Long;
use warnings;
use Class::Struct;
use Modules::Map_Helper;

##########################
# Script variables (run) #
##########################

my $maintainers = "All";
my %output      = ();
my %finaloutput = ();

########################
# command-line options #
########################

my $help;       # Help perldoc
my $test;       # Test mode
my $debug;      # Debug mode, output only goes to one user
my $verbose;    # verbose mode, more command line outout
my $acename;    # specify a custom acefile
my $store;      # specify a frozen configuration file
my $no_parse;   # don't parse the acefile

GetOptions(
    'debug=s'   => \$debug,
    'verbose'   => \$verbose,
    'test'      => \$test,
    'help'      => \$help,
    'acefile=s' => \$acename,
    'store=s'   => \$store,
    'no_parse'  => \$no_parse
);

# Display help if required
&usage("Help") if ($help);

############################
# recreate configuration   #
# ##########################
my $wb;
if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug, -test => $test, ) }

###########################################
# Variables Part II (depending on $wb)    #
###########################################

$test  = $wb->test  if $wb->test;     # Test mode
$debug = $wb->debug if $wb->debug;    # Debug mode, output only goes to one user

# Use debug mode?
if ($debug) {
    print "DEBUG = \"$debug\"\n\n";
    ( $maintainers = $debug . '\@sanger.ac.uk' );
}

#############
# Paths etc #
#############
my $tace    = $wb->tace;                                                    # tace executable path
my $dbdir   = $wb->autoace;                                                 # Database path
my $gffdir  = $wb->gff_splits;                                              # GFF splits directory
my $acefile = $acename ? $acename : "$dbdir/acefiles/Oligo_mappings.ace";

my @chromosomes = $test
  ? qw( III )
  : qw( I II III IV V X );    # chromosomes to parse (TEST_BUILD should be III)

################
# Structs      #
################
use Class::Struct;
struct( Exon => [ start => '$', stop => '$', type => '$', id => '$' ] );
struct( Gene => [ start => '$', stop => '$', orientation => '$',parent => '$', exons => '@' ] );

################
# Open logfile #
################
my $log = Log_files->make_build_log($wb);

#####################################################
# get exons and Oligo_set info out of the gff files #
#####################################################

foreach my $chromosome (@chromosomes) {
    $log->write_to("Processing chromosome $chromosome\n");
    print "\nProcessing chromosome $chromosome\n" if ($verbose);
    my %genes = ();
    my %Oligo = ();

    # Get Oligo set info from split GFF file
    open( GFF, "<$gffdir/CHROMOSOME_${chromosome}_Oligo_set.gff" )
      || die "Failed to open $gffdir/CHROMOSOME_${chromosome}_Oligo_set.gff : $!\n\n";
    while (<GFF>) {
        chomp;
        s/\#.*//;
        next unless /Oligo_set/;
#for WS156	next unless /exon/
        my @f = split /\t/;

        my ($name) = ( $f[8] =~ /Oligo_set \"(.*)\"$/ );
        unless ($name) {
            $log->write_to("ERROR: Cant get name from $f[8]\n");
            next;
        }
        $Oligo{$name} = [ $f[3], $f[4] ];
        print "Oligo_set : '$name'\n" if ($verbose);
    }
    close(GFF);
#######################
    # Get exon info from split UTR GFF files
    Map_Helper::get_from_gff( "$gffdir/CHROMOSOME_${chromosome}_UTR.gff", 'Transcript', 'UTR', \%genes );

    # Get exon info from split exon GFF files
    Map_Helper::get_from_gff( "$gffdir/CHROMOSOME_${chromosome}_curated.gff", 'CDS', qw{exon}, \%genes );

    # Get exon info from split pseudogene exon GFF files
    Map_Helper::get_from_gff( "$gffdir/CHROMOSOME_${chromosome}_Pseudogene.gff", 'Pseudogene', qw{exon}, \%genes );

    # Get exon info from split transcript exon GFF file
    Map_Helper::get_from_gff( "$gffdir/CHROMOSOME_${chromosome}_Non_coding_transcript.gff", 'Transcript', qw{exon}, \%genes );

    # Get exon info from split ncRNA exon GFF file
    Map_Helper::get_from_gff( "$gffdir/CHROMOSOME_${chromosome}_ncRNA.gff", 'Transcript', qw{exon}, \%genes );
#######################
    print "Finished GFF loop\n" if ($verbose);

    ###################
    # make indexlists #
    ###################

    print " sorting genes\n" if ($verbose);

    my @sorted_genes =
      sort { $genes{$a}->start <=> $genes{$b}->start || $genes{$a}->stop <=> $genes{$b}->stop } keys %genes;

    ##########
    # map it #
    ##########
    Map_Helper::map_it( \%output, \%Oligo, \@sorted_genes, \%genes );

}

########################
# produce output files #
########################
open( OUTACE, ">$acefile" ) || die "cannot open file $acefile\n";

foreach my $oligoname ( keys %output ) {
    print "mapped $oligoname\tto ", join ' ', ( map { $_->id } @{ $output{$oligoname} } ), "\n"
      if ($verbose);

    foreach my $exon ( @{ $output{$oligoname} } ) {
        if ( $exon->type eq "CDS" ) {
            print OUTACE "CDS : \"", $exon->id, "\"\n";
            print OUTACE "Corresponding_oligo_set \"$oligoname\"\n\n";
        }
        elsif ( $exon->type eq "Pseudogene" ) {
            print OUTACE "Pseudogene : \"", $exon->id, "\"\n";
            print OUTACE "Corresponding_oligo_set \"$oligoname\"\n\n";
        }
        elsif ( $exon->type eq "Transcript" ) {
            print OUTACE "Transcript : \"", $exon->id, "\"\n";
            print OUTACE "Corresponding_oligo_set \"$oligoname\"\n\n";
        }
    }
}
close(OUTACE);

##############################
# read acefiles into autoace #
##############################

unless ($test||$no_parse) {

    my $command = "pparse $acefile\nsave\nquit\n";

    open( TACE, "| $tace -tsuser map_Oligo_products $dbdir" )
      || die "Couldn't open tace connection to $dbdir\n";
    print TACE $command;
    close(TACE);
}

###############
# hasta luego #
###############

$log->mail( "$maintainers", "BUILD REPORT: $0" );

exit(0);

##############################################################
##### Subroutines #####
##############################################################

sub usage {
    my $error = shift;

    if ( $error eq "Help" ) {

        # Normal help menu
        system( 'perldoc', $0 );
        exit(0);
    }
}

################################
# hit to exon converter
# including: adding types based on source/features
sub to_exon {
    my $hit = shift;
    my $type;

    if    ( $hit->{feature} eq 'curated'               && $hit->{source} eq 'exon' ) { $type = 'CDS' }
    elsif ( $hit->{feature} eq 'Pseudogene'            && $hit->{source} eq 'exon' ) { $type = 'Pseudogene' }
    elsif ( $hit->{feature} eq 'Non_coding_transcript' && $hit->{source} eq 'exon' ) { $type = 'Transcript' }
    elsif ( $hit->{feature} eq 'Coding_transcript'     && $hit->{source} eq 'exon' ) { $type = 'Transcript' }
    
    else { $type = "feature:" . $hit->{feature} . " source:" . $hit->{source} }
    my $exon = Exon->new( id => get_id( $hit->{fluff} ), start => $hit->{start}, stop => $hit->{stop}, type => $type );

    return \$exon;
}

__END__

=pod

=head2 NAME - map_Oligo_products.pl

=head1 USAGE

=over 4

=item map_Oligo_products.pl [-options]

=back

map_Oligo_products.pl calculates the overlap between the mapped Oligo_products
and the CDS, transcript and pseudogene coordinates in the WS database release. 
It will generate an acefile for this data and upload into /wormsrv2/autoace

map_Oligo_products mandatory arguments:

=over 4

=item none

=back

map_Oligo_products optional arguments:

=over 4

=item -debug username, debug mode email goes to user specified by -debug

=item -verbose, toggles extra output

=item -test, Test mode, generate the acefile but do not upload 

=item -acefile, custom name for the acefile

=item -help, Help pages

=item -store, specify configuration file

=back

=cut
