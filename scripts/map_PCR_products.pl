#!/usr/local/bin/perl
# map_PCR_products
#
# Usage: map_PCR_products [-options]
#
# 010927 :  dl : modified output print line to include speech marks. this prevents any acefile
#                parsing problems
# 010717 :  kj : modified to (hopefully) run faster and produce the final output files in one go
#
# 051114 :  mh6: major rewrite to get ready for mapping module
# 060331 :  mh6: moved to new mapping + small cleanup of the surrounding code
#
#############################################################################################

use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Class::Struct;
use Modules::GFF_sql;

##########################
# Script variables (run)

my $maintainers = "All";
my %output      = ();

########################
# command-line options

my $help;       # Help perdoc
my $test;       # Test mode
my $debug;      # Debug mode, output only goes to one user
my $verbose;    # verbose mode, more command line outout
my $acefile;    # output ace file
my $store;      # specify a frozen configuration file
my $gffdir;     # specify some GFF_Split dir
my $no_parse;   # don't parse the acefile

GetOptions(
    "debug=s"   => \$debug,
    "verbose"   => \$verbose,
    "test"      => \$test,
    "help"      => \$help,
    "acefile=s" => \$acefile,
    'store=s'   => \$store,
    'gffdir=s'  => \$gffdir
);

# Display help if required
&usage("Help") if ($help);

############################
# recreate configuration

my $wb;
if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug, -test => $test, ) }

###########################################
# Variables Part II (depending on $wb)

$test  = $wb->test  if $wb->test;     # Test mode
$debug = $wb->debug if $wb->debug;    # Debug mode, output only goes to one user

# Use debug mode?
if ($debug) {
    print "DEBUG = \"$debug\"\n\n";
    $maintainers = $debug . '\@sanger.ac.uk';
}

#############
# Paths etc

my $tace  = $wb->tace;                # tace executable path
my $dbdir = $wb->autoace;             # Database path
$gffdir = $gffdir ? $gffdir : $wb->gff_splits;    # GFF splits directory
my @chromosomes = $test    ? qw ( I ) : qw( I II III IV V X );                 # chromosomes to parse
my $outace      = $acefile ? $acefile : $wb->acefiles . "/PCR_mappings.ace";

# Struct for historical reasons
struct( Exon => [ start => '$', stop => '$', type => '$', id => '$' ] );

# make a new log
my $log = Log_files->make_build_log($wb);

###########################################
# get exons and PCRs out of the gff files

my $map = GFF_sql->new( { -build => 1 } );                                     # connect to build database

foreach my $chromosome (@chromosomes) {
    $log->write_to("Processing chromosome $chromosome\n");
    print "\nProcessing chromosome $chromosome\n" if ($verbose);

    ################
    # GFF database part
    $map->clean("CHROMOSOME_$chromosome");                                     # reset the chromosome table

    foreach my $end ('curated', 'Non_coding_transcript', 'Pseudogene' ) {
        my $file = "$gffdir/CHROMOSOME_${chromosome}_${end}.gff";
        $map->generate_tags($file);
        $map->load_gff( $file, "CHROMOSOME_$chromosome" );
    }

    # Get PCR_product info from split GFF file
    &get_PCRs_from_GFF( $chromosome, 'Orfeome',   \%output, $map );
    &get_PCRs_from_GFF( $chromosome, 'GenePairs', \%output, $map );

}

######################
# %output magic to remove redundancy
print "Remove duplicates \n" if ($verbose);

foreach my $pcr ( keys %output ) {
        # remove duplicate primary connections
	my %genes2pcr;
        @{ $output{$pcr} } = grep { ( ( !$genes2pcr{ $_->id . $_->type } ) && ( $genes2pcr{ $_->id . $_->type } = $pcr ) ) } @{ $output{$pcr} };
}


########################
# produce output files

open( OUTACE, ">$outace" ) || die "Couldn't write to PCR_mappings.ace\n";

foreach my $mapped ( sort keys %output ) {
    print "mapped $mapped\tto ", join ' ', ( map { $_->id } @{ $output{$mapped} } ), "\n" if ($verbose);

    foreach my $exon ( @{ $output{$mapped} } ) {
        if ( $exon->type eq "CDS" ) {
            print OUTACE "CDS : \"", $exon->id, "\"\n";
            print OUTACE "Corresponding_PCR_product \"$mapped\"\n\n";
        }
        elsif ( $exon->type eq "Pseudogene" ) {
            print OUTACE "Pseudogene : \"", $exon->id, "\"\n";
            print OUTACE "Corresponding_PCR_product \"$mapped\"\n\n";
        }
        elsif ( $exon->type eq "Transcript" ) {
            print OUTACE "Transcript : \"", $exon->id, "\"\n";
            print OUTACE "Corresponding_PCR_product \"$mapped\"\n\n";
        }
    }
}
close(OUTACE);

##############################
# read acefiles into autoace

$wb->load_to_database( $wb->autoace, $outace, 'map_PCR_products', $log ) unless ($test||$no_parse);
$log->mail();

exit(0);

##############################################################
# Subroutines

sub usage {
    my $error = shift;
    if ( $error eq "Help" ) {

        # Normal help menu
        system( 'perldoc', $0 );
        exit(0);
    }
}

####################################
# returns second word without "
sub get_id {
    my $fluff = shift;
    $fluff =~ s/\"//g;
    my @fields = split " ", $fluff;
    return $fields[1];
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

############################
# unified function to iterate through the gff
# $log and $gffdir are not passed along :-(
sub get_PCRs_from_GFF {
    my ( $chromosome, $method, $pcr2gene, $map ) = @_;

    open GFF_in, "<$gffdir/CHROMOSOME_${chromosome}_$method.gff" or $log->log_and_die("Failed to open PCR_product gff file\n\n");

    while (<GFF_in>) {
        chomp;
        s/\#.*// if !(/\".+\#.+\"/);
        next unless /PCR_product/;
        my @f = split /\t/;

        my ($name) = ( $f[8] =~ /PCR_product \"(.*)\"$/ );
        unless ($name) {
            $log->write_to("WARNING: Cant get name from $f[8] , line $. in CHROMOSOME_${chromosome}_$method.gff\n");
            next;
        }

        my @hits = $map->get_chr( "CHROMOSOME_$chromosome", { start => $f[3], stop => $f[4] } );
        foreach my $hit (@hits) {
            my $exon = ${ &to_exon($hit) };
            push @{ $$pcr2gene{$name} }, $exon;
        }

        print "PCR_product : '$name'\n" if ($verbose);
    }
    close GFF_in;
}

__END__

=pod
use Hum::EMBL;

=head2 NAME - map_PCR_products.pl

=head1 USAGE

=over 4

=item map_PCR_products.pl [-options]

=back

map_PCR_products.pl calculates the overlap between the mapped PCR_products
and the CDS, transcript and pseudogene coordinates in the WS database release. 
It will generate an acefile for this data and upload into /wormsrv2/autoace

map_PCR_products mandatory arguments:

=over 4

=item none

=back

map_PCR_products optional arguments:

=over 4

=item -debug, Verbose/Debug mode

=item -test, Test mode, generate the acefile but do not upload 

=item -help, Help pages

=item -store specify configuration file

=back

=cut
