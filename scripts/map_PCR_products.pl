#!/usr/local/bin/perl
#
# map_PCR_products
#
# Cronjob integrity check controls for generic ACEDB database.
#
# Usage: map_PCR_products [-options]
#
# 010927 :  dl : modified output print line to include speech marks. this prevents any acefile
#                parsing problems
# 010717 :  kj : modified to (hopefully) run faster and produce the final output files in one go
#
# 051114 :  mh6: major rewrite to get ready for mapping module
#
#############################################################################################

#############
# variables #
#############

use strict;
use warnings;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Class::Struct;
use Modules::Map_Helper;

##########################
# Script variables (run) #
##########################

my $maintainers = "All";
my %output      = ();
my $log         = Log_files->make_build_log();

########################
# command-line options #
########################

my $help;       # Help perdoc
my $test;       # Test mode
my $debug;      # Debug mode, output only goes to one user
my $verbose;    # verbose mode, more command line outout
my $acefile;    # output ace file

GetOptions(
    "debug=s"   => \$debug,
    "verbose"   => \$verbose,
    "test"      => \$test,
    "help"      => \$help,
    "acefile=s" => \$acefile
);

# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if ($debug) {
    print "DEBUG = \"$debug\"\n\n";
    ( $maintainers = $debug . '\@sanger.ac.uk' );
}

#############
# Paths etc #
#############

my $tace   = &tace;                                       # tace executable path
my $dbdir  = '/wormsrv2/autoace';                         # Database path
my $gffdir = $test ? "$dbdir/GFF_SPLITS/WS150" : "$dbdir/GFF_SPLITS/GFF_SPLITS";   # GFF splits directory
my @chromosomes = $test ? qw ( I ) : qw( I II III IV V X );                  # chromosomes to parse
my %genetype;                                             # gene type hash
my $outace = $acefile ? $acefile : "$dbdir/acefiles/PCR_mappings.ace";

################
# Structs      #
###############
struct( Exon => [ start => '$', stop => '$', type => '$', id => '$' ] );
struct( Gene => [ start => '$', stop => '$', exons => '@' ] );

###########################################
# get exons and PCRs out of the gff files #
###########################################

foreach my $chromosome (@chromosomes) {
    $log->write_to("Processing chromosome $chromosome\n");
    print "\nProcessing chromosome $chromosome\n" if ($verbose);
    my %genes;
    my %pcr;

    # Get PCR_product info from split GFF file
    open( GFF_in, "<$gffdir/CHROMOSOME_${chromosome}.PCR_products.gff" )
      || die "Failed to open PCR_product gff file\n\n";
    while (<GFF_in>) {
        chomp;
        s/\#.*// if ! (/\".+\#.+\"/);
        next unless /\S/;
        my @f = split /\t/;

        my ($name) = ( $f[8] =~ /PCR_product \"(.*)\"$/ );
        unless ($name) {
            $log->write_to("ERROR: Cant get name from $f[8]\n");
            next;
        }
        $pcr{$name} = [ $f[3], $f[4] ];
        print "PCR_product : '$name'\n" if ($verbose);
    }
    close(GFF_in);

    #################
    # read GFFs     #
    #################

    # Get exon info from split exon GFF files
    Map_Helper::get_from_gff( "$gffdir/CHROMOSOME_${chromosome}.exon.gff",
        'CDS', qw{\S}, \%genes );
    # Get exon info from split pseudogene exon GFF files
    Map_Helper::get_from_gff( "$gffdir/CHROMOSOME_${chromosome}.exon_pseudogene.gff",
        'Pseudogene', qw{\S}, \%genes );
    # Get exon info from split transcript exon GFF file
    Map_Helper::get_from_gff( "$gffdir/CHROMOSOME_${chromosome}.exon_noncoding.gff",
        'Transcript', qw{\S}, \%genes );
    # Get exon info from split UTR GFF files
    Map_Helper::get_from_gff( "$gffdir/CHROMOSOME_${chromosome}.UTR.gff",
        'Transcript', 'UTR', \%genes );

    print "Finished GFF loop\n" if ($verbose);

    ###################
    # make indexlists #
    ###################
    print "sorting genes\n" if ($debug);
    my @sorted_genes = sort {
             $genes{$a}->start <=> $genes{$b}->start
          || $genes{$a}->stop <=> $genes{$b}->stop
    } keys %genes;

    ##########
    # map it #
    ##########
    print "Find overlaps for PCR_product\n" if ($debug);
    # sub map_it(%output,%pcr,@sorted_genes,%genes)
   Map_Helper::map_it(\%output,\%pcr,\@sorted_genes,\%genes);
}

########################
# produce output files #
########################

open( OUTACE, ">$outace" )
  || die "Couldn't write to PCR_mappings.ace\n";

foreach my $mapped ( sort keys %output ) {
    print "mapped $mapped\tto ", join ' ',
      ( map { $_->id } @{ $output{$mapped} } ), "\n"
      if ($verbose);

    foreach my $exon ( @{ $output{$mapped} } ) {
        if ( $exon->type eq "CDS" ) {
            print OUTACE "CDS : \"",$exon->id,"\"\n";
            print OUTACE "Corresponding_PCR_product \"$mapped\"\n\n";
        }
        elsif ( $exon->type eq "Pseudogene" ) {
            print OUTACE "Pseudogene : \"",$exon->id,"\"\n";
            print OUTACE "Corresponding_PCR_product \"$mapped\"\n\n";
        }
        elsif ( $exon->type eq "Transcript" ) {
            print OUTACE "Transcript : \"",$exon->id,"\"\n";
            print OUTACE "Corresponding_PCR_product \"$mapped\"\n\n";
        }
    }
}
close(OUTACE);

##############################
# read acefiles into autoace #
##############################

unless ($test) {

    my $command =
      "pparse /wormsrv2/autoace/acefiles/PCR_mappings.ace\nsave\nquit\n";

    open( TACE, "| $tace -tsuser map_PCR_products $dbdir" )
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

__END__

=pod

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

=back

=cut
