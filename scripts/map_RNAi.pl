#!/nfs/team71/worm/mh6/bin/perl
# map_RNAi.pl
# Add information to RNAi objects based on overlaps in GFF files
#
# by Kerstin Jekosch
#
# Version: $Version: $
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2006-03-31 12:56:33 $

use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Modules::GFF_sql;
use Wormbase;
use Getopt::Long;
use Ace;

###############
# variables   #
###############

my $maintainers = "All";
my %rnai2genes;    # for RNAi
my %rnai2exp;      # for Expression profiles

################################
# command-line options         #
################################

my $help;          # Help perldoc
my $test;          # Test mode
my $debug;         # Debug mode, verbose output to user running script
my $verbose;
my $load;          # use to automatically load file to autoace
my $ace;           # specify where to put the output ace file
my $store;         # specify a frozen configuration file
my $chrom;         # specify a chromosome

my $gffdir;        # use a specific GFF_SPLIT directory
my $dbdir;         # connect to a different acedb

GetOptions(
    "debug=s"      => \$debug,
    "verbose"      => \$verbose,
    "test"         => \$test,
    "help"         => \$help,
    "load"         => \$load,
    "acefile=s"    => \$ace,
    'store=s'      => \$store,
    'chromosome=s' => \$chrom,
    'gffdir=s'     => \$gffdir,
    'dbdir=s'      => \$dbdir
);

# Display help if required
&usage("Help") if ($help);

############################
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

# Use debug mode?
if ($debug) {
    print "DEBUG = \"$debug\"\n\n";
    ( $maintainers = $debug . '\@sanger.ac.uk' );
}

my $log = Log_files->make_build_log($wb);

##############
# Paths etc. #
##############

my $tace = $wb->tace;    # tace executable path
$dbdir  = $wb->autoace    if ( !$dbdir );     # Database path
$gffdir = $wb->gff_splits if ( !$gffdir );    # GFF_SPLITS directory
my @chromosomes = $test ? qw ( IV ) : qw( I II III IV V X );    # chromosomes
@chromosomes = ($chrom) if $chrom;
my $acefile = $ace ? $ace : $wb->acefiles . "/RNAi_mappings.ace";
################
# Structs      #
################
use Class::Struct;
struct( Exon => [ start => '$', stop => '$', type => '$', id => '$' ] );

##########################
# MAIN BODY OF SCRIPT
##########################

###########################################################
# get exons, RNAis and Expr_profiles out of the gff files #
###########################################################

my $map = GFF_sql->new( { -build => 1 } );    # connect to build database

foreach my $chromosome (@chromosomes) {

    print "Processing chromosome $chromosome\n" if $verbose;
    $log->write_to("Processing chromosome $chromosome\n");

    ################
    # GFF database part
    $map->clean("CHROMOSOME_$chromosome");    # reset the chromosome table

    foreach my $end ( 'Expr_profile', 'curated', 'Non_coding_transcript', 'Pseudogene' ) {
        my $file = "$gffdir/CHROMOSOME_${chromosome}_${end}.gff";
        $map->generate_tags($file);
        $map->load_gff( $file, "CHROMOSOME_$chromosome" );
    }

    my %genes;
    my %exon;
    my %expression;

    # loop through the split GFF RNAi file
    # New RNAi lines : CHROMOSOME_I    RNAi_primary    RNAi_reagent    1681680 1683527 .       .       .       Target "RNAi:WBRNAi00004820" 1 1848

    print "Loop through primary RNAi GFF file CHROMOSOME_${chromosome}\n" if ($verbose);
    open( GFF, "<$gffdir/CHROMOSOME_${chromosome}_RNAi_primary.gff" ) || die "Failed to open RNAi gff file CHROMOSOME_${chromosome}_RNAi_primary.gff :$!\n\n";
    while (<GFF>) {
        chomp;
        s/^\#.*//;
        next unless /RNAi_primary\s+RNAi_reagent/;
        my @line = split /\t/;
        my %rnai_tmp;

        my ($name) = ( $line[8] =~ /\"RNAi:(\S+.+)\"\s+\d+\s+\d+$/ );

        # NB store type of RNAi (primary) with hits
        $rnai_tmp{"${name}_p"} = [ $line[3], $line[4], "primary" ];
        my @hits = $map->get_chr( "CHROMOSOME_$chromosome", { start => $line[3], stop => $line[4] } );
        foreach my $hit (@hits) {
            my $exon = ${ &to_exon($hit) };
            if ( $exon->type eq 'Expr_profile' ) { push @{ $rnai2exp{"${name}_p"} }, $exon }
            elsif ( $exon->type ) { push @{ $rnai2genes{"${name}_p"} }, $exon }
        }

    }
    close(GFF);

    # note which is secondary by adding "_s" to the RNAi
    print "Loop through secondary RNAi GFF file CHROMOSOME_${chromosome}\n" if ($verbose);
    open( GFF, "<$gffdir/CHROMOSOME_${chromosome}_RNAi_secondary.gff" )
      || die "Failed to open RNAi gff file CHROMOSOME_${chromosome}_RNAi_secondary.gff:$!\n\n";
    while (<GFF>) {
        chomp;
        s/^\#.*//;
        next unless /RNAi_secondary\sRNAi_reagent/;
        my @line = split /\t/;
        my %rnai_tmp;

        my ($name) = ( $line[8] =~ /\"RNAi:(\S+.+)\"\s+\d+\s+\d+$/ );

        my @hits = $map->get_chr( "CHROMOSOME_$chromosome", { start => $line[3], stop => $line[4] } );
        foreach my $hit (@hits) {
            my $exon = ${ &to_exon($hit) };
            if ( $exon->type eq 'Expr_profile' ) { push @{ $rnai2exp{"${name}_s"} }, $exon }
            elsif ( $exon->type ) { push @{ $rnai2genes{"${name}_s"} }, $exon }
        }

    }
    close(GFF);
}

# store some statistics
# count number of primary RNAis (after removing duplicate mapping to same gene)
my $no_of_rnais_primary = 0;

# count number of secondary RNAis (after removing duplicate mappings to same gene)
my $no_of_rnais_secondary = 0;

############################
# sort the output for RNAi #
############################

print "Remove duplicates for RNAi->Gene\n" if ($verbose);

foreach my $RNAiID ( keys %rnai2genes ) {
    if ( $RNAiID =~ /(.*)_p$/ ) {
        my %genes2rnaip;
        my %genes2rnais;

        # basic idea is assign a id+type to a hash to use for reassigning the rnai2gene connections

        # remove duplicate primary connections
        @{ $rnai2genes{$RNAiID} } =
          grep { ( ( !$genes2rnaip{ $_->id . $_->type } ) && ( $genes2rnaip{ $_->id . $_->type } = $RNAiID ) ) } @{ $rnai2genes{$RNAiID} };
        next if ( !defined ${ $rnai2genes{"$1_s"} }[0] );

        # for secondary connections remove duplicate secondaries as well as a secondary to an existing primary
        @{ $rnai2genes{"$1_s"} } =
          grep { ( ( !$genes2rnais{ $_->id . $_->type } ) && ( $genes2rnais{ $_->id . $_->type } = $RNAiID ) ) } @{ $rnai2genes{"$1_s"} };

        foreach my $index ( 0 .. scalar( @{ $rnai2genes{"$1_s"} } ) - 1 ) {
            delete ${ $rnai2genes{"$1_s"} }[$index]
              if $genes2rnaip{ ${ $rnai2genes{"$1_s"} }[$index]->id . ${ $rnai2genes{"$1_s"} }[$index]->type };    # eek ... slightly iffy
        }
    }
}

########################
# produce output files #
########################

print "Produce output file\n" if ($verbose);

open( DELACE, ">${acefile}_del" );

# Produce connections for RNAi->Genes
# remove existing connections
foreach my $mapped ( keys %rnai2genes ) {
    $mapped =~ s/_[ps]$//;
    print DELACE "RNAi : $mapped\n";
    print DELACE "-D Inhibits\n\n";
}
close DELACE;

open( OUTACE, ">$acefile" ) || die "Failed to open $acefile file\n";

# open autoace connection with aceperl
my $db = Ace->connect( -path => $dbdir, -program => $tace ) || die "Couldn't connect to $dbdir\n", Ace->error;
my %inverse;
foreach my $mapped ( keys %rnai2genes ) {

    $mapped =~ /(.*)_([ps])$/;
    my $id = $1;
    my $type = $2 eq 'p' ? 'primary' : 'secondary';    # $type holds whether the RNAi mapping to this gene is primary or secondary
    my $seq;
    my $worm_gene;    # CDS, Transcript, or Pseudogene name

    foreach my $exon ( @{ $rnai2genes{$mapped} } ) {
        next if !( defined $exon );

        print OUTACE "// $type RNAi '$mapped' is mapped to ", $exon->id, " (", $exon->type, ")\n" if ($verbose);
        next if !( $exon->type eq 'CDS' || $exon->type eq 'Pseudogene' || $exon->type eq 'Transcript' );

        if ( $type eq "primary" ) {
            $no_of_rnais_primary++;
        }
        else {
            $no_of_rnais_secondary++;
        }

        print OUTACE "RNAi : \"$id\"\n";

        # Does this CDS have a Gene object?
        if ( $exon->type eq "CDS" ) {
            $seq = $db->fetch( -class => 'CDS', -name => $exon->id );
            if ( defined $seq ) {
                print OUTACE "Predicted_gene \"", $exon->id, "\" Inferred_automatically \"RNAi_$type\"\n";
                if ( defined( $seq->at('Visible.Gene') ) ) {
                    my ($gene) = ( $seq->get('Gene') );
                    print OUTACE "Gene \"$gene\" Inferred_automatically \"RNAi_$type\"\n";
                    push @{ $inverse{$gene} }, [ $id, $type ];
                }
            }
            else {
                print '*** WARNING - skipping missing gene ', $exon->id, "\n";
            }
        }

        # Does this Pseudogene have a Gene object?
        elsif ( $exon->type eq "Pseudogene" ) {
            $seq = $db->fetch( -class => 'Pseudogene', -name => $exon->id );
            if ( defined $seq ) {
                print OUTACE "Pseudogene \"", $exon->id, "\" Inferred_automatically \"RNAi_$type\"\n";
                if ( defined( $seq->at('Visible.Gene') ) ) {
                    my ($gene) = ( $seq->get('Gene') );
                    print OUTACE "Gene \"$gene\" Inferred_automatically \"RNAi_$type\"\n";
                    push @{ $inverse{$gene} }, [ $id, $type ];
                }
            }
            else {
                print '*** WARNING - skipping missing gene ', $exon->id, "\n";
            }
        }

        # Does this transcript have a Gene object?
        elsif ( $exon->type eq "Transcript" ) {
            $seq = $db->fetch( -class => 'Transcript', -name => $exon->id );
            if ( defined $seq ) {
                print OUTACE "Transcript \"", $exon->id, "\" Inferred_automatically \"RNAi_$type\"\n";
                if ( defined( $seq->at('Visible.Gene') ) ) {
                    my ($gene) = ( $seq->get('Gene') );
                    print OUTACE "Gene \"$gene\" Inferred_automatically \"RNAi_$type\"\n";
                    push @{ $inverse{$gene} }, [ $id, $type ];
                }
            }
            else {
                print "*** WARNING - skipping missing gene ", $exon->id, "\n";
            }
        }
        print OUTACE "\n";
    }
}

$db->close;

print OUTACE "\n\n//Expression profiles\n";

# Produce connections for RNAi->Expr_profile
my %printedexp;
foreach my $mapped ( keys %rnai2exp ) {
    $mapped =~ /(.*)_[ps]$/;
    my $name = $1;
    foreach my $exon ( @{ $rnai2exp{$mapped} } ) {
        if ( !$printedexp{ $exon->id } ) {
            print OUTACE "\nRNAi : \"$name\"\n";
            print OUTACE "Expr_profile ", $exon->id, "\n";
            $printedexp{ $name . $exon->id } = 1;
        }
    }
}

# Write primary/secondary type evidence for RNAi tag in ?Gene model
# (This is done explicitly because you can't (easily) make a XREF between #Evidence tags in the models)
foreach my $gene ( keys %inverse ) {
    for ( my $n = 0 ; $n < ( scalar @{ $inverse{$gene} } ) ; $n++ ) {
        my $RNAi = $inverse{$gene}->[$n][0];
        my $type = $inverse{$gene}->[$n][1];
        print "Gene: $gene RNAi: $RNAi type: $type\n" if ($verbose);
        print OUTACE "\nGene : \"$gene\"\n";
        print OUTACE "Experimental_info RNAi_result  \"$RNAi\" Inferred_automatically \"RNAi_$type\"\n";
    }
}

close(OUTACE);

#########################################################
# read acefiles into autoace (unless running test mode) #
#########################################################
$wb->load_to_database( $dbdir, $acefile, "RNAi_mappings" ) if ($load);

####################################
# print some statistics to the log #
####################################

# count secondary RNAis which map to a gene which is already mapped by that RNAi as a primary
$log->write_to("\n\nStatistics\n");
$log->write_to("----------\n\n");
$log->write_to("No. of primary   RNAi links to genes written to database: $no_of_rnais_primary\n");
$log->write_to("No. of secondary RNAi links to genes written to database: $no_of_rnais_secondary\n\n");

$log->mail( "$maintainers", "BUILD REPORT: $0" );

##############################################################
# Subroutines
##############################################################

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
    elsif ( $hit->{feature} eq 'Expr_profile' ) { $type = 'Expr_profile' }
    else { $type = "feature:" . $hit->{feature} . " source:" . $hit->{source} }
    my $exon = Exon->new( id => get_id( $hit->{fluff} ), start => $hit->{start}, stop => $hit->{stop}, type => $type );

    return \$exon;
}

############################################

__END__

=pod

=head2 NAME - map_RNAi.pl

=head1 USAGE

=over 4

=item map_RNAi.pl [-options]

=back

map_RNAi.pl calculates the overlap between the genomic regions used in RNAi
experiments and the CDS, transcript and pseudogene coordinates in the WS
database release. It will generate an acefile which will remove any existing
connections and make new ones. It will check the current database and make Gene
connections where valid and attach expression_profiles as needed. This acefile
is then loaded into autoace if -load is specified.

map_RNAi mandatory arguments:


=over 4

=item none

=back

map_RNAi optional arguments:

=over 4

=item -debug, debug mode

=item -verbose, additional debugging stuff

=item -test, Test mode, generate the acefile but do not upload them, use only one chromosome 

=item -help, Help pages

=item -acefile, write to a specific acefile

=item -load, loads file to autoace

=item -store specifiy configuration file

=back

=cut
