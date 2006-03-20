#!/nfs/team71/worm/mh6/bin/perl
# map_RNAi.pl
# Add information to RNAi objects based on overlaps in GFF files
#
# by Kerstin Jekosch
#
# Version: $Version: $
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2006-03-20 11:59:36 $

use strict;
use warnings;
use Modules::Map_Helper;
use lib $ENV{'CVS_DIR'};
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
my $exp_arg;       # only do expression profile mapping (for debugging)
my $store;         # specify a frozen configuration file
my $chrom;    # specify a chromosome

GetOptions(
    "debug=s"    => \$debug,
    "verbose"    => \$verbose,
    "test"       => \$test,
    "help"       => \$help,
    "load"       => \$load,
    "acefile=s"  => \$ace,
    'expression' => \$exp_arg,
    'store=s'    => \$store,
    'chromosome=s'=> \$chrom
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

my $tace        = $wb->tace;                                           # tace executable path
my $dbdir       = $wb->autoace;                                        # Database path
my $gffdir      = $wb->gff_splits;                                     # GFF_SPLITS directory
my @chromosomes = $test ? qw ( IV ) : qw( I II III IV V X );            # chromosomes
@chromosomes = ( $chrom ) if $chrom;
my $acefile     = $ace ? $ace : $wb->acefiles."/RNAi_mappings.ace";

################
# Structs      #
################
use Class::Struct;
struct( Exon => [ start => '$', stop => '$', type => '$', id => '$' ] );
struct( Gene => [ start => '$', stop => '$', exons => '@' ] );

##########################
# MAIN BODY OF SCRIPT
##########################

###########################################################
# get exons, RNAis and Expr_profiles out of the gff files #
###########################################################

foreach my $chromosome (@chromosomes) {

    print "Processing chromosome $chromosome\n" if $verbose;
    $log->write_to("Processing chromosome $chromosome\n");

    my %genes;
    my %exon;
    my %expression;


## here could go the new mapping part

    #map_it2

##
    # loop through the split GFF Expr_profile file (Expr_profile|experimental_result_region) -> Expr_profile
    print "Loop through Expr_profile GFF file CHROMOSOME_${chromosome}\n" if ($verbose);
    Map_Helper::get_from_gff( "$gffdir/CHROMOSOME_${chromosome}_Expr_profile.gff",
        'Expr_profile', qw{\S}, \%expression );
    print 'Loaded ', scalar( keys %expression ), " expression profiles\n" if $verbose;
    if ( !$exp_arg ) {

        # loop through the split GFF exon file (curated|exon) -> CDS
        print "Loop through exon GFF file CHROMOSOME_${chromosome}\n" if ($verbose);
        Map_Helper::get_from_gff( "$gffdir/CHROMOSOME_${chromosome}_curated.gff", 'CDS', qw{exon}, \%genes );

        # loop through the split GFF exon_pseudogene file (Pseudogene|exon) -> Pseudogene)
        print "Loop through exon_pseudogene GFF file CHROMOSOME_${chromosome}\n" if ($verbose);
        Map_Helper::get_from_gff( "$gffdir/CHROMOSOME_${chromosome}_Pseudogene.gff", 'Pseudogene', qw{exon}, \%genes );

        # loop through the split GFF exon_noncoding file  (Non_coding_transcript|exon) -> Transcript
        print "Loop through exon_noncoding GFF file CHROMOSOME_${chromosome}\n" if ($verbose);
        Map_Helper::get_from_gff( "$gffdir/CHROMOSOME_${chromosome}_Non_coding_transcript.gff", 'Transcript', qw{exon},
            \%genes );
    }
    print "Finished GFF loop\n" if ($verbose);

    ###################
    # make indexlists #
    ###################

    print "Index exons,genes,RNAi and Expr_profiles\n" if ($verbose);

    my @sorted_genes =
      sort { $genes{$a}->start <=> $genes{$b}->start || $genes{$a}->stop <=> $genes{$b}->stop } keys %genes;

    my @sorted_exp =
      sort { $expression{$a}->start <=> $expression{$b}->start || $expression{$a}->stop <=> $expression{$b}->stop }
      keys %expression;

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

        # NB store type of RNAi (primary) with the start/end positions
        $rnai_tmp{"${name}_p"} = [ $line[3], $line[4], "primary" ];

    	# map RNAis to genes
        Map_Helper::map_it( \%rnai2genes, \%rnai_tmp, \@sorted_genes, \%genes );
        # map RNAis to Expression profiles
        Map_Helper::map_it( \%rnai2exp, \%rnai_tmp, \@sorted_exp, \%expression );

    }
    close(GFF);

    # add the seondary RNAi hits to the same data structure
    # note which is secondary by adding "secondary" to the gene mapped to
    print "Loop through secondary RNAi GFF file CHROMOSOME_${chromosome}\n" if ($verbose);
    open( GFF, "<$gffdir/CHROMOSOME_${chromosome}_RNAi_secondary.gff" ) || die "Failed to open RNAi gff file CHROMOSOME_${chromosome}_RNAi_secondary.gff:$!\n\n";
    while (<GFF>) {
        chomp;
        s/^\#.*//;
        next unless /RNAi_secondary\sRNAi_reagent/;
        my @line = split /\t/;
	my %rnai_tmp;

        my ($name) = ( $line[8] =~ /\"RNAi:(\S+.+)\"\s+\d+\s+\d+$/ );

        # NB store type of RNAi (secondary) with the start/end positions
        $rnai_tmp{"${name}_s"} = [ $line[3], $line[4], "secondary" ];
    	# map RNAis to genes
        Map_Helper::map_it( \%rnai2genes, \%rnai_tmp, \@sorted_genes, \%genes );
        # map RNAis to Expression profiles
        Map_Helper::map_it( \%rnai2exp, \%rnai_tmp, \@sorted_exp, \%expression );

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
#        map { $genes2rnaip{ $_->id . $_->type } = $RNAiID } @{ $rnai2genes{$RNAiID} };
	@{$rnai2genes{$RNAiID}}= grep { ((! $genes2rnaip{ $_->id . $_->type }) && ($genes2rnaip{ $_->id.$_->type} = $RNAiID)) } @{ $rnai2genes{$RNAiID} };
        next if (! defined ${$rnai2genes{"$1_s"}}[0]);
        foreach my $index ( 0 .. scalar( @{$rnai2genes{"$1_s"}} ) - 1 ) {
            delete ${ $rnai2genes{"$1_s"} }[$index]
              if $genes2rnaip{ ${ $rnai2genes{"$1_s"} }[$index]->id . ${ $rnai2genes{"$1_s"} }[$index]->type };    #eek
        }
    }
}

########################
# produce output files #
########################

print "Produce output file\n" if ($verbose);

open( OUTACE, ">$acefile" ) || die "Failed to open RNAi_mappings.ace file\n";

# Produce connections for RNAi->Genes
# remove existing connections
foreach my $mapped ( keys %rnai2genes ) {
    $mapped =~ s/_[ps]$//;
    print OUTACE "RNAi : $mapped\n";
    print OUTACE "-D Inhibits\n\n";
}

# open autoace connection with aceperl
my $db = Ace->connect( -path => $dbdir, -program => $tace ) || die "Couldn't connect to $dbdir\n", Ace->error;
my %inverse;
foreach my $mapped ( keys %rnai2genes ) {

    $mapped =~ /(.*)_([ps])$/;
    my $id   = $1;
    my $type =
      $2 eq 'p' ? 'primary' : 'secondary';   # $type holds whether the RNAi mapping to this gene is primary or secondary
    my $seq;
    my $worm_gene;                           # CDS, Transcript, or Pseudogene name

    foreach my $exon ( @{ $rnai2genes{$mapped} } ) {
	next if !(defined $exon);
	    
        print "$type RNAi '$mapped' is $type mapped to ", $exon->id, " (", $exon->type, ")\n" if ($verbose);

        if ( $type eq "primary" ) {
            $no_of_rnais_primary++;
        }
        else {
            $no_of_rnais_secondary++;
        }

        print OUTACE "\nRNAi : \"$id\"\n";

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
$wb->load_to_database($dbdir, $acefile,"RNAi_mappings") if ($load);

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
is then loaded into /wormsrv2/autoace

map_RNAi mandatory arguments:


=over 4

=item none

=back

map_RNAi optional arguments:

=over 4

=item -debug, Verbose/Debug mode

=item -test, Test mode, generate the acefile but do not upload them 

=item -help, Help pages

=item -acefile, write to a specific acefile

=item -expression, map only expression profiles (for debugging)

=item -load, loads file to autoace

=item -store specifiy configuration file

=back

=cut
