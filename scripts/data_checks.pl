#!/usr/local/bin/perl5.8.7 -w
#
# data_checks.pl
#
# by Keith Bradnam
#
# This is a example of a good script template
#
# Last updated by: $Author: pad $
# Last updated on: $Date: 2015-02-11 10:42:18 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

#use Ace;
#use Sequence_extract;
#use Coords_converter;

######################################
# variables and command-line options #
######################################

my ( $help, $debug, $test, $verbose, $store, $wormbase, $database );
my ( $ace,  $gff );

GetOptions(
    "help"       => \$help,
    "debug=s"    => \$debug,
    "test"       => \$test,
    "verbose"    => \$verbose,
    "store:s"    => \$store,
    "database:s" => \$database,
    "ace"        => \$ace,
    "gff"        => \$gff
);

if ($store) {
    $wormbase = retrieve($store) or croak("Can't restore wormbase from $store\n");
}
else {
    $wormbase = Wormbase->new(
        -debug => $debug,
        -test  => $test,
    );
}

# Display help if required
&usage("Help") if ($help);

# establish log file.
my $log  = Log_files->make_build_log($wormbase);
my $tace = $wormbase->tace;                        # TACE PATH

##########################
# MAIN BODY OF SCRIPT
##########################
$database ||= $wormbase->autoace;
my $autoace = Ace->connect( -path => $database ) or $log->log_and_die( Ace->error . "\n" );
my @queries;

# acedb queries, where the results of specified queries are compared against expected values
if ($ace) {
    my @queries = &read_acedb_queries;
    foreach my $query (@queries) {
	my($desc,$test,$expect) = @{$query};

        $log->write_to( "\nTEST (tace query): $test  ($desc )");
        my $count = $autoace->count( -query => $test );;
        $log->write_to(" . . . ok") if ( &pass_check( $expect, $count ) == 0 );
    }
}

# GFF queries, where the number of lines in GFF files is compared to what is expected based on the number of objects in the database.
if ($gff) {
    my $gff_dir = $database . "/CHROMOSOMES";
    &read_GFF_queries;
    foreach my $query (@queries) {
        $log->write_to( "\nTEST   (grep GFF): " . $$query{'GFF'} . ' (' . $$query{'DESC'}  . ')');
        my $actual = `cat $gff_dir/CHROMOSOME*.gff3 | grep -c \'$$query{'GFF'}\'`;
        my $expect = $$query{'EXPECT'} ? $$query{'EXPECT'} : $autoace->count( -query => $$query{'QUERY'} );
        $log->write_to(" . . . ok") if ( &pass_check( $expect, $actual ) == 0 );
    }
}

$log->mail();
exit(0);

##############################################################
#
# Subroutines
#
##############################################################

sub pass_check {
    my $expect = shift;
    my $count  = shift;
    if ( $expect == 0 ) {
        if ( $count != 0 ) {
            $log->write_to( "\nERROR: expect - " . $expect . ", actual $count\n" );
            $log->error;
            return 1;
        }
    }
    else {
        my $diff = $count / $expect;
        if ( $count == 0 or $diff > 1.1 or $diff < 0.9 ) {
            $log->write_to( "\nERROR: expect - " . $expect . ", actual $count\n" );
            $log->error;
            return 1;
        }
    }
    return 0;
}

sub read_acedb_queries { 
    my @queries;
    my $species = $wormbase->species;
    if($species eq 'elegans'){    
	@queries = (
	    ["The number of RNAi experiments with more than one associated Gene", 'find rnai COUNT gene > 1', 15734],
	    ["The number of RNAi results with connections to genes", 'find RNAi Gene', 104817],
	    ["The number of microarray results with connections to genes", 'find microarray_results gene', 340041],
	    ["PCR products overlapping CDS", "find PCR_product Overlaps_CDS", 62852],
	    ["The number of wormpep without pep_homol", 'find wormpep !pep_homol', 957],
	    ["tRNAs not attached to parent properly", 'Transcript AND NEXT AND NOT NEXT', 0],
	    ["Homol_data without waba", 'find Homol_data *waba !DNA_homol', 0],
	    ["Homol_data without Pep_homol", 'find Homol_data *wublastx* !Pep_homol', 0],
	    ["Inverted repeat Feature_data without features", 'find Feature_data *inverted !feature', 0],
	    ["TRF repeat Feature_data without features", 'find Feature_data *TRF !Feature', 0],
	    ["Oligo_sets with overlapping_CDS", 'find Oligo_Set Overlaps_CDS', 290113],
	    ["operons without genes", 'find operon !contains Gene', 0],
	    ["variation gene connection", 'find Variation Gene', 1205821],
	    ["genes with structured description", 'find Gene Structured_description', 137137],
	    ["genes with GO_annotation", 'find Gene GO_annotation', 12975],
	    ["CDSs with no source_exons", 'find CDS !Source_exons, method', 0],
	    ["Operons without parent ", 'find Operon CEO* !History AND !Canonical_parent',  0],
	    ["GO_term without Name or Definition", 'find GO_term !(Name or Definition)',  7],
	    ["Homol mapped Expression Patterns", 'find Expr_pattern where Homol_homol', 4506],
	    ["Transposon Objects mapped in the database", 'find Transposon where Sequence', 11256],
	    );
    }
    elsif( $species eq 'japonica'){  
	@queries = (
		   ["The number of wormpep without pep_homol", 'find protein JA* !pep_homol', 747],
		   ["tRNAs not attached to parent properly", 'Transcript AND NEXT AND NOT NEXT', 0],
		   ["Homol_data without Pep_homol", 'find Homol_data *wublastx* !Pep_homol', 0],
		   ["Inverted repeat Feature_data without features", 'find Feature_data *inverted !feature', 437],
		   ["TRF repeat Feature_data without features", 'find Feature_data *TRF !Feature', 0],
		   ["genes with GO_term", 'find Gene GO_term', 7708],
		   ["CDSs with no source_exons", 'find CDS CJA !Source_exons, method', 0],
		   ["Operons without parent ", 'find Operon !History AND !Canonical_parent',  0],
		   ["Proteins without peptides ", 'find Protein JA* !Peptide',  0],
		   ["CDSs without transcripts ", 'find CDS CJA* !Corresponding_transcript',  0],
		   );
    }    
    elsif( $species eq 'remanei'){  
	@queries = (
		   ["The number of wormpep without pep_homol", 'find protein RP* !pep_homol', 1362],
		   ["tRNAs not attached to parent properly", 'Transcript AND NEXT AND NOT NEXT', 0],
		   ["Homol_data without Pep_homol", 'find Homol_data *wublastx* !Pep_homol', 0],
		   ["Inverted repeat Feature_data without features", 'find Feature_data *inverted !feature', 437],
		   ["TRF repeat Feature_data without features", 'find Feature_data *TRF !Feature', 0],
		   ["genes with GO_term", 'find Gene GO_term', 7708],
		   ["CDSs with no source_exons", 'find CDS CRE !Source_exons, method', 0],
		   ["Operons without parent ", 'find Operon !History AND !Canonical_parent',  0],
		   ["Proteins without peptides ", 'find Protein RP*RP !Peptide',  0],
		   ["CDSs without transcripts ", 'find CDS CRE* !Corresponding_transcript',  0],
		   );
    }
    return @queries;
}



sub read_GFF_queries {
    @queries = ();
    my $i = 0;

    $queries[$i]{'DESC'}  = "ncRNA genes";
    $queries[$i]{'GFF'}   = "ncRNA\tncRNA"; # ncRNA genes including 21uRNA
    $queries[$i]{'QUERY'} = 'find elegans_RNA_genes; method = ncRNA OR method = 7kncRNA'; # this finds ncRNA genes including 21uRNA

    $i++;
    $queries[$i]{'DESC'}  = "tRNAs";
    $queries[$i]{'GFF'}   = "tRNA\ttRNA";
    $queries[$i]{'QUERY'} = 'find elegans_RNA_genes tRNA';
    
    $i++;
    $queries[$i]{'DESC'}  = "miRNAs";
    $queries[$i]{'GFF'}   = "miRNA_mature\tmiRNA";
    $queries[$i]{'QUERY'} = 'find elegans_RNA_genes; method = miRNA';

    $i++;
    $queries[$i]{'DESC'}  = "miRNAs - primary";
    $queries[$i]{'GFF'}   = "miRNA_primary_transcript\tmiRNA_primary_transcript";
    $queries[$i]{'QUERY'} = 'find elegans_RNA_genes; method = miRNA_primary_transcript';

    $i++;
    $queries[$i]{'DESC'}  = "snoRNAs";
    $queries[$i]{'GFF'}   = "snoRNA\tsnoRNA";
    $queries[$i]{'QUERY'} = 'find elegans_RNA_genes snoRNA';

    $i++;
    $queries[$i]{'DESC'}  = "snRNAs";
    $queries[$i]{'GFF'}   = "snRNA\tsnRNA";
    $queries[$i]{'QUERY'} = 'find elegans_RNA_genes snRNA';

    $i++;
    $queries[$i]{'DESC'}  = "rRNAs";
    $queries[$i]{'GFF'}   = "rRNA\trRNA";
    $queries[$i]{'QUERY'} = 'find elegans_RNA_genes rRNA';

    $i++;
    $queries[$i]{'DESC'}  = "scRNAs";
    $queries[$i]{'GFF'}   = "scRNA\tscRNA";
    $queries[$i]{'QUERY'} = 'find elegans_RNA_genes scRNA';

    $i++;
    $queries[$i]{'DESC'}  = "piRNAs";
    $queries[$i]{'GFF'}   = "piRNA\tpiRNA";
    $queries[$i]{'QUERY'} = 'find elegans_RNA_genes piRNA';

    $i++;
    $queries[$i]{'DESC'}  = "asRNAs";
    $queries[$i]{'GFF'}   = "asRNA\tantisense_RNA";
    $queries[$i]{'QUERY'} = 'find elegans_RNA_genes asRNA';

    $i++;
    $queries[$i]{'DESC'}  = "lincRNAs";
    $queries[$i]{'GFF'}   = "lincRNA\tlincRNA";
    $queries[$i]{'QUERY'} = 'find elegans_RNA_genes lincRNA';

    $i++;
    $queries[$i]{'DESC'}  = "stRNAs";
    $queries[$i]{'GFF'}   = "stRNA\tstRNA";
    $queries[$i]{'QUERY'} = 'find elegans_RNA_genes stRNA';

    $i++;
    $queries[$i]{'DESC'}  = "circRNAs";
    $queries[$i]{'GFF'}   = "circRNA\tcircular_ncRNA";
    $queries[$i]{'QUERY'} = 'find elegans_RNA_genes ncRNA AND method = "circRNA"';

    $i++;
    $queries[$i]{'DESC'}  = "Deletion and insertion alleles";
    $queries[$i]{'GFF'}   = "complex_substitution";
    $queries[$i]{'QUERY'} = 'find Variation flanking_sequences AND method = "Deletion_and_insertion_allele"';

    $i++;
    $queries[$i]{'DESC'}  = "Deletion alleles";
    $queries[$i]{'GFF'}   = "Allele\tdeletion";
    $queries[$i]{'QUERY'} = 'find Variation flanking_sequences AND method = "Deletion_allele"';

    $i++;
    $queries[$i]{'DESC'}  = "Substitution alleles";
    $queries[$i]{'GFF'}   = "\tsubstitution\t";
    $queries[$i]{'QUERY'} = 'find Variation flanking_sequences AND method = "Substitution_allele"';

    $i++;
    $queries[$i]{'DESC'}  = "Allele";
    $queries[$i]{'GFF'}   = "\tAllele\tsequence_alteration";
    $queries[$i]{'QUERY'} = 'find Variation Flanking_sequences AND (method = "Allele" OR method = "Engineered_allele")';

    $i++;
    $queries[$i]{'DESC'}  = "Transposon_insertion";
    $queries[$i]{'GFF'}   = "Allele\ttransposable_element_insertion_site\t";
    $queries[$i]{'QUERY'} = 'find Variation flanking_sequences AND method = "Transposon_insertion"';

    $i++;
    $queries[$i]{'DESC'}   = "RNAi primary locations";
    $queries[$i]{'GFF'}    = "RNAi_primary";
    $queries[$i]{'EXPECT'} = '285000';

    $i++;
    $queries[$i]{'DESC'}   = "RNAi secondary locations";
    $queries[$i]{'GFF'}    = "RNAi_secondary";
    $queries[$i]{'EXPECT'} = 75000;

    $i++;
    $queries[$i]{'DESC'}  = "Vancouver fosmids";
    $queries[$i]{'GFF'}   = "Vancouver_fosmid";
    $queries[$i]{'QUERY'} = 'find Sequence "WRM*" AND S_parent';    #bit of a cheat but much faster!

    $i++;
    $queries[$i]{'DESC'}  = "Coding_transcripts";
    $queries[$i]{'GFF'}   = "Coding_transcript\tmRNA";
    $queries[$i]{'QUERY'} = 'find Transcript Method=Coding_transcript';

    $i++;
    $queries[$i]{'DESC'}  = "All PCR products";
    $queries[$i]{'GFF'}   = "PCR_product";
    $queries[$i]{'QUERY'} = 'find PCR_product Canonical_parent';

    $i++;
    $queries[$i]{'DESC'}  = "Expression profiles";
    $queries[$i]{'GFF'}   = "Expr_profile";
    $queries[$i]{'QUERY'} = 'find Expr_profile S_parent';

    $i++;
    $queries[$i]{'DESC'}  = "cDNA for RNAi";
    $queries[$i]{'GFF'}   = "cDNA_for_RNAi";
    $queries[$i]{'QUERY'} = 'find RNAi Method; Homol_homol; follow Sequence';

# doesn't work as Oligo_sets don't map uniquely
#    $i++;
#    $queries[$i]{'DESC'}  = "mapped Oligo_set";
#    $queries[$i]{'GFF'}   = "Oligo_set";
#    $queries[$i]{'QUERY'} = 'find Oligo_set SMap';

    $i++;
    $queries[$i]{'DESC'}  = "mapped Operons";
    $queries[$i]{'GFF'}   = "operon";
    $queries[$i]{'QUERY'} = 'find Operon CEO*'; # small cheat

    $i++;
    $queries[$i]{'DESC'}  = "SL1 features";
    $queries[$i]{'GFF'}   = "SL1_acceptor_site";
    $queries[$i]{'QUERY'} = 'find Feature method = "SL1"';

    $i++;
    $queries[$i]{'DESC'}  = "SL2 features";
    $queries[$i]{'GFF'}   = "SL2_acceptor_site";
    $queries[$i]{'QUERY'} = 'find Feature method = "SL2"';

    $i++;
    $queries[$i]{'DESC'}  = "polyA sites";
    $queries[$i]{'GFF'}   = "polyA_site";
    $queries[$i]{'QUERY'} = 'find Feature method = "polyA_site"';

    $i++;
    $queries[$i]{'DESC'}  = "Non_coding_transcript isoforms";
    $queries[$i]{'GFF'}   = "non_coding_transcript\tnc_primary_transcript";
    $queries[$i]{'QUERY'} = 'find Transcript; method = non_coding_transcript';

    $i++;
    $queries[$i]{'DESC'}   = "Mass Spectrometry peptides";
    $queries[$i]{'GFF'}    = "mass_spec_genome";
    $queries[$i]{'EXPECT'} = 136433;
    
    $i++;
    $queries[$i]{'DESC'}   = "Expr_pattern mapped Expression Patterns";
    $queries[$i]{'GFF'}    = "Expr_pattern\treagent";
    $queries[$i]{'EXPECT'} = 2934;

    $i++;
    $queries[$i]{'DESC'}   = "Chronogram mappings";
    $queries[$i]{'GFF'}    = "Chronogram\treagent";
    $queries[$i]{'EXPECT'} = 1974;    
    
    #	$i++;
    #	$queries[$i]{'DESC'}  = "";
    #	$queries[$i]{'GFF'}   = "";
    #	$queries[$i]{'QUERY'} = '';

}

# Add perl documentation in POD format
# This should expand on your brief description above and
# add details of any options that can be used with the program.
# Such documentation can be viewed using the perldoc command.

__END__

=pod

=head2 NAME - data_checks.pl

=head1 USAGE

=over 4

=item data_checks.pl  -ace -gff

=back

This script runs checks to ensure that data is correct. 

There are two types of checks; one that runs a query and checks the number of objects returned and 
one that compares the number of objects of a type with how many are in the GFF files.



script_template.pl MANDATORY arguments:

=over 4

=item None at present.

=back

data_checks.pl  OPTIONAL arguments:

=over 4

=item -ace

runs defined queries against the specified database and compares results to expected values

=back

=over 4

=item -gff

compare the number of GFF lines meeting certain criteria with number expected from database query eg the number of operons should be the same in the GFF as acedb

=back

=over 4

=item -database

which database to run these checks against.  If -gff requires chromosome gff files in CHROMOSOMES dir

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode uses TEST_BUILD

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

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
