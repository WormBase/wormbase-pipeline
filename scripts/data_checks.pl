#!/usr/local/bin/perl5.8.7 -w
#
# data_checks.pl                        
# 
# by Keith Bradnam                         
#
# This is a example of a good script template
#
# Last updated by: $Author: ar2 $     
# Last updated on: $Date: 2006-04-06 13:47:05 $      

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

my ($help, $debug, $test, $verbose, $store, $wormbase, $database);
my ($ace, $gff);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
		      "test"       => \$test,
	    		"verbose"    => \$verbose,
	     		"store:s"    => \$store,
	     		"database:s" => \$database,
	     		"ace"        => \$ace,
	     		"gff"        => \$gff
	    	);

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# establish log file.
my $log = Log_files->make_build_log($wormbase);
my $tace            = $wormbase->tace;        # TACE PATH

##########################
# MAIN BODY OF SCRIPT
##########################
$database = $database or $wormbase->autoace;
my $autoace = Ace->connect( -path => $database) or $log->log_and_die(Ace->error."\n");
my @queries;

# acedb queries, where the results of specified queries are compared against expected values
if( $ace ) {
	&read_acedb_queries;
	foreach my $query (@queries) {
		$log->write_to("\nTEST: ".$$query{'QUERY'});
		my $count= $autoace->count(-query => $$query{'QUERY'});
		my $expect = $$query{'RESULT'};
		$log->write_to(" . . . ok") if (&pass_check($expect, $count) == 0);
	}
}

# GFF queries, where the number of lines in GFF files is compared to what is expected based on the number of objects in the database.
if( $gff ){
	my $gff_dir = $database."/CHROMOSOMES";
	&read_GFF_queries;
	foreach my $query (@queries) {
		$log->write_to("\nTEST: ".$$query{'DESC'});
		my $actual = `cat $gff_dir/CHROMOSOME*.gff | grep -c $$query{'GFF'}`;
		my $expect = $$query{'EXPECT'} ?  $$query{'EXPECT'} : $autoace->count(-query => $$query{'QUERY'});
		$log->write_to(" . . . ok") if (&pass_check($expect, $actual) == 0);
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
	if( $expect ==0 ) {
		if( $count != 0 ) {
			$log->write_to("\nERROR: expect - ".$expect.", actual $count\n");
			$log->error;
			return 1;
		}
	}
	else {
		my $diff = $count / $expect;
		if ($count == 0 or $diff > 1.1 or $diff < 0.9 ){
			$log->write_to("\nERROR: expect - ".$expect.", actual $count\n");
			$log->error;
			return 1;
		}
	}	
	return 0;
}


sub read_acedb_queries {

	my $i = 0;
	$queries[$i]{'DESC'}   = "The number of RNAi experiments with more than one associated Gene";
	$queries[$i]{'QUERY'}  = 'find rnai COUNT gene > 1 AND uniquely_mapped';
	$queries[$i]{'RESULT'} = 1716;

	$i++;
	$queries[$i]{'DESC'}   = "The number of RNAi results with connections to genes";
	$queries[$i]{'QUERY'}  = 'find RNAi Gene';
	$queries[$i]{'RESULT'} = 61844 ;
	
	$i++;
	$queries[$i]{'DESC'}   = "The number of microarray results with connections to genes";
	$queries[$i]{'QUERY'}  = 'find microarray_results gene';
	$queries[$i]{'RESULT'} = 58800 ;
			
	$i++;
	$queries[$i]{'DESC'}   = "PCR products overlapping CDS";
	$queries[$i]{'QUERY'}  = "find PCR_product Overlaps_CDS";
	$queries[$i]{'RESULT'} = 62852 ;
			
	$i++;
	$queries[$i]{'DESC'}   = "The number of wormpep without pep_homol";
	$queries[$i]{'QUERY'}  = 'find wormpep !pep_homol';
	$queries[$i]{'RESULT'} = 45 ;

	$i++;
	$queries[$i]{'DESC'}   = "tRNAs not attached to parent properly";
	$queries[$i]{'QUERY'}  = 'Transcript AND NEXT AND NOT NEXT';
	$queries[$i]{'RESULT'} = 0 ;
	
	$i++;
	$queries[$i]{'DESC'}   = "Homol_data without waba";
	$queries[$i]{'QUERY'}  = 'find Homol_data *waba !DNA_homol';
	$queries[$i]{'RESULT'} = 0 ; 
	
	$i++;
	$queries[$i]{'DESC'}   = "Homol_data without Pep_homol";
	$queries[$i]{'QUERY'}  = 'find Homol_data *wublastx* !Pep_homol';
	$queries[$i]{'RESULT'} = 5207 ; 			
	$i++;
	
	$queries[$i]{'DESC'}   = "Inverted repeat Feature_data without features";
	$queries[$i]{'QUERY'}  = 'find Feature_data *inverted !feature';
	$queries[$i]{'RESULT'} = 505 ; 
	
	$i++;
	$queries[$i]{'DESC'}   = "TRF repeat Feature_data without features";
	$queries[$i]{'QUERY'}  = 'find Feature_data *TRF !Feature';
	$queries[$i]{'RESULT'} = 0 ; 
		
	$i++;
	$queries[$i]{'DESC'}   = "Oligo_sets with overlapping_CDS";
	$queries[$i]{'QUERY'}  = 'find Oligo_Set Overlaps_CDS';
	$queries[$i]{'RESULT'} = 74615 ; 
		
	$i++;
	$queries[$i]{'DESC'}   = "operons without genes";
	$queries[$i]{'QUERY'}  = 'find operon !contains Gene';
	$queries[$i]{'RESULT'} = 0 ; 

	$i++;
	$queries[$i]{'DESC'}   = "variation gene connection";
	$queries[$i]{'QUERY'}  = 'find Variation Gene';
	$queries[$i]{'RESULT'} =  21489; 

	$i++;
	$queries[$i]{'DESC'}   = "genes with structured description";
	$queries[$i]{'QUERY'}  = 'find Gene Structured_description';
	$queries[$i]{'RESULT'} = 4157 ; 
	
	$i++;
	$queries[$i]{'DESC'}   = "genes with GO_term";
	$queries[$i]{'QUERY'}  = 'find Gene GO_term';
	$queries[$i]{'RESULT'} = 9897 ; 

#	$i++;
#	$queries[$i]{'DESC'}   = "";
#	$queries[$i]{'QUERY'}  = '';
#	$queries[$i]{'RESULT'} =  ; 	
}

sub read_GFF_queries {
	@queries = ();
	my $i = 0;
	
	$queries[$i]{'DESC'}  = "Deletion and insertion alleles";
	$queries[$i]{'GFF'}   = "complex_change_in_nucleotide_sequence";
	$queries[$i]{'QUERY'} = 'find Variation flanking_sequences AND method = "Deletion_and_insertion_allele"';

	$i++;	
	$queries[$i]{'DESC'}  = "Deletion alleles";
	$queries[$i]{'GFF'}   = "deletion";
	$queries[$i]{'QUERY'} = 'find Variation flanking_sequences AND method = "Deletion_allele"';

	$i++;	
	$queries[$i]{'DESC'}  = "Substitution alleles";
	$queries[$i]{'GFF'}   = "substitution";
	$queries[$i]{'QUERY'} = 'find Variation flanking_sequences AND method = "Substitution_allele"';

	$i++;	
	$queries[$i]{'DESC'}  = "RNAi primary";
	$queries[$i]{'GFF'}   = "RNAi_primary";
	$queries[$i]{'EXPECT'}= '134246';

	$i++;	
	$queries[$i]{'DESC'}  = "RNAi secondary";
	$queries[$i]{'GFF'}   = "RNAi_secondary";
	$queries[$i]{'EXPECT'}= '14165';
	
	$i++;	
	$queries[$i]{'DESC'}  = "Alleles";
	$queries[$i]{'GFF'}   = "sequence_variant";
	$queries[$i]{'QUERY'} = 'find Variation flanking_sequences AND method = "Allele"';

	$i++;
	$queries[$i]{'DESC'}  = "Vancouver fosmids";
	$queries[$i]{'GFF'}   = "Vancouver_fosmid";
	$queries[$i]{'QUERY'} = 'find Sequence "WRM*"'; #bit of a cheat but much faster!
		
	$i++;
	$queries[$i]{'DESC'}  = "Coding_transcripts";
	$queries[$i]{'GFF'}   = "protein_coding_primary_transcript";
	$queries[$i]{'QUERY'} = 'find Coding_transcripts"';
			
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
					
	$i++;
	$queries[$i]{'DESC'}  = "mapped Oligo_set";
	$queries[$i]{'GFF'}   = "Oligo_set";
	$queries[$i]{'QUERY'} = 'find Oligo_set';
						
	$i++;
	$queries[$i]{'DESC'}  = "mapped Operons";
	$queries[$i]{'GFF'}   = "operon";
	$queries[$i]{'QUERY'} = 'find Operon';
							
	$i++;
	$queries[$i]{'DESC'}  = "SL1";
	$queries[$i]{'GFF'}   = "SL1_acceptor_site";
	$queries[$i]{'QUERY'} = 'find Feature method = "SL1"';
							
	$i++;
	$queries[$i]{'DESC'}  = "SL2";
	$queries[$i]{'GFF'}   = "SL2_acceptor_site";
	$queries[$i]{'QUERY'} = 'find Feature method = "SL2"';
								
	$i++;
	$queries[$i]{'DESC'}  = "polyA sites";
	$queries[$i]{'GFF'}   = "polyA_site";
	$queries[$i]{'QUERY'} = 'find Feature method = "polyA_site"';
								
	$i++;
	$queries[$i]{'DESC'}  = "Non_coding_transcripts";
	$queries[$i]{'GFF'}   = "Non_coding_transcript";
	$queries[$i]{'QUERY'} = 'find Transcript; method = non_coding_transcript';
															
	$i++;
	$queries[$i]{'DESC'}  = "tRNAs";
	$queries[$i]{'GFF'}   = "tRNA_primary_transcript";
	$queries[$i]{'QUERY'} = 'find elegans_RNA_genes tRNA';
															
	$i++;
	$queries[$i]{'DESC'}  = "miRNAs";
	$queries[$i]{'GFF'}   = "miRNA_primary_transcript";
	$queries[$i]{'QUERY'} = 'find elegans_RNA_genes miRNA';

															
	$i++;
	$queries[$i]{'DESC'}  = "snoRNAs";
	$queries[$i]{'GFF'}   = "snoRNA_primary_transcript";
	$queries[$i]{'QUERY'} = 'find elegans_RNA_genes snoRNA';

																							
	$i++;
	$queries[$i]{'DESC'}  = "snRNAs";
	$queries[$i]{'GFF'}   = "snRNA_primary_transcript";
	$queries[$i]{'QUERY'} = 'find elegans_RNA_genes snRNA';

																							
	$i++;
	$queries[$i]{'DESC'}  = "rRNAs";
	$queries[$i]{'GFF'}   = "rRNA_primary_transcript";
	$queries[$i]{'QUERY'} = 'find elegans_RNA_genes rRNA';

																							
	$i++;
	$queries[$i]{'DESC'}  = "scRNAs";
	$queries[$i]{'GFF'}   = "scRNA_primary_transcript";
	$queries[$i]{'QUERY'} = 'find elegans_RNA_genes scRNA';

																							
	$i++;
	$queries[$i]{'DESC'}  = "stRNAs";
	$queries[$i]{'GFF'}   = "stRNA_primary_transcript";
	$queries[$i]{'QUERY'} = 'find elegans_RNA_genes stRNA';

																							
	$i++;
	$queries[$i]{'DESC'}  = "ncRNAs";
	$queries[$i]{'GFF'}   = "ncRNA_primary_transcript";
	$queries[$i]{'QUERY'} = 'find elegans_RNA_genes ncRNA';

								
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

script_template.pl MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

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
