#!/usr/bin/env perl
#
# genestats.pl
#
# dl
#
# Usage : genestats.pl [-options]
#
# Last edited by: $Author: klh $
# Last edited on: $Date: 2013-10-10 09:56:39 $
 

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use IO::Handle;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);

GetOptions ('help'       => \$help,
            'debug=s'    => \$debug,
	    'test'       => \$test,
	    'verbose'    => \$verbose,
	    'store:s'    => \$store,
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

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

#################################
# Set up some useful paths      #
#################################

# Set up top level base directories (these are different if in test mode)
my $ace_dir         = $wormbase->autoace;     # AUTOACE DATABASE DIR
my $tace            = $wormbase->tace;        # TACE PATH
my $reports_dir     = $wormbase->reports;     # AUTOACE REPORTS


#######################
# get data from ACEDB #
#######################

my $coding_cmd_prefix    = "query find Gene where Species = \"*elegans\" AND Live AND Corresponding_CDS\nspush\n";
my $noncoding_cmd_prefix = "query find Gene where Species = \"*elegans\" AND Live AND NOT Corresponding_CDS AND Molecular_info\nspush\n";
my $uncloned_cmd_prefix = "query find Gene where Species = \"*elegans\" AND Live AND NOT Molecular_info\nspush\n";

my $rest = "query where Concise_description\n"
    . "spop\nspush\n"
    . "query where Automated_description\n"
    . "spop\nspush\n"
    . "query where Disease_info\n"
    . "spop\nspush\n"
    . "query where Reference\n"
    . "spop\nspush\n"
    . "query where CGC_name\n"
    . "spop\nspush\n"
    . "query where RNAi_result\n"
    . "spop\nspush\n"
    . "query where Expr_pattern\n"
    . "spop\nspush\n"
    . "query where Interaction\n"
    . "spop\nspush\n"
    . "query where Microarray_results\n"
    . "spop\nspush\n"
    . "query where Allele\n"
    . "quit\n";

open (my $out_fh, ">$reports_dir/genedata") || die "Failed to open output file\n";

# database connection
my $coding_query = $coding_cmd_prefix . $rest;
my $nc_query = $noncoding_cmd_prefix . $rest;
my $uncloned_query = $uncloned_cmd_prefix . $rest;

my $values1 = &tace_query($coding_query);
my $values2 = &tace_query($nc_query);
my $values3 = &tace_query($uncloned_query);


printf $out_fh "\nC. elegans gene data (%d genes in total)\n", $values1->[0] + $values2->[0] + $values3->[0];
print $out_fh "----------------------------------------------\n";

printf $out_fh "\nProtein-coding (%d genes):\n", $values1->[0];;
&report_gene_data($out_fh, $values1);
printf $out_fh "\nNon-coding RNA and pseudogene (%d genes):\n", $values2->[0];
&report_gene_data($out_fh, $values2);
printf $out_fh "\nUncloned (%d genes):\n", $values3->[0];
&report_gene_data($out_fh, $values3);

$log->mail();


###########################################################

sub tace_query {
  my ($q) = @_;

  my @results;

  open (my $tace_fh,"echo '$q' | $tace $ace_dir |");
  while (<$tace_fh>) {
    chomp;
    if (/\/\/ Found (\d+) objects/) {
      push (@results,$1);
    }
  }
  close($tace_fh);

  return \@results;
}


sub report_gene_data {
  my ($fh, $values) = @_;

  my $live_genes          = $values->[0];
  my $concise_desc        = $values->[1];
  my $auto_desc           = $values->[2];
  my $disease_assoc       = $values->[3];
  my $reference           = $values->[4];
  my $CGC_name            = $values->[5];
  my $RNAi_result         = $values->[6];
  my $expr_pat            = $values->[7];
  my $interactions        = $values->[8];
  my $microarray_result   = $values->[9];
  my $variations          = $values->[10];
  
  my $percent_concise_desc        = sprintf("(%2.1f%%)", ($concise_desc / $live_genes) * 100);
  my $percent_auto_desc           = sprintf("(%2.1f%%)", ($auto_desc / $live_genes) * 100);
  my $percent_disease_assoc       = sprintf("(%2.1f%%)", ($disease_assoc / $live_genes) * 100);
  my $percent_reference           = sprintf("(%2.1f%%)", ($reference / $live_genes) * 100);
  my $percent_CGC_name            = sprintf("(%2.1f%%)", ($CGC_name / $live_genes) * 100);
  my $percent_RNAi_result         = sprintf("(%2.1f%%)", ($RNAi_result / $live_genes) * 100);
  my $percent_expr_pat            = sprintf("(%2.1f%%)", ($expr_pat / $live_genes) * 100);
  my $percent_interaction         = sprintf("(%2.1f%%)", ($interactions / $live_genes) * 100);
  my $percent_microarray_result   = sprintf("(%2.1f%%)", ($microarray_result / $live_genes) * 100);
  my $percent_variation           = sprintf("(%2.1f%%)", ($variations / $live_genes) * 100);

  printf $fh "  Curated description        %5d   %8s\n",  $concise_desc, $percent_concise_desc;
  printf $fh "  Automated description      %5d   %8s\n",  $auto_desc, $percent_auto_desc;
  printf $fh "  Human disease association  %5d   %8s\n",  $disease_assoc, $percent_disease_assoc;
  printf $fh "  Approved Gene name         %5d   %8s\n",  $CGC_name, $percent_CGC_name;
  printf $fh "  Reference                  %5d   %8s\n",  $reference, $percent_reference;
  printf $fh "  RNAi results               %5d   %8s\n",  $RNAi_result, $percent_RNAi_result;
  printf $fh "  Microarray results         %5d   %8s\n",  $microarray_result, $percent_microarray_result;
  printf $fh "  Expression patterns        %5d   %8s\n",  $expr_pat,  $percent_expr_pat;
  printf $fh "  Variations                 %5d   %8s\n",  $variations, $percent_variation;
  printf $fh "  Interaction data           %5d   %8s\n",  $interactions,  $percent_interaction;
}



__END__

