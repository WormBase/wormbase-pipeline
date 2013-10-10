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

my @values;
my $command;

$command  = "query find Gene where Species = \"*elegans\" AND Live\nspush\n"
          . "query where Molecular_info\n"
          . "spop\nspush\n"
          . "query where Concise_description\n"
          . "spop\nspush\n"
          . "query where GO_term\n"
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

# database connection

open (TACE,"echo '$command' | $tace $ace_dir |");
while (<TACE>) {
    chomp;
    if (/\/\/ Found (\d+) objects/) {
	push (@values,$1);
    }
}
close TACE;

###########################################################
# assign values (this is verbose for ease of maintenance) #
###########################################################

my $live_genes          = $values[0];
my $molecular_info      = $values[1];
my $concise_desc        = $values[2];
my $go_term             = $values[3];
my $disease_assoc       = $values[4];
my $reference           = $values[5];
my $CGC_name            = $values[6];
my $RNAi_result         = $values[7];
my $expr_pat            = $values[8];
my $interactions        = $values[9];
my $microarray_result   = $values[10];
my $variations          = $values[11];

my $percent_molecular_info      = ($molecular_info / $live_genes) * 100;
my $percent_concise_desc        = ($concise_desc / $live_genes) * 100;
my $percent_go_term             = ($go_term / $live_genes) * 100;
my $percent_disease_assoc       = ($disease_assoc / $live_genes) * 100;
my $percent_reference           = ($reference / $live_genes) * 100;
my $percent_CGC_name            = ($CGC_name / $live_genes) * 100;
my $percent_RNAi_result         = ($RNAi_result / $live_genes) * 100;
my $percent_expr_pat            = ($expr_pat / $live_genes) * 100;
my $percent_interaction         = ($interactions / $live_genes) * 100;
my $percent_microarray_result   = ($microarray_result / $live_genes) * 100;
my $percent_variation           = ($variations / $live_genes) * 100;

##################
# report to file #
##################
open (OUT, ">$reports_dir/genedata") || die "Failed to open output file\n";
print OUT "Gene data set (Live C. elegans genes $values[0])\n";
print OUT "------------------------------------------\n";
printf OUT "Molecular info             %5d   (%2.1f%%)\n",  $molecular_info, $percent_molecular_info;
printf OUT "Concise description        %5d   (%2.1f%%)\n",  $concise_desc, $percent_concise_desc;
printf OUT "Human disease association  %5d   (%2.1f%%)\n",  $disease_assoc, $percent_disease_assoc;
printf OUT "Approved Gene name         %5d   (%2.1f%%)\n",  $CGC_name, $percent_CGC_name;
printf OUT "Reference                  %5d   (%2.1f%%)\n",  $reference, $percent_reference;
printf OUT "RNAi results               %5d   (%2.1f%%)\n",  $RNAi_result, $percent_RNAi_result;
printf OUT "Microarray results         %5d   (%2.1f%%)\n",  $microarray_result, $percent_microarray_result;
printf OUT "Expression patterns        %5d   (%2.1f%%)\n",  $expr_pat,  $percent_expr_pat;
printf OUT "Variations                 %5d   (%2.1f%%)\n",  $variations, $percent_variation;
printf OUT "Interaction data           %5d   (%2.1f%%)\n",  $interactions,  $percent_interaction;
close OUT;

##################
# Check the files
##################

$wormbase->check_file("$reports_dir/genedata", $log,
minlines => 12,
maxlines => 12,
line1 => '^Gene data set \(Live C. elegans genes \d+\)',
line2 => '^\-+',
lines => ['^\S+.+\s+\d+\s+\(\d+(\.\d)*\%\)', '^\S+\s+\S+\s+\S+\s+\d+\s+\(\d+(\.\d)*\%\)'],
);


$log->mail();
print "Finished.\n" if ($verbose);


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
  }
}

##########################################




__END__

=pod

=head1 NAME 

=over 4

=item genestats.pl

=back

=head1 USAGE

=over 4

=item genestats.pl 

=back

Queries autoace to extract a number of counts for the Live gene class. Computes percentage values to 1 decimel place and
write the report to autoace/REPORTS.

=back

=over 4

=head1 REQUIREMENTS

=over 4

=item None.

=back

=head1 AUTHOR

=over 4

=item Daniel Lawson (dl1@sanger.ac.uk)

=back

=cut
