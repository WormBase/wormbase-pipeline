#!/usr/bin/env perl
#
# genestats.pl
#
# dl
#
# Usage : genestats.pl [-options]
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2012-06-27 12:20:59 $
 

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
          . "query where Human_disease_relevance\n"
          . "spop\nspush\n"
          . "query where Reference\n"
          . "spop\nspush\n"
          . "query where CGC_name\n"
          . "spop\nspush\n"
          . "query where RNAi_result\n"
          . "spop\nspush\n"
          . "query where Microarray_results\n"
          . "spop\nspush\n"
          . "query where SAGE_tag\n"
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
my $concise_description = $values[2];
my $human_disease_rel   = $values[3];
my $reference           = $values[4];
my $CGC_name            = $values[5];
my $RNAi_result         = $values[6];
my $microarray_result   = $values[7];
my $SAGE_transcript     = $values[8];

my $percent_molecular_info      = (int ( ( ($molecular_info / $live_genes) * 1000) + 0.5) / 10); 
my $percent_concise_description = (int ( ( ($concise_description / $live_genes) * 1000) + 0.5) / 10);
my $percent_human_disease_rel   = (int ( ( ($human_disease_rel / $live_genes) * 1000) + 0.5) / 10);
my $percent_reference           = (int ( ( ($reference / $live_genes) * 1000) + 0.5) / 10); 
my $percent_CGC_name            = (int ( ( ($CGC_name / $live_genes) * 1000) + 0.5) / 10); 
my $percent_RNAi_result         = (int ( ( ($RNAi_result / $live_genes) * 1000) + 0.5) / 10); 
my $percent_microarray_result   = (int ( ( ($microarray_result / $live_genes) * 1000) + 0.5) / 10); 
my $percent_SAGE_transcript     = (int ( ( ($SAGE_transcript / $live_genes) * 1000) + 0.5) / 10); 

##################
# report to file #
##################
open (OUT, ">$reports_dir/genedata") || die "Failed to open output file\n";
print OUT "Gene data set (Live C. elegans genes $values[0])\n";
print OUT "------------------------------------------\n";
print OUT "Molecular_info              " . $values[1] . " (" . $percent_molecular_info . "%)\n";
print OUT "Concise_description         " . $values[2] . "  (" . $percent_concise_description . "%)\n";
print OUT "Human_disease_relevance     " . $values[3] . "    (" . $percent_human_disease_rel . "%)\n";
print OUT "Reference                   " . $values[4] . " (" . $percent_reference . "%)\n";
print OUT "WormBase_approved_Gene_name " . $values[5] . " (" . $percent_CGC_name . "%)\n";
print OUT "RNAi_result                 " . $values[6] . " (" . $percent_RNAi_result . "%)\n";
print OUT "Microarray_results          " . $values[7] . " (" . $percent_microarray_result . "%)\n";
print OUT "SAGE_transcript             " . $values[8] . " (" . $percent_SAGE_transcript . "%)\n";
close OUT;

##################
# Check the files
##################

$wormbase->check_file("$reports_dir/genedata", $log,
minlines => 10,
maxlines => 10,
line1 => '^Gene data set \(Live C. elegans genes \d+\)',
line2 => '^\-+',
lines => ['^\S+\s+\d+\s+\(\d+(\.\d)*\%\)', '^\S+\s+\S+\s+\S+\s+\d+\s+\(\d+(\.\d)*\%\)'],
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
