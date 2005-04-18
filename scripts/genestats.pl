#!/usr/local/bin/perl
#
# genestats.pl
#
# dl
#
# Usage : genestatsr.pl [-options]
#
# Last edited by: $Author: dl1 $
# Last edited on: $Date: 2005-04-18 09:54:41 $
 
use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Getopt::Long;

our $database   = "/wormsrv2/autoace";
our $tace       =  &tace;

#######################
# get data from ACEDB #
#######################

my @values;
my $command;

$command  = "query find Gene where Species = \"*elegans\" AND Live\nspush\n";
$command .= "query where Molecular_info\n";
$command .= "spop\nspush\n";
$command .= "query where Concise_description\n";
$command .= "spop\nspush\n";
$command .= "query where Reference\n";
$command .= "spop\nspush\n";
$command .= "query where CGC_name\n";
$command .= "spop\nspush\n";
$command .= "query where RNAi_result\n";
$command .= "spop\nspush\n";
$command .= "query where Microarray_results\n";
$command .= "spop\nspush\n";
$command .= "query where SAGE_transcript\n";
$command .= "quit\n";

# database connection

open (TACE,"echo '$command' | $tace $database |");
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
my $reference           = $values[3];
my $CGC_name            = $values[4];
my $RNAi_result         = $values[5];
my $microarray_result   = $values[6];
my $SAGE_transcript     = $values[7];

my $percent_molecular_info      = (int ( ( ($molecular_info / $live_genes) * 1000) + 0.5) / 10); 
my $percent_concise_description = (int ( ( ($concise_description / $live_genes) * 1000) + 0.5) / 10); 
my $percent_reference           = (int ( ( ($reference / $live_genes) * 1000) + 0.5) / 10); 
my $percent_CGC_name            = (int ( ( ($CGC_name / $live_genes) * 1000) + 0.5) / 10); 
my $percent_RNAi_result         = (int ( ( ($RNAi_result / $live_genes) * 1000) + 0.5) / 10); 
my $percent_microarray_result   = (int ( ( ($microarray_result / $live_genes) * 1000) + 0.5) / 10); 
my $percent_SAGE_transcript     = (int ( ( ($SAGE_transcript / $live_genes) * 1000) + 0.5) / 10); 

##################
# report to file #
##################

open (OUT, ">/wormsrv2/autoace/REPORTS/genedata") || die "Failed to open output file\n";
print OUT "Gene data set (Live C.elegans genes $values[0])\n";
print OUT "------------------------------------------\n";
print OUT "Molecular_info              "  . $values[1] . " (" . $percent_molecular_info . "%)\n";
print OUT "Concise_description          " . $values[2] . " (" . $percent_concise_description . "%)\n";
print OUT "Reference                    " . $values[3] . " (" . $percent_reference . "%)\n";
print OUT "CGC_approved Gene name       " . $values[4] . " (" . $percent_CGC_name . "%)\n";
print OUT "RNAi_result                 "  . $values[5] . " (" . $percent_RNAi_result . "%)\n";
print OUT "Microarray_results          "  . $values[6] . " (" . $percent_microarray_result . "%)\n";
print OUT "SAGE_transcript             "  . $values[7] . " (" . $percent_SAGE_transcript . "%)\n";
close OUT;

# hasta luego
exit(0);

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
write the report to /wormsrv2/autoace/REPORTS.

=back

=over 4

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Daniel Lawson (dl1@sanger.ac.uk)

=back

=cut
