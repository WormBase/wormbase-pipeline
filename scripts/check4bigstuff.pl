#!/software/bin/perl -w
#
# check4bigstuff.pl
#
# by Paul Davis
#
# Script to look for oversized objects.
#
# Last updated by: $Author: pad $
# Last updated on: $Date: 2011-11-11 17:14:21 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Storable;
use Log_files;


######################################
# variables and command-line options #
######################################

my ($help, $debug, $test, $store, $file, $term, $size, $gff_dir);

GetOptions (
            'help'         => \$help,   # help documentation.
            'test'         => \$test,   # test build
            'debug=s'      => \$debug,  # debug option for email
	    'store:s'      => \$store,
	    'file:s'       => \$file,   # check a single file
	    'term:s'       => \$term,   # option to restrict to a specific term
	    'size:s'       => \$size,   # how big to allow things to be
	    'gff_dir:s'    => \$gff_dir, # an alternate gff_dir to the build/test_build 
           );

my $wormbase;
if ($store) {
  $wormbase = retrieve($store) or croak ("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug => $debug,
                             -test => $test,
                           );
}


# establish log file.
my $log = Log_files->make_build_log($wormbase);


##########################
# MAIN BODY OF SCRIPT
##########################
my $WS_name         = $wormbase->get_wormbase_version_name(); # e.g. WS132
my $flag;

unless (defined $gff_dir) {
  $gff_dir         = $wormbase->gff;         # AUTOACE GFF
}

if (!defined $size) {
$size = "100000";
}

if (!defined $term){
$term = "CHROM";
}

my $gff_file = "/nfs/wormpub/tmp/${WS_name}_partial_${term}.gff";

# Probe the build gff data.
if (!defined $file){
#  if (!-e ${gff_dir}."/CHROMOSOME_I.gff"){
 #   $log->log_and_die("Can't open input files\n");
  #}

  # Set a flag to remove the temp file
  $flag = "";

  # extract data but get rid of Link objects, large "region" lines, genetic map anomalies, CGH alleles, ttn-1, nhr-27 which contains a very large RST confirmed intron and a few other known biggies.
  $wormbase->run_command("cat $gff_dir/*.gff | grep $term | grep -v region | grep -v map | grep -v CGH_allele | grep -v W06H8.8 | grep -v ttn | grep -v nhr-27 | grep -v F16H9.2b | grep -v CEOP1017 | grep -v WBsf041392 | grep -v WBsf041396 | grep -v mGene_pred_17713 | grep -v WBGene00006436 | grep -v  WBGene00008901 | grep -v Aff_Y116F11.ZZ33 | grep -v Aff_Y105E8A.H > $gff_file", $log);
  $file = "$gff_file";
}

$log->write_to ("Checking $file for \"Big Stuff\" anything larger than $size will be flagged\nttn-1 and other known large objects have been pre-filtered from the input data\n----------------------------------------\n");

if (!-e $file){
  $log->log_and_die("Can't open input file\n");
}

open (IN,  "<$file") or die "Cannot open $file \n";
my $count = "0";
while (<IN>) {
  chomp;
  $count++;
  #CHROMOSOME_I    gene    gene    10413   16842   .       +       .       Gene "WBGene00022276"
  if (/^\S+\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s+\S+/) {
    my $span = $2 - $1;
    if ($span > $size) {
      $log->write_to ("WARNING: Line $count contains a span (${span}) greater than $size [$_]\n");
    }
    else {
      next;
    }
  }
  else {$log->write_to ("ERROR: $_ does not appear to be a properly formatted gff line.\n");
      }
}

if (defined $flag) {
  $wormbase->run_command("rm $file", $log);
}

if ($count < 500) {
  $log->write_to ("Warning: Only checked $count Lines of gff, were you expecting more\n");
}
else {
  $log->write_to ("Checked $count Lines of gff\n");
}

$log->mail;
print "Finished.\n" if ($debug);
close IN;
exit(0);

__END__

=pod

=head2 NAME - check4bigstuff.pl

=over 3

=head1 USAGE

=over 4

=item check4bigstuff.pl  [-options]

=back

This script looks for oversized objects in the GFF files or the entire chromosome_*.gff files.

Check4bigstuff.pl MANDATORY arguments:

=over 4

=item None at present.

=back

Check4bigstuff.pl MANDATORY arguments:

=over 4

=item None at present.

=back

Check4bigstuff.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4

=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.

=back

=over 4

=item -gff_dir, Specify a directory containing CHROMOSOME*.gff files eg. you can check the previous build release.

=back

=over 4

=item -size, Specify the maximum allowable length of objects in the database.

=back

=over 4

=item -term, This allows you to specify a data type based on a GFF_source/GFF_method or a specific ID type.

=back

=over 4

=item -file, This option allows you to specify a single file as apposed to checking the CHROMOSOME.gff files in a build/defined directory.

=back

=over 4

=item -test, Test mode, run the script, and checks the test build gff data by defult.

=back

=over 4

=item -verbose, output extra info.

=back

=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Paul Davis (pad@sanger.ac.uk)

=back

=cut
