#!/usr/local/bin/perl5.8.0 -w 
#
#   remap_gff_between_releases.pl                 
# 
# by Gary Williams                         
#
# This remaps the clone positions and genomic
# sequence of two relreases
#
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2013-06-12 16:01:34 $      

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

use Modules::Remap_Sequence_Change;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase,$species);
my ($release1, $release2, $version, $gff, $output, $genome_diffs_dir);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
            "verbose"    => \$verbose,
            "store:s"    => \$store,
            "gff=s"      => \$gff,
            "output=s"   => \$output,
            "release1=i" => \$release1,
            "release2=i" => \$release2,
            'species=s'  => \$species,
            'genomediff=s' => \$genome_diffs_dir,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species,
			     );
}

$genome_diffs_dir = $wormbase->genome_diffs if not defined $genome_diffs_dir;

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


if (! defined $release1 || ! defined $release2) {
  die "Specify the release numbers to use\n";
}
if (! defined $gff || ! defined $output) {
  die "Specify the input and output files\n";
}

##########################
# read in the mapping data
##########################

my $assembly_mapper = Remap_Sequence_Change->new($release1, $release2, $wormbase->species, $genome_diffs_dir);

##########################
# MAIN BODY OF SCRIPT
##########################

my ($indel, $change, $start_del, $end_del);

open (OUT, "> $output") || die "Can't open $output";
open (GFF, "< $gff") || die "Can't open GFF file $gff\n";

while (my $line = <GFF>) {
  chomp $line;
  
  next if $line !~ /\S/;
  next if $line =~ /^\#/;

  my @f = split /\t+/, $line;

  my ($chromosome, $start, $end, $sense) = ($f[0], $f[3], $f[4], $f[6]);

  $chromosome =~ s/^CHROMOSOME_//;
  if ($wormbase->species eq 'elegans') {
    $chromosome = "CHROMOSOME_${chromosome}";
  }

  my ($obj_name) = ($f[8] =~ /^\S+\s+\"(\S+)\"/);
  
  # some checks for malformed GFF files
  if (not defined $chromosome or
      not defined $start or 
      not defined $end or
      not defined $sense) {
    $log->log_and_die("Malformed line (missing values? spaces instead of TABs?): $line\n");
  }

  if ($start !~ /^\d+$/) {
    $log->log_and_die("Malformed line (non-numeric start?): $line\n");
  }
  if ($end !~ /^\d+$/) {
    $log->log_and_die("Malformed line (non-numeric end?): $line\n");
  }
  if ($sense !~ /^[\+|\-]$/) {
    $log->log_and_die("Malformed line (invalid sense?): $line\n");
  }

  print "chrom, start, end=$chromosome, $start, $end\n" if ($verbose);
  ($f[3], $f[4], $f[6], $indel, $change, $start_del, $end_del) = $assembly_mapper->remap_gff($chromosome, $start, $end, $sense);
  
  if ($indel or $change or $start_del or $end_del) {
    $log->write_to("There was a complication in mapping $obj_name (change=$change indel=$indel startdel=$start_del enddel=$end_del) " . 
                   "($chromosome $start $end mapped to $chromosome $f[3] $f[4])\n");
  } 
  
  $line = join "\t", @f;
  print OUT $line,"\n";
}

close (GFF);
close (OUT);

# Close log files and exit
$log->write_to("Finished.\n");

$log->mail();
exit(0);






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
    exit (0);
  }
}

##########################################

# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - remap_gff_between_releases.pl

=head1 USAGE

=over 4

=item remap_gff_between_releases.pl [options]

=back

This script reads in a GFF file and the numbers of two releases of
wormbase and maps the chromosomal locations of the first release to
the second release.

script_template.pl MANDATORY arguments:

=over 4

=item -release1 The first (earlier) database to convert from e.g. 140

=back

=item -release2 The second (later) database to convert to e.g. 155

=back

=item -gff the name of the GFF file to read in and convert

=back

=item -outfile the name of the converted GFF file to write out

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

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

=item Gary Williams

=back

=cut
