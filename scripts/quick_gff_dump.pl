#!/software/bin/perl -w
#
# quick_dump_gff.pl
# 
# A script to quickly dump features from a given sequence, allowing the user to specify a method if required
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2014-03-04 15:19:31 $
#====================

use strict;
use Getopt::Long;
use lib $ENV{'CVS_DIR'};

my ($database, $method, $file, $sequence, $version, $new);

GetOptions (
	    'database=s' => \$database,
            'method=s'   => \$method,
	    'file=s'     => \$file,
	    'sequence=s' => \$sequence,
	    'version=s'  => \$version,
	    'new'        => \$new
	   );

unless ($database) {die "No database defined\n";}
unless ($sequence) {die "No sequence defined\n";}
unless ($version) {print "Defaulting to gff2";$version = "2";}
unless ($file) {$file = "$sequence.gff$version";
		print "$file being used as output\n";
	      }
my $giface;
if ($new) {$giface = "/software/worm/acedb/current/bin/giface";} 
else {$giface = "/software/worm/acedb/old_versions/giface";}

my $cmd;
if ($method) {
  $cmd = sprintf("gif seqget $sequence +method $method; seqfeatures -version $version -file $file");
}
else {
  $cmd = sprintf("gif seqget $sequence; seqfeatures -version $version -file $file");
}

open (DB,"| $giface $database") or die ("failed to open giface connection to $database\n");
print DB $cmd;
close(DB) or die("Could not successfully close command '$cmd'\n");
print "Dumped gff$version data for $sequence to $file\n";
exit(0);

__END__

=pod

=head1 NAME:quick_dump.gff

=head1 USAGE:

=over 4
 
=item perl quick_gff_dump.pl [-options]

=back

quick_gff_dump.pl is a wrapper with options to quickly dump features from a given sequence, allowing the user to specify a method if required

=head1 Mandatory arguments:

=over 4

=item -database

The script needs to know the location of a database containing the speciefied sequence

=back

=over 4

=item -sequence

The clone/chromosome/scaffnew you want data for.

=back

=head1 Optional arguments:

=over 4

=item -version

standard gff2 or gff3

=back

=over 4

=item -file

Specify the patch to an output file....defaults to <sequence>.gff<2/3>

=back

=over 4

=item -method

Specify a method to dump, if none is specified you get all gff features from the specified sequence

=back

=over 4

=item -new

This flips between the new version of giface (current) and the old 4.9.39 giface binary
Always defaults to old as this doesn't require the ace2so data class which might not be in the
database you are using.

=back

=head1 EXAMPLES:

=over 4

=item quick_gff_dump.pl -database ~/wormpub/ovolvulus_curation -sequence OVOC_MITOCHONDRIAL -method curated -version 3 -file my_mito_output.gff3

Dumps the curated CDS features from the Mitochondrial genome of O. volvulus

=back

=head1 AUTHOR - Paul Davis

Email pad@ebi.ac.uk

=cut
