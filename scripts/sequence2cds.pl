#!/usr/local/bin/perl5.8.0 -w
#
# sequence2cds.pl
# 
# by Keith Bradnam                         
#
# A script to take ?Sequence objects and make ?CDS objects
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2003-11-17 14:46:23 $     

use strict;
use Getopt::Long;


##################
# variables etc. #
##################

my $file; # where input file is located

# pattern to match timestamps
my $ts = "-O \"(\\d{4}\-\\d{2}\-\\d{2}_\\d{2}:\\d{2}:\\d{2}\.?\\d?_\\S+|original)\"";


GetOptions ("file=s"   => \$file);


# Open output/input streams
open(IN,"<$file") || die "Can't open input file\n";
open(SEQ,">$file.sequence") || die "Couldn't open output file\n";
open(CDS,">$file.cds") || die "Couldn't open output file\n";


# reset input line separator
$/ = "";


# process input file, changing lines as appropriate and redirecting output
# to two files

while(<IN>){

  # Misc tidy up of bits of models
  s/From\s+$ts\s+Source_Exons/Source_exons/g;
  s/From_Laboratory/From_laboratory/g;
  s/From_Database/From_database/g;
  s/From_Author/From_author/g;


  # convert things which have CDS tag to ?CDS objects
  if (m/Properties\s+$ts\s+Coding\s+$ts\s+CDS/){
    # Convert to new CDS class
    s/^Sequence :/CDS :/;

    # Need to add tags for SMap
    s/Structure\s+$ts\s+From\s+$ts\s+Source/Sequence/;

    # Get rid of this line (now removed in camace)
    s/Properties\s+$ts\s+Status\s+$ts\s+Annotated\s+$ts\s+\d{4}-\d{2}-\d{2}\s+$ts\s//g;

    # Change Has_allele tag
    s/Has_allele/Allele/;
    
    # output to ?CDS file
    print CDS;
  }
  # Make changes in Parent sequence objects that might link to CDS objects
  else{
    # change Subsequence tag for parent clones, but only where used for subsequence objects
    # i.e. objects that have a '.' in their object name
    s/Structure\s+$ts\s+Subsequence\s+$ts\s+(\"[\w\d]+\.[\d\w:]+\")\s+$ts\s+(\d+)\s+$ts\s+(\d+)\s+$ts/CDS $3 $5 $7/g;

    s/Visible\s+$ts\s+Matching_Genomic/Matching_CDS/g;
    
    # Transcript_child tag now just called Transcript
    s/Transcript_child/Transcript/;

    # print to ?Sequence file
    print SEQ;
  }

}

close(IN);
close(SEQ);
close(CDS);


# reset input line separator
$/ = "\n";


exit(0);

__END__

=pod

=head2   NAME - sequence2cds.pl

=head1 USAGE

=over 4

=item sequence2cds.pl --file <file>

=back

=head1 DESCRIPTION

This script will process a dumped ace file of ?Sequence objects (containing
timestamps) and produce two output files, one for the ?Sequence class and
one for the new ?CDS class.

The script tidies up some tags (e.g. Sourc_Exons => Source_exons) and changes
other tags to be compliant with the new ?CDS model.

=back

=head1 MANDATORY arguments: --file

=over 4

=item --file <name of valid acefile>

File must be a ?Sequence class dump with timestamps.  Output files will use this
name and append .sequence and .cds for the two new files.

=back


=head1 AUTHOR Keith Bradnam (krb@sanger.ac.uk) 

=back

=cut
