#!/usr/local/bin/perl5.8.0 -w
#
# newgene.pl
#
# by Keith Bradnam
#
# simple script for creating new (sequence based) Gene objects 
#
# Last edited by: $Author: krb $
# Last edited on: $Date: 2004-07-15 15:12:39 $

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;

my $seq;    # cds/transcript/pseudogene name (including C:, P: or T: prefix)
my $class;  # will store class (C, P, or T)
my $person; # ID of WBPerson curator
my $id;     # new gene ID
my $load;   # option to load directly to geneace after making acefile

GetOptions ("seq=s"   => \$seq,
	    "who=i"   => \$person,
            "id=i"    => \$id,
	    "load"    => \$load);

die "no -seq option\n" if (!$seq);
die "no -who option\n" if (!$person);
die "no -id option\n"  if (!$id);

$seq =~ m/(.):(.*)/;
$class = $1;
$seq = $2;

# write information to acefile

open(OUT, ">/wormsrv1/geneace/fix.ace") || die "Can't write to output file\n";

print OUT "Gene WBGene000${id}\n";
print OUT "Live\n";
print OUT "Version 1\n";
print OUT "Sequence_name $seq\n";
print OUT "Public_name $seq\n";
print OUT "Species \"Caenorhabditis elegans\"\n";
print OUT "History Version_change 1 now WBPerson${person} Event Created\n";
print OUT "CDS $seq\n"        if ($class eq "C");
print OUT "Transcript $seq\n" if ($class eq "T");
print OUT "Pseudogene $seq\n" if ($class eq "P");
print OUT "Method Gene\n";

close(OUT);

# load information to geneace if -load is specified
if ($load){
  my $tace = &tace;
  my $command = "pparse /wormsrv1/geneace/fix.ace\nsave\nquit\n";
  open (GENEACE,"| $tace -tsuser \"krb\" /wormsrv1/geneace") || die "Failed to open pipe to /wormsrv1/geneace\n";
  print GENEACE $command;
  close GENEACE;
}


exit(0);

                                                                                                                                         
__END__
                                                                                                                                         
=pod
 
=head2   NAME - newgene.pl
 
=head1 USAGE
 
=over 4
 
=item newgene.pl -[options]
 
=back
 
=head1 DESCRIPTION
 
Very simple script designed to create new gene objects to load into geneace.  Mainly written to
save time from adding all the mandatory tags that each new object needs.  Just supply
a sequence name, class, person ID of curator providing the information and a new Gene object
ID.  Resulting acefile will be made in /wormsrv1/geneace/fix.ace

E.g.

newgene.pl -seq C:AH6.24 -who 1971 -id 23428 -load

The 'C:' prefix indicates it is a CDS, use 'P' for Pseudogenes and 'T' for Transcripts

This would produce the following acefile:

Gene WBGene00023428
Live
Version 1
Sequence_name AH6.24
Public_name AH6.24
Species "Caenorhabditis elegans"
History Version_change 1 now WBPerson1971 Event Created
CDS AH6.24
Method Gene

Note that this doesn't handle any new genes based on CGC names, only those which have been created
by WashU and Sanger curators.


The -load option will attempt to load the acefile into geneace


=head1 AUTHOR Keith Bradnam (krb@sanger.ac.uk)
 
=back
 
=cut
