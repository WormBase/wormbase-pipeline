#!/usr/local/bin/perl -w
#
# check_protein_ID.pl
# v0.1
# dl
# 2000-07-12
#
# Consistency script to check camace for 'Protein_id' tags.
#  (a) warns if Predicted_gene does not have a Protein_id tag
#  (b) warns if Protein_id tag is corrupted
#
# Usage : check_protein_ID.pl 
#
# History:
#
# v0.1
# 000712 : dan : PP version
# 020116 : dan : moved to /wormsrv1/camace as acedb database

use strict;

#####################################################################################################
# Declare variables                                                                                 # 
#####################################################################################################

my $CDS = "";                                           # CDS name
my $cds_count = 0;                                      # No. of CDS in database
my $cds_pid_count = 0;                                  # No. of CDS with Protein_ID in database
my $PID_count = 0;                                      # No. of Protein_ID tags in database
my $errors = 0;                                         # No. of CDS with corrupted Protein_ID tags
my @LOG = ();                                           # List of CDS with corrupted Protein_ID tags
my $tace = glob("~acedb/RELEASE.SUPPORTED/bin.ALPHA_4/tace"); # tace executable path
my $db = "/wormsrv1/camace";                            # Database path
my $exec="$tace $db";     
my $command1=<<EOF;
find Predicted_gene
list -a
quit
EOF
my $command2=<<EOF;
find Predicted_gene
show -a Protein_id
quit
EOF

#####################################################################################################
# MAIN LOOP                                                                                         # 
#####################################################################################################

#####################################################################################################
# Count the No. of CDS objects in the database                                                      # 
#####################################################################################################

open (output1, ">/tmp/gene_list");
open(textace1, "echo '$command1' | $exec -| ");
while (<textace1>) {
    next if ($_ eq "");
    next if (/acedb/);
    next if (/\/\//);
    if (/^Sequence : \"(\S+)\"/) {
       print output1 "$1\n";
	$cds_count++;
    }
}    
close (textace1);     
close (output1);

#####################################################################################################
# Count the No. of CDS objects with Protein_ID tags                                                 # 
#####################################################################################################

open (output2, ">/tmp/protein_ID_list");
open(textace2, "echo '$command2' | $exec -| ");
while (<textace2>) {
    chomp;
    next if ($_ eq "");
    next if (/acedb/);
    next if (/\/\//);
    s/\"//g;
    if (/^Sequence : (\S+)/) {
	$CDS = $1;
	print output2 "$CDS\n";
	$cds_pid_count++;
	next;
    }
    if (/Protein_id\s+(\S+)\s+(\S+)\s+(\d+)/) {
	$PID_count++;
    }
    else {
	push (@LOG, "Problem with $CDS : [$_]");
    }
}    
close textace2;     
close (output2);

#####################################################################################################
# Report numbers and list of errors / discrepencies                                                 # 
#####################################################################################################

$errors = scalar (@LOG);

print "No. of CDS objects          : '$cds_count'\n";
print "No. of CDS with Protein_IDs : '$cds_pid_count'\n";
print "No. of Protein_IDs          : '$PID_count'\n";
print "No. of CDS with problems    : '$errors'\n\n";

if ($errors > 0) {
    print "Putative errors in the database\n\n";
    foreach (@LOG) {
	next if ($_ eq "");
	print "$_\n";
    }
    print "\n";
}

if ($cds_count != $cds_pid_count) {
    print "Mismatch between No, of CDS objects and No. of objects with Protein_ID's\n\n";
    print "CDS Objects without Protein_ID's are :\n";
    system ("comm -3 /tmp/gene_list /tmp/protein_ID_list");
}

#####################################################################################################
# Tidy up                                                                                           # 
#####################################################################################################

unlink '/tmp/gene_list';
unlink '/tmp/protein_ID_list';

exit(0);

__END__

=pod

=head1 NAME - check_protein_ID.pl

=head2 AUTHOR

Dan Lawson dl1@sanger.ac.uk

=head2 USAGE

check_protein_ID.pl 

check_protein_ID.pl uses tace to query an ACeDB database for predicted 
gene objects (default database is camace). It will then check the 
'Protein_id' field for the presence/absence of the tag and perform a
simple consistency check (3 tag format - '?Sequence Text Int').

The script will report objects with corrupted tags and in the 
case of a discrepency between the number of prediced CDS and the
number with Protein_id tags (i.e. missing data), those CDS which
do not yet have a Protein_ID in the ACeDB database.

=head2 MANDATORY ARGUMENTS

=over 2

=item *
none

=back
=cut

=head2 REVISIONS

=over 2

=item *
v0.1
2000-07-12 : PP version

=back
=cut



