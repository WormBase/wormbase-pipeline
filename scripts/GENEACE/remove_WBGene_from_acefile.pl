#!/software/bin/perl -w
#
# remove_WBGene_from_acefile.pl
# 
# by Gary Williams
#
# ACE editing script for MA
#
# Last updated by: $Author: mt3 $     
# Last updated on: $Date: 2007-10-15 13:34:57 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Getopt::Long;
use Carp;

######################################
# variables and command-line options # 
######################################

my ($acefile, $remove_list, $out);

GetOptions ("acefile:s"       => \$acefile,
            "remove_list:s"   => \$remove_list,
	    "out:s"           => \$out,
	    );


# read remove list
my @remove = &read_remove($remove_list);

# read/write ace file
my $remove = 0;
open (OUT, ">$out") || die "Can't open $out";
open (A, "<$acefile") || die "Can't open $acefile";
while (my $line = <A>) {
  if ( my ($id) = $line =~ /Gene\s+:\s+(\S+)/) {
    $remove = 0;
    print "ID = $id\n";
    if (grep /$id/, @remove) {
      print "Remove $id\n";
      $remove = 1;
    } else {
      print OUT "$line";
    }
    while (my $line = <A>) {
      if ($line =~ /^\s*$/) {
	print OUT "\n";
	last;
      } elsif (!$remove) {
	print OUT "$line";
      }
    }
  }
}
close(A);
close (OUT);


print "Finished.\n";
exit(0);






##############################################################
#
# Subroutines
#
##############################################################


sub read_remove {
  my ($file) = @_;
  my @list;
  open (R, "<$file") || die "\n";
  while (my $line = <R>) {
#Gene : "WBGene00023575"
    if ($line =~ /^Gene\s+:\s+\"(\S+)\"/) {
      push @list, $1;
      print "$1\n";
    }
  }
  close(R);
  return @list;
}


__END__

=pod

=head2 NAME - remove_WBGene_from_acefile.pl

=head1 USAGE

=over 4

=item remove_WBGene_from_acefile.pl  [-options]

=back

This script removes WBGene objects (as listed in remove.list) from an ace file of WBGene objects.

remove_WBGene_from_acefile.pl MANDATORY arguments:

=over 4

=item -acefile

=item -remove_list

=item -out

=back

remove_WBGene_from_acefile.pl  OPTIONAL arguments:

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

=item Keith Bradnam (krb@sanger.ac.uk)

=back

=cut
