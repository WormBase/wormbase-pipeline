#!/usr/local/bin/perl5.6.1 -w
#
# new_wormpep_entries.pl
#
# by Keith Bradnam
#
# Script to find new wormpep entries in the last release of wormpep
# and produce a fasta file of these sequence to load into the wormprot
# mysql database prior for the pre-build pipeline.
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2002-08-27 16:10:09 $


#################################################################################
# variables                                                                     #
#################################################################################

$|=1;
use IO::Handle;
use Getopt::Std;
use Cwd;
use lib '/wormsrv2/scripts';
use Wormbase;
use strict;

 ##############################
 # command-line options       #
 ##############################

my $opt_h = "";      # Help/Usage page
getopts ('h');
&usage if ($opt_h);


 ##############################
 # Script variables (run)     #
 ##############################
my $release_name = &get_wormbase_version_name;
my $release      = &get_wormbase_version;

our $dbfile = "/wormsrv2/WORMPEP/wormpep${release}/wp.fasta${release}";

# make a hash of wormpep IDs -> peptide sequence
our %wormpep;
&make_hash;

 ###############################
 # Main Loop                   #
 ###############################

open (LOG,  ">/wormsrv2/WORMPEP/wormpep${release}/new_entries.$release_name") || die "Couldn't write output file\n";;
open (DIFF, "</wormsrv2/WORMPEP/wormpep${release}/wormpep.diff${release}") || die "Couldn't read from diff file\n";
while (<DIFF>) {    
  my ($new_acc,$seq);
    if (/changed:\t\S+\s+\S+ --> (\S+)/) {
      $new_acc = $1;
      print LOG ">$new_acc\n$wormpep{$new_acc}";
    }
    if (/new:\t\S+\s+(\S+)/) {
      $new_acc = $1;
      print LOG ">$new_acc\n$wormpep{$new_acc}";

    }
}
close(DIFF);
close(LOG);
exit(0);


##################
# Usage
##################

sub usage {
    system ('perldoc',$0);
    exit;
}

#################################

sub make_hash {
  my $id;
  my $seq;
  open (WP, "<$dbfile") || die "Couldn't read from file\n";
  while (<WP>) {
    if (/^>(\S+)/) {
      chomp;
      my $new_id = $1;
      if ($id) {
	$id =~ /CE0*([1-9]\d*)/ ; 
	$wormpep{$id} = $seq;
      }
      $id = $new_id ; $seq = "" ;
    } 
    elsif (eof) {
      if ($id) {
	$id =~ /CE0*([1-9]\d*)/ ;
	$seq .= $_ ;
	$wormpep{$id} = $seq;
      }
    } 
    else {
      $seq .= $_ ;
    }
  }

  close (WP);
  
}
exit(0);


__END__

=pod

=head1 NAME - new_wormpep_entries.pl

=head2 USAGE

This file looks at the wormpep.diff file from the last release of
wormpep and extracts any new peptide sequences which it writes to
a file inside the respective wormpep directory.  This fasta file
is needed by the prebuild process to populate the wormprot mysql
database with the new protein entries.


Optional arguments:

=over 4

=item *

-h  this help

=back

=cut
