#!/usr/local/bin/perl5.6.1 -w
# 
# cgc_allele_parser.pl
#
# by Keith Bradnam
#
# Script to convert cgc allele/lab links into ace file for geneace
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2002-10-17 09:57:30 $

use strict;
use Ace;
use Getopt::Std;

##############################
# command-line options       #
##############################
our $opt_i = "";      # input file
our $opt_h = "";      # help
getopts ('i:h');
&usage if ($opt_h);


$|=1;

# hash to store allele 2 letter code and lab 2 letter codes
my %alleles_labs;

open(INPUT, "<$opt_i") || die "Couldn't open input file\n";
while(<INPUT>){
  my $line = $_;
  chomp($line);
  my @fields = split (/\s+/,$line);
  $alleles_labs{$fields[1]} = $fields[0];
}
close(INPUT);

# open a local database connection
my $db = Ace->connect(-path  =>  '/wormsrv1/geneace') || die "Couldn't connect to geneace\n";


foreach my $key (sort keys %alleles_labs){
#  print "ALLELE: $key LAB: $alleles_labs{$key}\n";
  for (my $i=0;$i<10;$i++){
    my @alleles  = $db->fetch("Allele","${key}${i}*");
    foreach my $i (@alleles){
      print "Allele : \"$i\"\n";
      print "Location \"$alleles_labs{$key}\"\n\n";
    }
  }
}

$db->close;

sub usage {
  system ('perldoc',$0);
  exit;       
}




exit(0);

__END__

=pod

=head1 NAME - cgc_allele_parser.pl

=back


=head1 USAGE

=over 4

=item cgc_allele_parser.pl -i <input file>

=back

This script will convert the CGC file of strain information into ace format.
It should be run against the file available at

http://www.cbs.umn.edu/CGC/Strains/gophstrn

The script will write two ace files to your current directory, one to be loaded 
into geneace, and a second to be archived in /wormsrv1/geneace which will have 
delete instructions for removing all the data you have just added.

=over 4

=item MANDATORY arguments:

-i Specify input file (CGC strain file)

=back

=over 4

=item OPTIONAL arguments:

-h (this help page)

=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk

=cut
