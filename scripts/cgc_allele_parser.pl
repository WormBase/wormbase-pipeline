#!/usr/local/bin/perl5.6.1 -w
# 
# cgc_allele_parser.pl
#
# by Keith Bradnam
#
# Script to convert cgc allele/lab links into ace file for geneace
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2002-12-10 10:16:04 $

use strict;
use lib "/wormsrv2/scripts/";                    
use Wormbase;
use Getopt::Long;
use Ace;

##############################
# command-line options       #
##############################
our ($help,$debug,$input);
my $maintainers = "All";

GetOptions ("help"      => \$help,
            "debug=s"   => \$debug,
	    "input=s"   => \$input);


# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}



$|=1;

# hash to store allele 2 letter code and lab 2 letter codes
my %alleles_labs;

open(INPUT, "<$input") || die "Couldn't open input file\n";
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



exit(0);


###########################################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

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

-input Specify input file (CGC strain file)

=back

=over 4

=item OPTIONAL arguments:

-help (this help page)

-debug <user> (debug mode)

=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk

=cut
