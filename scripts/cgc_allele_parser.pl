#!/usr/local/bin/perl5.6.1 -w
# 
# cgc_allele_parser.pl
#
# by Keith Bradnam
#
# Script to convert cgc allele/lab links into ace file for geneace
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2002-12-17 12:40:56 $

use strict;
use lib "/wormsrv2/scripts/";                    
use Wormbase;
use Getopt::Long;
use Ace;

##############################
# command-line options       #
##############################
our ($help,$debug,$input,$log);
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


&create_log_files;



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


close(LOG);
exit(0);

##############################################################
#
# Subroutines
#
##############################################################

sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate     = `date +%y%m%d`; chomp $rundate;
  $log        = "/wormsrv2/logs/$script_name.$rundate.$$";

  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name\n";
  print LOG "started at ",`date`,"\n";
  print LOG "=============================================\n";
  print LOG "\n";

}

##########################################


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

This script will convert the CGC file of Laboratory/Allele information into ace format.
It should be run against the file available at

http://www.cbs.umn.edu/CGC/Nomenclature/allele.htm

The script checks geneace for extra information and writes output to standard out.

=over 4

=item MANDATORY arguments:

-input Specify input file (CGC Lab/allele info)

=back

=over 4

=item OPTIONAL arguments:

-help (this help page)

-debug <user> (debug mode)

=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk

=cut
