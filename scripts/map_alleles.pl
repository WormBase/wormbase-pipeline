#!/usr/local/bin/perl5.6.0 -w                    # perl5.6.0 and -w flag
#
# map_alleles.pl                             
# 
# by Anthony Rogers                               
#
# This maps alleles to the genome based on their flanking sequence
#
# Last updated by: $Author: ar2 $                      # These lines will get filled in by cvs and helps us
# Last updated on: $Date: 2002-08-19 08:34:37 $        # quickly see when script was last changed and by whom


use strict;                                     
use lib "/wormsrv2/scripts/";                    
use Wormbase;
use Ace;

use Getopt::Std;
#######################################
# command-line options                #
#######################################

use vars qw($opt_d);
# $opt_d debug   -  all output goes to ar/allele_mapping

getopts ('d');

##############
# variables  #                                                                   
##############

# Most checking scripts should produce a log file that is a) emailed to us all 
# and b) copied to /wormsrv2/logs

my $maintainers = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
my $log;
my $fa_file;
my $psl_file;

if( defined($opt_d) ) {  
  $log        = "/wormsrv2/logs/map_alleles.$rundate";
  $psl_file = glob("~ar2/allele_mapping/alleles.psl");
  $fa_file = "/wormsrv2/autoace/BLATS/alleles.fa";
}
else {
  $fa_file = glob("~ar2/allele_mapping/alleles.fa");
  $psl_file = glob("~ar2/allele_mapping/alleles.psl");
  $log        = glob("~ar2/allele_mapping/map_alleles.$rundate");
}

open (LOG,">$log");

#get allele info from database
my $db = Ace->connect(-path  => '/wormsrv1/geneace/') || do { print LOG "Connection failure: ",Ace->error; die();};
my @alleles = $db->fetch(-query =>'Find Allele;flanking_sequences');

my %allele_data;
my $count = 0;
foreach my $allele (@alleles)
  {
    $count++;
    my $type = $allele->Allelic_difference;
    my $left = $allele->Flanking_sequences;
    my $right = $allele->Flanking_sequences->right;

    $allele_data{$allele}[0] = $type;
    $allele_data{$allele}[1] = $left;
    $allele_data{$allele}[2] = $right;

    print "$allele $type $left $right\n";
  }
print "$count alleles\n";

#make the .fa file
open (FA, ">$fa_file") or die "cant open $fa_file\n";

foreach (keys %allele_data )
  {
    print FA "\> $_\_5\n$allele_data{$_}[1]\n";
    print FA "\> $_\_3\n$allele_data{$_}[2]\n" 
  }

close FA;

$db->close;
my @fail_count;
for (my $tileSize = 18 ;$tileSize > 5; $tileSize-- )
  {
    print "Starting BLAT with tileSize $tileSize\n";
    `/nfs/disk100/wormpub/blat/blat -minMatch=1 -tileSize=$tileSize -dots=10 /wormsrv2/autoace/BLAT/autoace.fa $fa_file $psl_file` or die "BLAT failed\n";
    
    #read in the results
    open (BLATS, "<$psl_file") or die "cant open $psl_file\n";
    my @data;
    while (<BLATS>)
      {
	@data = split(/\s+/,$_);
	if( $data[9] =~ /(\w+)_([3,5])/ )
	  {
	    $allele_data{$1}[3] = $data[13];     #superlink
	    if( "$2" eq "3" ) {
	      $allele_data{$1}[5] = $data[15];     #start of 3'
	    }
	else{
	  $allele_data{$1}[4] = $data[16];     #end of 5'
	}
	  }
      }
    close BLATS;
    my $fails = glob("~ar2/allele_mapping/failed_to_map");
    open (FAIL,">$fails");
    foreach (keys %allele_data)
      {
	if( defined($allele_data{$_}[4]) and defined($allele_data{$_}[5]) ) {
	  print "allele $_ deletes from $allele_data{$_}[4] to $allele_data{$_}[5] of superlink $allele_data{$_}[3]\n";
	}
    else{ 
      print FAIL "allele $_ failed mapping";
      if( defined($allele_data{$_}[4]) ){
	print FAIL "5' = $allele_data{$_}[4]";
      }
      if( defined($allele_data{$_}[5]) ) {
	print FAIL ":  3' = $allele_data{$_}[5]";
      }
      print FAIL "\n\n";
    }
      }
    close FAIL;
    $fail_count[$tileSize] = `grep -c fail $fails`;
    print "tileSize $tileSize gives $fail_count[$tileSize] failures\n";
  }    

foreach (@fail_count)
  {
    print "tileSize $_ gives $fail_count[$_] failures\n";
  }



close LOG;



















# Always good to cleanly exit from your script
exit(0);


# Add perl documentation in POD format
# This should expand on your brief description above and add details of any options
# that can be used with the program.  Such documentation can be viewed using the perldoc
# command.


__END__

=pod

=head2 NAME - map_alleles.pl

=head1 USAGE

=over 4

=item map_alleles.pl  [-options]

=back

This script:


map_alleles.pl MANDATORY arguments:

=over 4

=item none

=back

map_alleles.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers ( ar2@sanger.ac.uk)

=back

=cut
