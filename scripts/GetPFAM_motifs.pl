#!/usr/local/bin/perl5.6.1 -w                   
#
# GetPFAM_motifs.pl 
# 
# by Anthony Rogers
#
# Gets latest PFAM motifs from sanger/pub and puts info in to ace file
#
# Last updated by: $Author: ar2 $                      # These lines will get filled in by cvs and helps us
# Last updated on: $Date: 2002-09-06 12:30:17 $                        # quickly see when script was last changed and by whom


use strict;                                     
use lib "/wormsrv2/scripts/";                    
use Wormbase;
use Getopt::Std;
#######################################
# command-line options                #
#######################################

use vars qw($opt_d);
# $opt_d debug   -  all output goes to ar/allele_mapping

getopts ('d');

# Try to keep different parts of code cleanly separated using comments...

##############
# variables  #                                                                   #
##############

# Most checking scripts should produce a log file that is a) emailed to us all 
# and b) copied to /wormsrv2/logs

my $maintainers = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
our $log        = "/wormsrv2/logs/GetPfam_motifs.$rundate";

open (LOG, ">$log") or die "cant open $log";

print LOG "$0\n";
print LOG "started at ",`date`,"\n";
print LOG "=============================================\n";
print LOG "\n";

unless (defined($opt_d))
  {
    #Get the latest version
    my $pfam_motifs_gz = "/wormsrv2/tmp/Pfam_motifs.gz";
    print LOG "Attempting to wget the latest version\n";
    print "Attempting to wget the latest version\n";
    `wget -O $pfam_motifs_gz ftp://ftp.sanger.ac.uk/pub/databases/Pfam/Pfam-A.full.gz` and die "$0 Couldnt get Pfam-A.full.gz \n";
    print LOG "...... got it!\nUnzipping . .";
    print "...... got it!\nUnzipping . .";
    `gunzip $pfam_motifs_gz` and die "gunzip failed\n";
    print LOG "DONE\n";
    print "DONE\n";
  }
else{
  $maintainers = "ar2\@sanger.ac.uk";
}

my $pfam_motifs = "/wormsrv2/tmp/Pfam_motifs";
print LOG "\n\nOpening file $pfam_motifs . . \n";
print "\n\nOpening file $pfam_motifs . . \n";
open (PFAM,"<$pfam_motifs") or die "cant open $pfam_motifs\n";

my $acefile;

if (defined ($opt_d)) {
    $acefile = "/wormsrv2/tmp/misc_pfam_motifs.ace";
  }
else {
  $acefile = "/wormsrv2/wormbase/misc/misc_pfam_motifs.ace";
}
open (PFAMOUT,">$acefile") or die "cant write misc_pfam_motifs.ace\n";

my $text;
my $pfam;

print LOG "\treading data . . . \n";
print "\treading data . . . \n";
while (<PFAM>)
  {
    chomp;
    if ($_ =~ /^\/\//)
      {
	if (defined $pfam)
	  {
	    print "$pfam went fine\n";
	    print PFAMOUT "Motif : \"PFAM:$pfam\"\n";
	    print PFAMOUT "Title \"$text\"\n";
	    print PFAMOUT "Database \"Pfam\" \"PFAM:$pfam\" \"$pfam\"\n";
	    print PFAMOUT "\n";
		undef $pfam;
	    $text = "";
	  }
	else{
	  print "gone thru a record with out picking up pfam\n";
	  die;
	}
      }
    #get the id
    if($_ =~ m/^\#=GF AC\s+(PF\d{5})/  )
      { $pfam = $1;}
    
    #get the description
    if($_ =~ m/^\#=GF DE\s+(.*$)/  ) 
      {
	$text = $1;
	$text =~ s/\"//g;
      }      
  }


print LOG "finsihed at ",`date`,"\n";
print "finsihed at ",`date`,"\n";
close PFAM;
close PFAMOUT;
close LOG;
#### use Wormbase.pl to mail Log ###########
my $name = "GetPFAM_motifs";
&mail_maintainer ($name,$maintainers,$log);
#########################################

exit(0);



__END__

=pod

=head2 NAME GetPFAM_motifs.pl

=head1 USAGE

=over 4

=item GetPFAM_motifs.pl

=back

This script:

wgets the latest version of ftp.sanger.ac.uk/pub/databases/Pfam/Pfam-A.full.gz
unzips it then parses it to produce an ace file of format

Motif : "PFAM:PF00351"

Title "Biopterin-dependent aromatic amino acid hydroxylase"

Database" "Pfam" "PFAM:PF00351" "PF00351"


writes to /wormsrv2/wormbase/misc/misc_pfam_motifs.ace

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
