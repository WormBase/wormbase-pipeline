#!/usr/local/bin/perl5.6.1 -w                   
#
# GetPFAM_motifs.pl 
# 
# by Anthony Rogers
#
# Gets latest PFAM motifs from sanger/pub and puts info in to ace file
#
# Last updated by: $Author: krb $                      
# Last updated on: $Date: 2002-12-09 16:53:55 $         


use strict;                                     
use lib "/wormsrv2/scripts/";                    
use Wormbase;
use Getopt::Long;


######################################
# variables and command-line options # 
######################################

my ($help, $debug);
our $log;
my $maintainers = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;

GetOptions (
            "help"      => \$help,
            "debug=s"   => \$debug
            );

# help 
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

&create_log_files;

#Get the latest version
my $pfam_motifs_gz = "/wormsrv2/tmp/Pfam_motifs.gz";
print LOG "Attempting to wget the latest version\n";
print "Attempting to wget the latest version\n";
`wget -O $pfam_motifs_gz ftp://ftp.sanger.ac.uk/pub/databases/Pfam/Pfam-A.full.gz` and die "$0 Couldnt get Pfam-A.full.gz \n";
print LOG "...... got it!\nUnzipping . .";
print "...... got it!\nUnzipping . .";
`gunzip -f $pfam_motifs_gz` and die "gunzip failed\n";
print LOG "DONE\n";
print "DONE\n";

my $pfam_motifs = "/wormsrv2/tmp/Pfam_motifs";
print LOG "\n\nOpening file $pfam_motifs . . \n";
print "\n\nOpening file $pfam_motifs . . \n";
open (PFAM,"<$pfam_motifs") or die "cant open $pfam_motifs\n";

my $acefile = "/wormsrv2/wormbase/misc/misc_pfam_motifs.ace";

open (PFAMOUT,">$acefile") or die "cant write misc_pfam_motifs.ace\n";

my $text;
my $pfam;

print LOG "\treading data . . . \n";
print "\treading data . . . \n";
my $pfcount = 0;
while (<PFAM>)
  {
    chomp;
    if ($_ =~ /^\/\//)
      {
	if (defined $pfam)
	  {
	    $pfcount++;
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

print LOG "added $pfcount PFAM motifs\n";
print LOG "finished at ",`date`,"\n";
print "finished at ",`date`,"\n";
close PFAM;
close PFAMOUT;
close LOG;

#### use Wormbase.pl to mail Log ###########
my $name = "GetPFAM_motifs";
&mail_maintainer ($name,$maintainers,$log);

exit(0);


##############################################################
#
# Subroutines
#
################################################################

sub create_log_files{

  my $rundate     = `date +%y%m%d`; chomp $rundate;

  # touch logfile for run details
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

  $log        = "/wormsrv2/logs/$1.$rundate.$$";
  open (LOG, ">$log") or die "cant open $log";
  print LOG "$0\n";
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
