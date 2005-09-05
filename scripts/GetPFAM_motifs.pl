#!/usr/local/bin/perl5.8.0 -w
#
# GetPFAM_motifs.pl 
# 
# by Anthony Rogers
#
# Gets latest PFAM motifs from sanger/pub and puts info in to ace file
#
# Last updated by: $Author: ar2 $                      
# Last updated on: $Date: 2005-09-05 13:50:40 $         


use strict;                                     
use lib "/wormsrv2/scripts/";                    
use Wormbase;
use Getopt::Long;


######################################
# variables and command-line options # 
######################################

my ($help, $debug);
my $load;      # option for loading resulting acefile to autoace
my $maintainers = "All";
my $rundate     = &rundate;
my $runtime     = &runtime;

GetOptions ("help"      => \$help,
            "debug=s"   => \$debug, 
	    "load"      => \$load
            );

# help 
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

# create log
my $log = Log_files->make_build_log();

#Get the latest version
my $pfam_motifs_gz = "/wormsrv2/tmp/Pfam_motifs.gz";
$log->write_to("Attempting to wget the latest version\n");
print "Attempting to wget the latest version\n";
`wget -O $pfam_motifs_gz ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz` and die "$0 Couldnt get Pfam-A.full.gz \n";

`gunzip -f $pfam_motifs_gz` and die "gunzip failed\n";

my $pfam_motifs = "/wormsrv2/tmp/Pfam_motifs";
$log->write_to("Opening file $pfam_motifs\n");
print "\n\nOpening file $pfam_motifs . . \n";
open (PFAM,"<$pfam_motifs") or die "cant open $pfam_motifs\n";


my $acefile = "/wormsrv2/autoace/acefiles/pfam_motifs.ace";

open (PFAMOUT,">$acefile") or die "cant write to /wormsrv2/autoace/acefiles/pfam_motifs.ace\n";

my $text;
my $pfam;

print "\treading data . . . \n";
my $pfcount = 0;
while (<PFAM>){
  chomp;
  if ($_ =~ /^\/\//){
    if (defined $pfam){
      $pfcount++;
      print "$pfam went fine\n";
      print PFAMOUT "Motif : \"PFAM:$pfam\"\n";
      print PFAMOUT "Title \"$text\"\n";
      print PFAMOUT "Database \"Pfam\" \"Pfam_ID\" \"$pfam\"\n";
      print PFAMOUT "\n";
      undef $pfam;
      $text = "";
    }
    else{
      print "gone through a record with out picking up pfam\n";
      die;
    }
  }
  #get the id
  if($_ =~ m/^\#=GF AC\s+(PF\d{5})/  ){ 
    $pfam = $1;
  }
  
  #get the description
  if($_ =~ m/^\#=GF DE\s+(.*$)/  ) {
    $text = $1;
    $text =~ s/\"//g;
  }         
}
  
$log->write_to("added $pfcount PFAM motifs\n");

print "finished at ",`date`,"\n";
close PFAM;
close PFAMOUT;

# load file to autoace if -load specified
if($load){
  my $command = "autoace_minder.pl -load /wormsrv2/autoace/acefiles/pfam_motifs.ace -tsuser pfam_motifs";
  my $status = system($command);
  if(($status >>8) != 0){
    die "ERROR: Loading pfam_motifs.ace file failed \$\? = $status\n";
  }
}


# tidy up and die
$log->mail("$maintainers");
exit(0);


##############################################################
#
# Subroutines
#
################################################################


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

Database" "Pfam" "Pfam_ID" "PF00351"


writes to /wormsrv2/wormbase/misc_dynamic/misc_pfam_motifs.ace

=head4 OPTIONAL arguments:

=over 4
  
=item -load
 
if specified will load resulting acefile to autoace
 
=back
 

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
