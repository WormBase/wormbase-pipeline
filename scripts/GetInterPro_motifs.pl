#!/usr/local/bin/perl5.6.0 -w                   
#
# make_Interpro2GO_mapping.pl 
# 
# by Anthony Rogers
#
# Gets latest Interpro:GO mappings from XXXX and puts info in to ace file
#
# Last updated by: $Author: ar2 $                      # These lines will get filled in by cvs and helps us
# Last updated on: $Date: 2002-09-16 14:51:49 $                        # quickly see when script was last changed and by whom


use strict;                                     
use lib "/wormsrv2/scripts/";                    
use Wormbase;

# Try to keep different parts of code cleanly separated using comments...

##############
# variables  #                                                                   #
##############

# Most checking scripts should produce a log file that is a) emailed to us all 
# and b) copied to /wormsrv2/logs

my $maintainers = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
our $log        = "/wormsrv2/logs/GetInterPro_motifs.$rundate";

open (LOG, ">$log") or die "cant open $log";

print LOG "$0\n";
print LOG "started at ",`date`,"\n";
print LOG "=============================================\n";
print LOG "\n";

#Get the latest version
my $motifs = "/wormsrv2/tmp/interpro_motifs.html";
print LOG "Attempting to wget the latest version\n";
`wget -O $motifs http://www.ebi.ac.uk/interpro/listing.html` and die "$0 Couldnt get listing.html\n";
print LOG "...... got it!\n";

print LOG "\n\nOpening file $motifs . . \n";
open (I2G,"<$motifs") or die "cant open $motifs\n";
open (IPDESC,">/wormsrv2/wormbase/misc/misc_interpro_motifs.ace") or die "cant write misc_interpro_motifs.ace\n";
my %interpro_des;   #IPR000018 => "P2Y4 purinoceptor"
my $text;
my $ip;

print LOG "\treading data . . . \n";
# <a href="http://www.ebi.ac.uk/interpro/IEntry?ac=IPR000159">IPR000159</a> RA domain<br>
while (<I2G>)
  {
    if( $_ =~ m/IEntry\?ac=(IPR\d{6})/ ) 
      {
	$ip = $1;
	if( $_ =~ m/a>\s(.*)<br>?/ ) #for some reason they all start with a space
	  {
	    $text = $1;
	    #Motif : "INTERPRO:IPR000006"
	    #Title    "Class I metallothionein"
	    #Database         "INTERPRO" "IPR000006" "IPR000006"
	    print IPDESC "Motif : \"INTERPRO:$ip\"\n";
	    print IPDESC "Title  \"$text\"\n";
	    print IPDESC "Database  \"INTERPRO\" \"$ip\"\n";
	    print IPDESC "\n";
	  }
	else {
	  print LOG "Cant find a description for $ip ($_)\n";}
      }
  }

############################################################
print LOG " adding patch for IPR005560 (in the first version, this one was missing the <br> tag in the HTML file so was missed)\n if its no longer picked out above this can be removed from script\n\n";
# in the first version, this one was missing the <br> tag in the HTML file so was missed
print IPDESC "Motif : \"INTERPRO:IPR005560\"\n";
print IPDESC "Title  \"Domain of Unknown Function DUF326\"\n";
print IPDESC "Database  \"INTERPRO\" \"IPR005560\"\n";
############################################################


print LOG "finsihed at ",`date`,"\n";
close IPDESC;
close I2G;
close LOG;
#### use Wormbase.pl to mail Log ###########
my $name = "GetInterPro_motifs";
#$maintainers = "ar2\@sanger.ac.uk";
&mail_maintainer ($name,$maintainers,$log);
#########################################

exit(0);



__END__

=pod

=head2 NAME make_InterproGO_mapping.pl

=head1 USAGE

=over 4

=item make_InterproGO_mapping.p  [-options]

=back

This script:

wgets the latest version of http://www.ebi.ac.uk/interpro/listing.html
then parses it to produce an ace file of format


Motif : "INTERPRO:IPR000018"

Title    "P2Y4 purinoceptor"

Database         "INTERPRO" "IPR000018" "IPR000018"


=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
