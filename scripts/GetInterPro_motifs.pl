#!/usr/local/bin/perl5.8.0 -w
#
# make_Interpro2GO_mapping.pl 
# 
# by Anthony Rogers
#
# Gets latest Interpro:GO mappings from XXXX and puts info in to ace file
#
# Last updated by: $Author: krb $                      
# Last updated on: $Date: 2004-09-01 13:50:57 $         

use strict;                                     
use lib "/wormsrv2/scripts/";                    
use Wormbase;
use Getopt::Long;

my $help;
my $load;   # option for loading resulting acefile into autoace

GetOptions ("help"      => \$help,
            "load"      => \$load);


# Display help if required
&usage("Help") if ($help);


##############
# variables  #                                                                   
##############

# Most checking scripts should produce a log file that is a) emailed to us all 
# and b) copied to /wormsrv2/logs

my $maintainers = "All";
my $rundate     = &rundate;
my $runtime     = &runtime;
our $log        = "/wormsrv2/logs/GetInterPro_motifs.$rundate";

open (LOG, ">$log") or die "cant open $log";

print LOG "$0\n";
print LOG &runtime, ": script started\n\n";

#Get the latest version
my $motifs = "/tmp/interpro_motifs.html";
print LOG &runtime, ": Running 'wget' to get interpro file from EBI\n";
`/usr/local/bin/wget -O $motifs ftp://ftp.ebi.ac.uk/pub/databases/interpro/entry.list` and die "$0 Couldnt get listing.html\n";


open (I2G,"<$motifs") or die "cant open $motifs\n";
open (IPDESC,">/wormsrv2/wormbase/acefiles/interpro_motifs.ace") or die "cant write /wormsrv2/autoace/acefiles/interpro_motifs.ace\n";
my %interpro_des;   #IPR000018 => "P2Y4 purinoceptor"
my $text;
my $ip;

# IPR000177 Apple domain

print LOG &runtime, ": Writing acefile to /wormsrv2/autoace/acefiles/interpro_motifs.ace\n";
while (<I2G>){
  chomp;
  my @info = split;
  $ip = shift @info;
  next if( (!defined $ip) or ("$info[0]" eq "entries") ); # header lines
  #Motif : "INTERPRO:IPR000006"
  #Title    "Class I metallothionein"
  #Database         "INTERPRO" "INTERPRO_ID" "IPR000006"
  print IPDESC "Motif : \"INTERPRO:$ip\"\n";
  print IPDESC "Title  \"@info\"\n";
  print IPDESC "Database  \"INTERPRO\" \"INTERPRO_ID\" \"$ip\"\n";
  print IPDESC "\n";
}


print LOG &runtime, ": script finished\n";
close IPDESC;
close I2G;
close LOG;

# load file if -load was specified
if($load){
  print LOG &runtime, ": Loading acefile to autoace\n";
  my $command = "autoace_minder.pl -load /wormsrv2/autoace/acefiles/interpro_motifs.ace -tsuser interpro_motifs";
 
  my $status = system($command);
  if(($status >>8) != 0){
    print LOG "ERROR: Loading interpro_motifs.ace file failed \$\? = $status\n";
  }
}

# mail Log file
&mail_maintainer ("GetInterPro_motifs.pl",$maintainers,$log);


exit(0);

###############################

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

Database "INTERPRO" "INTERPRO_ID" "IPR000018"


=head1 OPTIONAL arguments:

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
