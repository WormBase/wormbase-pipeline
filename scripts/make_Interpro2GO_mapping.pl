#!/usr/local/bin/perl5.8.0 -w
#
# make_Interpro2GO_mapping.pl 
# 
# by Anthony Rogers
#
# Gets latest Interpro:GO mappings from XXXX and puts info in to ace file
#
# Last updated by: $Author: krb $                      # These lines will get filled in by cvs and helps us
# Last updated on: $Date: 2003-12-01 11:54:26 $                        # quickly see when script was last changed and by whom


use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Data::Dumper;

# Try to keep different parts of code cleanly separated using comments...

##############
# variables  #                                                                   #
##############

# Most checking scripts should produce a log file that is a) emailed to us all 
# and b) copied to /wormsrv2/logs

my $maintainers = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
our $log        = "/wormsrv2/logs/make_Interpro2GO_mapping.$rundate";
my $latest_version = "/wormsrv2/tmp/newip2gos";

open (LOG, ">$log") or die "cant open $log";

print LOG "$0\n";
print LOG "started at ",`date`,"\n";
print LOG "=============================================\n";
print LOG "\n";

my $get_latest = 1;
if( $get_latest == 1)
  {
    #Get the latest version
    print LOG "Attempting to FTP the latest version from ebi \n";
    `wget -O $latest_version ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro2go`;
   
  }
else {
  print LOG "Using the existing version of interpro2go mapping file (ie not FTPing latest)\n";
}

print LOG "\n\nOpening file $latest_version . . \n";
open (I2G,"<$latest_version") or die "cant open $latest_version\n";

my %interpro_des;   #IPR000018 => "P2Y4 purinoceptor"
my %interpro_GO;    #GO:0004930 GO:0016020
my $ip;
my @data;
my $description;
my $i;


print LOG "\treading data . . . \n";
while (<I2G>)
  {
    @data = split(/\s+/,$_);
    $description = "";
    if( $data[0] =~ m/:(IPR\d{6})/ )
      {
	$ip = $1;
      	unless( defined($interpro_des{$ip}) )
	  {
	    $i = 1;
	    while ($data[$i] ne ">")
	      {
		$description .= "$data[$i] ";
		$i++;
	      }
	    $interpro_des{$ip} = $description;
	  }
	if( $_ =~ m/(GO:\d{7})/g )
	  {
	    $interpro_GO{$ip} .= "$1 ";
	  }
      }
  }
close I2G;
print LOG "\tabout to write ace file  .  .  \n";

#now write .ace file

#Motif : "INTERPRO:IPR000018"
#Title    "P2Y4 purinoceptor"
#Database "INTERPRO" "INTERPRO_ID" "IPR000018"
#GO_term  "GO:0004930"
#GO_term  "GO:0005624"

open (I2GACE, ">/wormsrv2/tmp/interpro2go.ace") or die "cant write ace file\n";
foreach my $key (keys %interpro_des)
  {
    print I2GACE "Motif : \"INTERPRO:$key\"\n";
    print I2GACE "Database \"INTERPRO\" \"INTERPRO_ID\" \"$key\"\n";
    @data = split(/\s+/,"$interpro_GO{$key}");
    foreach (@data){
      print I2GACE "GO_term \"$_\"\n";
    }
    print I2GACE "\n";
  }
close I2GACE;


# write Data::Dumper has of interpro => GO mapping
open (IP2GO,">/wormsrv2/autoace/COMMON_DATA/interpro2go.dat") or die "cant open i2g\n";
print IP2GO Data::Dumper->Dump([\%interpro_GO]);
close IP2GO;



print LOG "$0 finished at ",`date`,"\n\n";
close LOG;

#### use Wormbase.pl to mail Log ###########
my $name = "make_Interpro2GO_mapping";
$maintainers = "ar2\@sanger.ac.uk";
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

wgets the latest version of ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro2go
then parses it to produce an ace file of format


Motif : "INTERPRO:IPR000018"

Title    "P2Y4 purinoceptor"

Database "INTERPRO" "INTERPRO_ID" "IPR000018"

GO_term  "GO:0004930"

GO_term  "GO:0005624"



=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
