#!/usr/local/bin/perl5.8.0 -w
#
# make_Interpro2GO_mapping.pl 
# 
# by Anthony Rogers
#
# Gets latest Interpro:GO mappings from XXXX and puts info in to ace file
#
# Last updated by: $Author: mh6 $                      
# Last updated on: $Date: 2014-02-27 10:28:18 $         

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

my ( $debug, $test, $store,$help,$noload );

GetOptions ("help"      => \$help,
            "noload"    => \$noload,
	    "test"      => \$test,
	    "debug:s"   => \$debug,
	    "store:s"   => \$store
)||die($@);

# Display help if required
&usage("Help") if ($help);

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

my $log = Log_files->make_build_log($wormbase);

# the interpro.xml file can be obtained from:
# ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz

# store it here
my $dir = "/tmp";
my $file = "$dir/interpro.xml_".$wormbase->species;

# get the interpro file from the EBI
unlink $file if -e $file;
get_interpro($file);
 
my $ip_ace_file = $wormbase->acefiles."/interpro_motifs.ace";
open (IPDESC,">$ip_ace_file") or $log->log_and_die("cant write $ip_ace_file\n");

# IPR000177 Apple domain
$log->write_to(": Writing acefile to $ip_ace_file\n");
open (XML, "< $file") || die "Failed to open file $file\n";
 
#Motif : "INTERPRO:IPR001624"
#Title  "Peptidase aspartic, active site"
#Database "INTERPRO" "INTERPRO_ID" "IPR001624"
#Database "INTERPRO" "short_name" "fluffinase 3"
#GO_term "GO:0003774"
#GO_term "GO:0005198"
#GO_term "GO:0001539"
#GO_term "GO:0009288"

while (my $line = <XML>) {
    if ($line =~ /<interpro id=\"(\S+)\"/) { # IPid - InterPro identifier
      my $IPid = $1;

      print IPDESC "\nMotif : \"INTERPRO:$IPid\"\n";
      print IPDESC "Database \"INTERPRO\" \"INTERPRO_ID\" \"$IPid\"\n";
      print IPDESC "Database \"INTERPRO\" \"short_name\" \"$1\"\n" if $line=~/short_name=\"(\S+)\"/;

    } elsif ($line =~ /\<name\>(.+)\<\/name\>/) { # IPdesc - short description
      my $IPdesc = $1;
      $IPdesc =~ s/\t/ /g;

      print IPDESC "Title  \"$IPdesc\"\n";

    } elsif ($line =~ /\<classification id=\"(\S+)\" class_type=\"GO\"\>/) { # GOterm
      my $GOterm = $1;

      print IPDESC "GO_term \"$GOterm\"\n";

    }
}
close XML;
close IPDESC;
unlink $file;

unless ($noload) {
    $wormbase->load_to_database( $wormbase->autoace, "$ip_ace_file",
                                 'interpo_motifs', $log);
}

# mail Log file
$log->mail;
exit(0);

###############################
# show the perldoc

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

#########################
# get the interpro file #

sub get_interpro {

  my $latest_version = $_[0];
 
  #Get the latest version
  print "Attempting to FTP the latest version of interpro.xml from ebi \n";
  `wget -q -O $latest_version.gz ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz`;
  `gunzip "${latest_version}.gz"`;
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

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
