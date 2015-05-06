#!/usr/bin/env perl
#
# GetInterPro_motifs.pl
# 
# Gets latest Interpro:GO mappings from XXXX and puts info in to ace file
#
# Last updated by: $Author: klh $                      
# Last updated on: $Date: 2015-05-06 12:03:48 $         

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

my $wb_go_terms = &get_wormbase_go_terms();

# the interpro.xml file can be obtained from:
# ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz
my $file = "/tmp/interpro.xml.$$.gz".$wormbase->species;

# get the interpro file from the EBI
unlink $file if -e $file;
&get_interpro($file);


my $ip_ace_file = $wormbase->acefiles."/interpro_motifs.ace";
$log->write_to("Writing acefile to $ip_ace_file\n");
open (my $out_fh,">$ip_ace_file") or $log->log_and_die("cant write $ip_ace_file\n");

open (my $in_xml, "gunzip -c $file |") or $log->log_and_die("Failed to open file $file\n");
 
#Motif : "INTERPRO:IPR001624"
#Title  "Peptidase aspartic, active site"
#Database "INTERPRO" "INTERPRO_ID" "IPR001624"
#Database "INTERPRO" "short_name" "fluffinase 3"
#GO_term "GO:0003774"
#GO_term "GO:0005198"
#GO_term "GO:0001539"
#GO_term "GO:0009288"

while (<$in_xml>) {
  if (/<interpro id=\"(\S+)\"/) { # IPid - InterPro identifier
    my $IPid = $1;

    print $out_fh "\nMotif : \"INTERPRO:$IPid\"\n";
    print $out_fh "Database \"INTERPRO\" \"INTERPRO_ID\" \"$IPid\"\n";
    if (/short_name=\"(\S+)\"/) {
      print $out_fh "Database \"INTERPRO\" \"short_name\" \"$1\"\n";
    }
    
  } elsif (/\<name\>(.+)\<\/name\>/) { # IPdesc - short description
    my $IPdesc = $1;
    $IPdesc =~ s/\t/ /g;
    
    print $out_fh "Title  \"$IPdesc\"\n";
    
  } elsif (/\<classification id=\"(\S+)\" class_type=\"GO\"\>/) { # GOterm
    my $GOterm = $1;
    if (exists $wb_go_terms->{$GOterm}) {
      print $out_fh "GO_term \"$GOterm\"\n";
    }
  }
}
close($in_xml);
close($out_fh);
unlink $file;

unless ($noload) {
  $wormbase->load_to_database( $wormbase->autoace, "$ip_ace_file",'interpo_motifs', $log);
}

$log->mail;
exit(0);


#########################
sub get_interpro {
  my ($lfile) = @_;
 
  #Get the latest version
  $log->write_to("Attempting to get latest version of Interpro XML file...\n");
  $wormbase->run_command("wget -q -O $lfile ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz");
}

#########################
sub get_wormbase_go_terms {
  my %go_terms;

  my $db = Ace->connect(-path => $wormbase->autoace ) or $log->log_and_die("Could not connect to autoace\n");
  my $q = $db->fetch_many(-query => "FIND GO_term Version");
  while(my $obj = $q->next) {
    $go_terms{$obj->name} = 1;
  }
  $db->close();

  return \%go_terms;
}

__END__

