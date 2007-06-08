#!/usr/local/bin/perl5.8.0 -w
#
# make_Interpro2GO_mapping.pl 
# 
# by Anthony Rogers
#
# Gets latest Interpro:GO mappings from XXXX and puts info in to ace file
#
# Last updated by: $Author: mh6 $                      
# Last updated on: $Date: 2007-06-08 12:42:30 $         

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

my $help;
my $load;   # option for loading resulting acefile into autoace
my ( $debug, $test, $store );

GetOptions ("help"      => \$help,
            "load"      => \$load,
	    "test"      => \$test,
	    "debug:s"   => \$debug,
	    "store:s"   => \$store
);


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

my $rundate     = $wormbase->rundate;
my $runtime     = $wormbase->runtime;

my $log = Log_files->make_build_log($wormbase);

#Get the latest version
my $motifs = "/tmp/interpro_motifs.html";
$log->write_to(": Running 'wget' to get interpro file from EBI\n");
`/usr/bin/wget -O $motifs ftp://ftp.ebi.ac.uk/pub/databases/interpro/entry.list` and $log->log_and_die("$0 Couldnt get listing.html\n");

my $ip_ace_file = $wormbase->acefiles."/interpro_motifs.ace";
open (I2G,"<$motifs") or $log->log_and_die("cant open $motifs\n");
open (IPDESC,">$ip_ace_file") or $log->log_and_die("cant write $ip_ace_file\n");
my %interpro_des;   #IPR000018 => "P2Y4 purinoceptor"
my $text;
my $ip;

# IPR000177 Apple domain
$log->write_to(": Writing acefile to $ip_ace_file");
while (<I2G>){
  chomp;
  my @info = split;
  $ip = shift @info;
  next if( (!defined $ip) or (/entries/) ); # header lines
  #Motif : "INTERPRO:IPR000006"
  #Title    "Class I metallothionein"
  #Database         "INTERPRO" "INTERPRO_ID" "IPR000006"
  print IPDESC "Motif : \"INTERPRO:$ip\"\n";
  print IPDESC "Title  \"@info\"\n";
  print IPDESC "Database  \"INTERPRO\" \"INTERPRO_ID\" \"$ip\"\n";
  print IPDESC "\n";
}


close IPDESC;
close I2G;

# load file if -load was specified
$wormbase->load_to_database( $wormbase->autoace,
			     "$ip_ace_file",
			     'interpo_motifs') if($load);

# mail Log file
$log->mail;
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

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
