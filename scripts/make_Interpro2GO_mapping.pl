#!/usr/local/bin/perl5.8.0 -w
#
# make_Interpro2GO_mapping.pl 
# 
# by Anthony Rogers
#
# Gets latest Interpro:GO mappings from XXXX and puts info in to ace file
#
# Last updated by: $Author: ar2 $                      
# Last updated on: $Date: 2006-02-17 11:32:47 $           


use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Data::Dumper;
use File::Copy;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"      => \$store,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


my $latest_version = "/tmp/newip2gos";


my $get_latest = 1;
if( $get_latest == 1)
  {
    #Get the latest version
    $log->write_to("Attempting to FTP the latest version from ebi \n");
    `wget -O $latest_version ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro2go`;
   
  }
else {
  $log->write_to("Using the existing version of interpro2go mapping file (ie not FTPing latest)\n");
}

$log->write_to("\n\nOpening file $latest_version . . \n");
open (I2G,"<$latest_version") or die "cant open $latest_version\n";

my %interpro_des;   #IPR000018 => "P2Y4 purinoceptor"
my %interpro_GO;    #GO:0004930 GO:0016020
my $ip;
my @data;
my $description;
my $i;


$log->write_to("\treading data . . . \n");
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
$log->write_to("\tabout to write ace file  .  .  \n");

#now write .ace file

#Motif : "INTERPRO:IPR000018"
#Title    "P2Y4 purinoceptor"
#Database "INTERPRO" "INTERPRO_ID" "IPR000018"
#GO_term  "GO:0004930"
#GO_term  "GO:0005624"

my $ace_dir = $wormbase->autoace;     # AUTOACE DATABASE DIR
my $common_data_dir = $wormbase->common_data; # AUTOACE COMMON_DATA

open (I2GACE, ">$ace_dir/acefiles/interpro2go.ace") or die "cant write to $ace_dir/acefiles/interpro2go.ace\n";
foreach my $key (keys %interpro_des){
  print I2GACE "Motif : \"INTERPRO:$key\"\n";
  print I2GACE "Database \"INTERPRO\" \"INTERPRO_ID\" \"$key\"\n";
  @data = split(/\s+/,"$interpro_GO{$key}");
  foreach (@data){
    print I2GACE "GO_term \"$_\"\n";
  }
  print I2GACE "\n";
}
close(I2GACE);


$log->write_to("Loading interpro2go.ace file to autoace\n");
my $command = "autoace_minder.pl -load $ace_dir/acefiles/interpro2go.ace -tsuser interpro2go_mappings";
$wormbase->run_script($command, $log);


# write Data::Dumper of interpro => GO mapping
open (IP2GO,">$common_data_dir/interpro2go.dat") or die "cant open i2g\n";
print IP2GO Data::Dumper->Dump([\%interpro_GO]);
close IP2GO;



$log->mail();
print "Finished.\n" if ($verbose);
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

=item None known.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
