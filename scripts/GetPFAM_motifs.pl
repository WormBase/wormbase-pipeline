#!/usr/bin/env perl
#
# GetPFAM_motifs.pl 
# 
# by Anthony Rogers
#
# Gets latest PFAM motifs from sanger/pub and puts info in to ace file
#
# Last updated by: $Author: klh $                      
# Last updated on: $Date: 2013-10-14 10:16:23 $         

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase,$noload, $acefile);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "noload"     => \$noload,
            "acefile=s"  => \$acefile,
            );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}
my $log = Log_files->make_build_log($wormbase);

$acefile = $wormbase->acefiles . "/pfam_motifs.ace" if not defined $acefile;
my $pfam_motifs = "/tmp/Pfam_motifs.".$wormbase->species.".gz";
$wormbase->run_command("wget -q -O $pfam_motifs ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz", $log);

open(my $pfam_fh, "gunzip -c $pfam_motifs |")  or 
    $log->log_and_die("Could not open gunzip stream to $pfam_motifs");

open (my $out_fh,">$acefile") or $log->log_and_die("Cannot open $acefile for writing\n");

my (%rec);

my $pfcount = 0;
while (<$pfam_fh>){

  if(/^\#=GF AC\s+(PF\d{5})/  ){ 
    $rec{acc} = $1;
    next;
  } 

  if(/^\#=GF ID\s+(.*$)/  ) {
    my $id = $1;
    $id =~ s/\"//g;
    $rec{id} = $id;
    next;
  }    

  if($_ =~ m/^\#=GF DE\s+(.*$)/  ) {
    my $title = $1;
    $title =~ s/\"//g;
    $rec{title} = $title;
    next;
  } 


  if (/^\/\//){
    if (exists $rec{acc} and exists $rec{id} and exists $rec{title}) {      
      print $out_fh "Motif : \"PFAM:$rec{acc}\"\n";
      print $out_fh "Title \"$rec{title}\"\n";
      print $out_fh "Database \"Pfam\" \"Pfam_ID\" \"$rec{acc}\"\n";
      print $out_fh "Database \"Pfam\" \"short_name\" \"$rec{id}\"\n";
      print $out_fh "\n";
      $pfcount++;
    } else{
      $log->write_to("Skipping $rec{acc} due to missing id/desc\n");
    }
    %rec = ();
  }
}
$log->write_to("added $pfcount PFAM motifs\n");
unlink $pfam_motifs;

close($out_fh) or $log->log_and_die("Could not close $acefile after writing\n");

unless ($noload) {
  $wormbase->load_to_database($wormbase->autoace, $acefile, 'pfam_motifs', $log);
}


$log->mail();
exit(0);


__END__
