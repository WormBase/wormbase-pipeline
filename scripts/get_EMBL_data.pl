#!/software/bin/perl -w
#
# get_embl_data
# 
# Usage : get_embl_data.pl 
#
# Reads protein ids and gets SwissProt IDs
#
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2010-12-21 16:48:51 $

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

##############################
# command-line options       #
##############################

my ($help, $debug, $test, $verbose, $store, $wormbase);

my $initial;             # run script in initial mode at start of build
my $final;               # for full run at end of build
my $table;
my $species = 'elegans';
my $mail;

GetOptions (
            "debug=s"    => \$debug,
            "test"       => \$test,
            "store:s"    => \$store,
            "species:s"  => \$species,
            "mail:s"	 => \$mail
           );



if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism=> $species
			     );
}


# establish log file.
my $log = Log_files->make_build_log($wormbase);

my %cds2gene_id = $wormbase->FetchData('cds2wbgene_id');
my %cds2wormpep = $wormbase->FetchData('cds2wormpep');
my %accession2clone   = $wormbase->FetchData('accession2clone');
my %uac2cds;  #uniprot acc to cds
my %uac2uid;  # acc to ID
&getIDs;

#open output
my $aceoutput = $wormbase->acefiles."/cds_embl_data.ace";
open(ACE,">$aceoutput") or $log->log_and_die("cant write $aceoutput : $!\n");

#read protein id mail
$mail = $mail ? $mail : $wormbase->wormpub."/protein_ID.mail";
open (MAIL,"<$mail") or $log->log_and_die("cant read $mail\n");

while(<MAIL>){
  #U10401  105     AAA19054        4       1408728171      T20B12.1        P41842  T20B12.1
  s/[\n\r]//g;
  my @data = split("\t",$_);
  
  next unless scalar(@data) == 8;
  my($cloneacc, $pid, $version, $cds, $uniprot) = ($data[0],$data[2],$data[3],$data[-1],$data[-2]);  
  
  next unless (defined $pid);
  print "Potential New Protein: $_" unless ($cloneacc and $pid and $version and $cds and $uniprot);
  
  next unless $accession2clone{$cloneacc}; #mail includes some mRNAs
  
  print ACE "\nCDS : \"$cds\"\n";
  print ACE "Protein_id ".$accession2clone{$cloneacc}." $pid $version\n";
  if($cds2wormpep{$cds}) {
    print ACE "\nProtein : \"WP:".$cds2wormpep{$cds}."\"\n";
    print ACE "Database UniProt UniProtAcc $uniprot\n" if $uniprot;
    print ACE "Database UniProt UniProtID ".$uac2uid{$uniprot}."\n" if ($uac2uid{$uniprot});
    
    print ACE "\nCDS : \"$cds\"\n";
    print ACE "Database UniProt UniProtAcc $uniprot\n" if $uniprot;
    print ACE "Database UniProt UniProtID ".$uac2uid{$uniprot}."\n" if ($uac2uid{$uniprot});
  }
  else {
    $log->write_to("no ".$wormbase->pepdir_prefix."pep for $cds\n");
  }
}	

close MAIL or $log->error("didn't close mail properly\n");

$wormbase->load_to_database($wormbase->autoace, $aceoutput, 'EMBL_ids',$log) unless $test;

$log->mail;

sub getIDs {
  my $query = "mfetch -d uniprot -i \"org:".$wormbase->full_name."\" -f \"div acc\"";
  open (IDS ,"$query |") or $log->log_and_die("mfetch id query failed : $!\n");
  #open(IDS,"</tmp/elegans_uniprot");
  #ID   ADD1_CAEEL              Reviewed;         732 AA.
  #AC   Q9U9K0; O44581; Q95X64; Q9U9J9;
  my $id;
  while(<IDS>){
    if(/^ID\s+(\w+)\s+/) {
      $id = $1;
      next;
    }
    foreach my $acc (/(\w+);/g) {
      $uac2uid{$acc} = $id;
    }
    undef $id;
  }
}