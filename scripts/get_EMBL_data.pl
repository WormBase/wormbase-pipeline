#!/software/bin/perl -w
#
# get_embl_data
# 
# Usage : get_embl_data.pl 
#
# Reads protein ids and gets SwissProt IDs
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2011-11-08 16:19:06 $

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

my ($help, $debug, $test, $verbose, $store, $wormbase, $species, $table, $infile, $autoace, $no_load, $acefile);

GetOptions (
  "debug=s"    => \$debug,
  "test"       => \$test,
  "store:s"    => \$store,
  "species:s"  => \$species,
  "infile:s"   => \$infile,
  "autoace=s"  => \$autoace,
  "noload"     => \$no_load,
  "acefile=s"  => \$acefile,
           );



if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism=> $species,
                             -autoace => $autoace,
			     );
}


# establish log file.
my $log = Log_files->make_build_log($wormbase);

my %cds2gene_id = $wormbase->FetchData('cds2wbgene_id');
my %cds2wormpep = $wormbase->FetchData('cds2wormpep');
my %accession2clone   = $wormbase->FetchData('accession2clone');

$acefile = $wormbase->acefiles."/cds_embl_data.ace" if not defined $acefile;
$infile = $wormbase->autoace .  "/protein_ID.mail" if not defined $infile;

my $acc2idmap = &get_acc_to_id_map();

#open output
open(ACE,">$acefile") or $log->log_and_die("cant write $acefile : $!\n");
open (IN,"<$infile") or $log->log_and_die("cant read $infile\n");

while(<IN>){
  #U10401  105     AAA19054        4       1408728171      T20B12.1        P41842  T20B12.1
  s/[\n\r]//g;
  my @data = split("\t",$_);
  
  next unless scalar(@data) == 8;
  my($cloneacc, $pid, $version, $cds, $uniprot) = ($data[0],$data[2],$data[3],$data[-1],$data[-2]);  
  
  next unless (defined $pid);
  print "Potential New Protein: $_\n" if $uniprot eq 'UNDEFINED';
  
  next unless $accession2clone{$cloneacc}; #mail includes some mRNAs
  
  print ACE "\nCDS : \"$cds\"\n";
  print ACE "Protein_id ".$accession2clone{$cloneacc}." $pid $version\n";
  if($cds2wormpep{$cds} and defined $uniprot and $uniprot ne 'UNDEFINED') {

    print ACE "\nProtein : \"WP:".$cds2wormpep{$cds}."\"\n";
    print ACE "Database UniProt UniProtAcc $uniprot\n" if $uniprot;
    if (exists $acc2idmap->{$uniprot}) {
      print ACE "Database UniProt UniProtID ".$acc2idmap->{$uniprot}."\n";
    }
    
    print ACE "\nCDS : \"$cds\"\n";
    print ACE "Database UniProt UniProtAcc $uniprot\n" if $uniprot;
    if (exists $acc2idmap->{$uniprot}) {
      print ACE "Database UniProt UniProtID ".$acc2idmap->{$uniprot}."\n";
    }
  }
  else {
    $log->write_to("no ".$wormbase->pepdir_prefix."pep for $cds\n");
  }
}	

close(ACE) or $log->log_and_die("Could not close $acefile properly\n");

unless ($no_load or $test) {
  $wormbase->load_to_database($wormbase->autoace, $acefile, 'EMBL_ids',$log);
}

$log->mail;


###################
sub get_acc_to_id_map {
  my $query = "mfetch -d uniprot -i \"org:".$wormbase->full_name."\" -f \"div acc\"";
  open (IDS ,"$query |") or $log->log_and_die("mfetch id query failed : $!\n");

  my ($id, %uac2uid);

  while(<IDS>){
    if(/^ID\s+(\w+)\s+/) {
      $id = $1;
      next;
    }
    if (/^AC/) {
      foreach my $acc (/(\w+);/g) {
        $uac2uid{$acc} = $id;
      }
    }
  }

  return \%uac2uid;
}
