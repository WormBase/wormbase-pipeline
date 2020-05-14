#!/usr/bin/env perl
#
# make_wormrna.pl
# 
# Builds a wormrna data set from the current autoace database
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2015-05-11 12:35:47 $

use strict;

use Getopt::Long;
use Storable;
use Ace;
use Bio::SeqIO;
use Digest::MD5 qw(md5 md5_hex md5_base64) ;

use lib $ENV{CVS_DIR};
use Wormbase;

my $rnacentral_md5s = "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/md5/md5.tsv.gz";

my ($help, $debug, $test, $store, $species, $noload);
my ($rnafile, $no_rnacentral, %seqs);

GetOptions (
  "help"         => \$help,
  "debug=s"      => \$debug,
  "test"         => \$test,
  "store:s"      => \$store,
  "species:s"    => \$species,
  "rnafile=s"    => \$rnafile,
  "nornacentral" => \$no_rnacentral,
  "noload"       => \$noload,
	   );

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
                             -organism => $species
			   );
}

my $log = Log_files->make_build_log($wormbase);


#######################################
# release data                        #
#######################################
my $dbdir     = $wormbase->autoace;
my $tace      = $wormbase->tace;

my $release   = $wormbase->get_wormbase_version;
my $new_wrdir = $wormbase->wormrna;

$rnafile = "$new_wrdir/".$wormbase->pepdir_prefix."rna$release.rna" 
    if not defined $rnafile;

###############################################
# retrieve the desired RNA sequence objects   #
###############################################
$log->write_to("Connecting to database...\n") if $debug;

my $db = Ace->connect (-path => $dbdir, -program => $tace) 
    or $log->log_and_die("Could not connect to $dbdir\n");

# Get RNA genes, but not other Transcript objects
my $query = "FIND Transcript WHERE Method != Coding_transcript"
    . " AND Method != history_transcript"
    . " AND Species = \"".$wormbase->full_name."\"";

$log->write_to("Fetching transcripts from db...\n") if $debug;
my $transcripts = $db->fetch_many (-query => $query);


###########################################################################
# get the rna sequence, write a rna.fasta file,
###########################################################################
my $outio = Bio::SeqIO->new(-format => 'fasta',
                            -file   => ">$rnafile");

$log->write_to("Writing transcript sequence...\n") if $debug;
my %skipped;
while( my $obj = $transcripts->next) {    
  
  my $gene = $obj->Gene;
  my $method = $obj->Method;

  if (not $gene) {
    $skipped{$method}++;
    next;
  }
  my $cgc_name = $gene->CGC_name;

  my $brief_id = $method->GFF_feature;
  if (not $brief_id) {
    $log->write_to("No type set for $obj - setting to ncRNA, but this should be fixed\n");
    $brief_id = "ncRNA";
  }

  my $dna = $obj->asDNA();
  if ((!defined ($dna)) || ($dna eq "")) {
    $log->error("ERROR: cannot extract dna sequence for $obj\n");
    # can't include in WormRNA if there is no DNA!
    next; 
  }
  $dna =~ s/\n/ /;
  $dna =~ s/\n//g;
  $dna =~ s/^>\S+\s+//; 
  $dna =~ s/\s//g; 
  $dna =~ tr/a-z/A-Z/;
  
  #  $seqs{$obj->name} = {
  #  md5  => md5_hex($dna),
  #  gene => $gene->name,
  #};
  $seqs{md5_hex($dna)} = {
    gene => $gene->name,
    transcript => $obj->name,
  };


  $dna =~ tr/T/U/; 
  my $desc .= "biotype=$brief_id gene=$gene";
  $desc    .= " locus=$cgc_name" if defined $cgc_name;

  my $seq = Bio::PrimarySeq->new(-seq => $dna,
                                 -id  => $obj,
                                 -desc => $desc);

  $outio->write_seq($seq);
}

$db->close;
$wormbase->check_files($log);

foreach my $meth (sort keys %skipped) {
  $log->write_to(sprintf("Ignored %d transcripts with method %s, because had no gene connection\n", $skipped{$meth}, $meth));
}

&add_rnacentral_xrefs() unless $no_rnacentral;

$log->mail;
exit(0);


##############################################################
sub add_rnacentral_xrefs {
  
  $log->write_to("Adding RNAcentral xrefs...\n") if $debug;

  my $lfile = "/tmp/rnacentral_md5s.wormbase.$$.gz";
  unlink $lfile if -e $lfile;

  $wormbase->run_command("wget -O $lfile $rnacentral_md5s",$log);

  my %ids_by_md5;
  my %genes;
  my %transcripts;
  my $acefile = $wormbase->acefiles . "/rnacentral_xrefs.ace";

  open(my $fh, "zcat $lfile |") or $log->log_and_die("Could not open stream to $lfile\n");
  while(<$fh>) {
    /^(\S+)\s+(\S+)/ and do {
      my ($id,$md5)=("$1","$2");

      if ($seqs{$md5}){
         $genes{$seqs{$md5}->{gene}}||=[];
         push @{$genes{$seqs{$md5}->{gene}}},$id;
         $transcripts{$seqs{$md5}->{transcript}}||=[];
         push @{$transcripts{$seqs{$md5}->{transcript}}},$id;
        }
      }
  }

  unlink $lfile;

  open(my $outfh, ">$acefile") or $log->log_and_die("Could not open $acefile for writing\n");

  while (my($k,$v)=each (%genes)){
	print $outfh "Gene : $k\n";
        map {print $outfh "Database \"RNAcentral\" \"URSid\" \"$_\"\n"} @$v;
	print $outfh "\n";
  }

  while (my($k,$v)=each (%transcripts)){
	print $outfh "Transcript : $k\n";
        map {print $outfh "Database \"RNAcentral\" \"URSid\" \"$_\"\n"} @$v;
	print $outfh "\n";
  }


  close($outfh) or $log->log_and_die("Could not close $acefile after writing\n");

  # To do: load to database
  # 
  if (not $noload) {
    $wormbase->load_to_database($wormbase->autoace, $acefile, "RNAcentral_xrefs", $log);
  }
}














