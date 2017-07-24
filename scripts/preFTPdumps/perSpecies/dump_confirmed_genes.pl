#!/usr/bin/env perl

use strict;

use IO::File;
use Storable;
use Getopt::Long;
use Bio::SeqIO;

use lib $ENV{CVS_DIR};

use Wormbase;
use Log_files;

my ($store,$debug,$test,$database,$species,$outfile);
GetOptions(
       'store=s' => \$store,
       'debug=s' => \$debug,
       'test'    => \$test,
       'species=s'  => \$species,
       'database=s' => \$database,
       'outfile=s'  => \$outfile,
)||die(@!);

my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n");
} 
else {
  $wormbase = Wormbase->new( -debug => $debug, 
                             -test => $test,
                             -organism => $species);
}

my $log = Log_files->make_build_log($wormbase);
$database = $wormbase->autoace if not defined $database;

# Establish a connection to the database.
$log->write_to("connecting to $database\n");
my $db = Ace->connect(-path => $database ) or $log->log_and_die(Ace->error);

$outfile = $wormbase->reports . '/confirmed_genes.fa'
    if not defined $outfile;
$log->write_to("writing to $outfile\n");

my $WS_version = $wormbase->get_wormbase_version();
my $wormdna = $wormbase->wormpep."/".$wormbase->pepdir_prefix."pep.dna${WS_version}";
my $wormpep = $wormbase->wormpep."/".$wormbase->pepdir_prefix."pep${WS_version}";

my $seqio_pep  = Bio::SeqIO->new('-file' => "$wormpep", '-format' => 'fasta');
my %conf_genes;
while(my $seq = $seqio_pep->next_seq) {
  my ($status) = $seq->desc =~ /status=(\S+)/;
  if ($status =~ /^confirmed$/i) {
    $conf_genes{$seq->id} = 1;
  }
}

my $seqio  = Bio::SeqIO->new('-file' => "$wormdna", '-format' => 'fasta');
my $seqio_out = Bio::SeqIO->new('-file' => ">$outfile", '-format' => 'fasta');

while(my $seq = $seqio->next_seq){
  if($conf_genes{$seq->id}) {
    $seqio_out->write_seq( $seq );
  }
}
$db->close;
$log->write_to("Finished extracting confirmed genes\n\n");

$log->mail;
exit(0);
