#!/usr/bin/env perl

use Getopt::Long;
use IO::File;
use Storable;

use Dumper;
use Wormbase;
use Log_files;

use strict;

my ($species,$format,$store,$debug,$test,$database);
GetOptions(
     'species=s'  => \$species,
     'format=s'   => \$format,
     'store=s'    => \$store,
     'debug=s'    => \$debug,
     'test'       => \$test,
     'database=s' => \$database,
)||die(@!);


my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")
}else {
  $wormbase = Wormbase->new( -debug => $debug, -test => $test,-autoace=> $database,-organism => $species)
}

my $log = Log_files->make_build_log($wormbase);

# Establish a connection to the database.
$log->write_to("connecting to ${\$wormbase->autoace}\n");
my $dbh = Ace->connect(-path => $wormbase->autoace)||$log->log_and_die(Ace->error);

my $file = $wormbase->reports . '/'.
   join('.',$wormbase->gspecies_name,$wormbase->ncbi_bioproject,'WSXXX.resource_gene_ids.txt');
my $of = IO::File->new($file,'w');
$log->write_to("writing to $file\n");

my $count;

print $of 
"# WormBase Gene identifiers\n".
'# WormBase version: '.$dbh->version."\n".
"# NCBI Taxonomy ID \t Species \t GeneID \t Public \t (Locus name) \t Sequence ID\n".
'# Generated: '.&get_date."\n";

my $i = $dbh->fetch_many(-query=>"find Gene Species \"${\$wormbase->long_name}\";Live");
while (my $gene = $i->next) {
    my $name = $gene->Public_name; # CGC or Sequence
    my $molecular_name = join('|',$gene->Molecular_name) ||'';
    my $seqName = $gene->Sequence_name ||'';

    print $of join("\t",$wormbase->ncbi_tax_id,$wormbase->long_name,$gene,$name,$molecular_name,$seqName),"\n";
    $count++;
}

$log->write_to("found $count genes\n");

$log->mail;
$of->close;
exit(0);
