#!/usr/bin/env perl
use IO::File;
use Storable;
use Getopt::Long;

use Wormbase;
use Log_files;

use strict;

my ($store,$debug,$test,$database,$species);
GetOptions(
       'store=s' => \$store,
       'debug=s' => \$debug,
       'test'    => \$test,
       'species=s'  => \$species,
       'database=s' => \$database,
)||die(@!);

my $wormbase;
if ($store) { $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")} 
else {$wormbase = Wormbase->new( -debug => $debug, -test => $test,-autoace=> $database,-organism => $species)}

my $log = Log_files->make_build_log($wormbase);

# Establish a connection to the database.
$log->write_to("connecting to ${\$wormbase->autoace}\n");
my $db = Ace->connect(-path => $wormbase->autoace )||die Ace->error;

my $file = $wormbase->reports . '/'.
   join('.',$wormbase->gspecies_name,$wormbase->ncbi_bioproject,'WSXXX.swissprot.txt');

my $of = IO::File->new($file,'w');
$log->write_to("writing to $file\n");

print $of join "\t",qw(#WormPep CDS  SwissProtAc  SwissProtID Description),"\n";

my $iterator = $db->fetch_many(-query=>"find Protein Species=\"${\$wormbase->long_name}\"; WormPep; Live",-filltag=>'Database');
while (my $p = $iterator->next) {
   my $swissprot_id  = $p->Database(0)->at('SwissProt.SwissProt_ID[1]');
   my $swissprot_acc = $p->Database(0)->at('SwissProt.SwissProt_AC[1]');
   foreach my $cds ($p->Corresponding_CDS){
     my $desc          = $cds->Brief_identification;
     print $of join("\t","$p",$cds,$swissprot_id,$swissprot_acc,$desc),"\n";
   }
}

$log->mail;
$of->close;
exit(0);
