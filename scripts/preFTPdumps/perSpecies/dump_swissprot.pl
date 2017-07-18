#!/usr/bin/env perl

use strict;

use IO::File;
use Storable;
use Getopt::Long;

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
if ($store) { $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")} 
else {$wormbase = Wormbase->new( -debug => $debug, -test => $test,-autoace=> $database,-organism => $species)}

my $log = Log_files->make_build_log($wormbase);

# Establish a connection to the database.
$log->write_to("connecting to ${\$wormbase->autoace}\n");
my $db = Ace->connect(-path => $wormbase->autoace )||die Ace->error;

$outfile = $wormbase->reports . '/swissprot.txt'
    if not defined $outfile;
my $of = IO::File->new($outfile,'w');
$log->write_to("writing to $outfile\n");

print $of "#" . join("\t",qw(WormPep CDS UniProtAc status Description)) . "\n";

my $iterator = $db->fetch_many(-query=>"find Protein Species=\"${\$wormbase->long_name}\"; Live",-filltag=>'Database');
while (my $p = $iterator->next) {
   my $swissprot_acc  = $p->Database?$p->Database(0)->at('SwissProt.UniProtAcc[1]'):'';
   my $uniprot_acc = $p->Database?$p->Database(0)->at('UniProt.UniProtAcc[1]'):'';
   foreach my $cds ($p->Corresponding_CDS){
     my $desc = $cds->Brief_identification;
     print $of join("\t",
                    "$p",
                    $cds,
                    $uniprot_acc,
                    ($swissprot_acc) ? "REVIEWED" : "UNREVIEWED",
                    $desc),"\n";
   }
}

$log->mail;
$of->close;
exit(0);
