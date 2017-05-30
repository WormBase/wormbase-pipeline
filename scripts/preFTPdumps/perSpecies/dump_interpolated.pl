#!/usr/bin/env perl

# This dumps out the table of interpolated genetic map positions for
# all of the clones.
# This query is too slow to do on the fly alas!
#
# will only work for C.elegans
#

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

unless(lc($wormbase->species) eq 'elegans'){
   $log->write_to("skipping ... as it only works for C.elegans\n");
   $log->mail;
   exit 0;
}

# Establish a connection to the database.
$log->write_to("connecting to ${\$wormbase->autoace}\n");
my $db = Ace->connect(-path => $wormbase->autoace )||die Ace->error;

my $file = $wormbase->reports . '/'.
   join('.',$wormbase->gspecies_name,$wormbase->ncbi_bioproject,'WSXXX.interpolated_clones.txt');

my $of = IO::File->new($file,'w');
$log->write_to("writing to $file\n");


my $query = 
'select g,g->Interpolated_map_position,g->Interpolated_map_position[2] from g '.
'in class "Genome_sequence" where exists_tag g->Interpolated_map_position order by :2 asc, :3 asc';

my @rows = $db->aql($query);

map {print $of,join("\t",@$_),"\n"} @rows;

$log->mail;
$of->close;
