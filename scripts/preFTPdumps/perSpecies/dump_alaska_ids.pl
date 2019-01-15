#!/usr/bin/env perl
# script to dump transcript data for the ALASKA project of CalTech
# Contact is Joseph Min (kmin@caltech.edu)

use strict;
use Getopt::Long;
use IO::File;
use Storable;

use lib $ENV{CVS_DIR};

use Wormbase;
use Log_files;

my ($species,$store,$debug,$test,$database,$outfile);
GetOptions(
     'species=s'  => \$species,
     'store=s'    => \$store,
     'debug=s'    => \$debug,
     'test'       => \$test,
     'database=s' => \$database,
     'outfile=s'  => \$outfile,
)||die(@!);


my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")
}else {
  $wormbase = Wormbase->new( -debug => $debug, 
                             -test => $test,
                             -organism => $species);
}

my $log = Log_files->make_build_log($wormbase);
$database ||= $wormbase->autoace;

$log->write_to("connecting to $database\n");
my $db = Ace->connect(-path => $database) || $log->log_and_die("Could not connect to $database (".Ace->error.")\n");

$outfile ||= $wormbase->reports . "/" . "alaska_ids.tsv";
my $of = IO::File->new($outfile,'w');
$log->write_to("writing to $outfile\n");

warn "finding transcripts...\n";
my @transcripts = $db->fetch(-query => "find Transcript Species=\"${\$wormbase->long_name}\";Sequence;Gene");

$log->write_to("found ",scalar(@transcripts)," transcripts\n");

# dump the data in alaska format columns
# Transcript_id	Gene_id	Public_name	Description	Biotype
for my $t (@transcripts) {
    my $line = join("\t",$t,$t->Gene,$t->Gene->Public_name,'"'.($t->Gene->Concise_description||$t->Gene->Automated_description||'n.a.').'"',$t->Method->GFF_feature);
    print $of "$line\n";
}

$of->close;
$db->close;
$log->mail;
exit(0);
