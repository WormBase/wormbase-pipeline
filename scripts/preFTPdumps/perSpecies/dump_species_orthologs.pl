#!/usr/bin/env perl

use strict;

use IO::File;
use Storable;
use Getopt::Long;

use lib $ENV{CVS_DIR};
use lib "$ENV{CVS_DIR}/preFTPdumps/perSpecies";

use Dumper;
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
  $wormbase = Wormbase->new( -debug   => $debug, 
                             -test    => $test,
                             -organism => $species);
}


my $log = Log_files->make_build_log($wormbase);

$database = $wormbase->autoace if not defined $database;

$outfile = $wormbase->reports . '/orthologs.txt'
    if not defined $outfile;
my $of = IO::File->new($outfile,'w');
$log->write_to("writing to $outfile\n");

my $full_name = $wormbase->long_name;

print $of 
"# $full_name orthologs\n".
    "# WormBase version: " . $wormbase->get_wormbase_version . "\n".
    '# Generated:'.get_date."\n".
    '# File is in record format with records separated by "=\n"'."\n".
    "#      Sample Record\n".
    '#      WBGeneID \t PublicName \n'."\n".
    '#      Species \t Ortholog \t Public_name \t MethodsUsedToAssignOrtholog \n'."\n".
    '# BEGIN CONTENTS'."\n";

$log->write_to("connecting to $database\n");
my $dbh = Ace->connect(-path => $database ) or $log->log_and_die("Connection failure: ".  Ace->error );

$log->write_to("Querying for genes...\n");
my @genes = $dbh->fetch(-query=>"find Gene WHERE Ortholog AND WBGene* AND Species=\"$full_name\"");
foreach my $gene (@genes) {
  my $nm = $gene->Public_name;
  my @orthologs = $gene->Ortholog;
  next unless @orthologs;

  print $of "=\n";
  print $of "$gene\t$nm\n";
  
  foreach my $o ($gene->Ortholog) {
    my @support = $o->col(3);
    @support = $o->col(2) unless @support;
    my $support_str = join(";", @support);

    next if not $o->Species;
    next if not $o->Public_name;
    
    print $of join("\t",$o->Species,$o,$o->Public_name,$support_str), "\n";
  }
}

$dbh->close();
$of->close;
$log->mail;

exit(0);
