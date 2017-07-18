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
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n");
} 
else {
  $wormbase = Wormbase->new( -debug => $debug, -test => $test,-autoace=> $database,-organism => $species);
}

my $log = Log_files->make_build_log($wormbase);

$outfile = $wormbase->reports . '/cdna2orf.fa'
    if not defined $outfile;
my $of = IO::File->new($outfile,'w');
$log->write_to("writing to $outfile\n");

my $tace = $wormbase->tace;
my $ace_dir = (defined $database) ? $database : $wormbase->autoace;
my $command=<<EOF;
Table-maker -p "$ace_dir/wquery/cDNA2CDS.def" quit
EOF

my $dir = "$ace_dir";

my %cDNA2orf;
open (TACE, "echo '$command' | $tace $dir | ") or $log->log_and_die("Couldn't access $dir\n");  
while (<TACE>){
  chomp;
  if (m/^\"/){
    s/\"//g;
    m/(.+)\s+(.+)/;     
    $cDNA2orf{$1} .= "$2 ";
  }
}
close(TACE);
  
foreach my $key (keys %cDNA2orf){
  print $of "$key,$cDNA2orf{$key}\n";
}
$of->close();

$log->mail;
exit(0);
