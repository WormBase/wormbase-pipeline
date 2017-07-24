#!/usr/bin/env perl

use strict;

use Dumper;
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
  $wormbase = Wormbase->new( -debug   => $debug, 
                             -test    => $test,
                             -organism => $species);
}


my $log = Log_files->make_build_log($wormbase);

$database = $wormbase->autoace if not defined $database;

$log->write_to("connecting to ${\$wormbase->autoace}\n");
my $dbh = Ace->connect(-path => $database ) or $log->log_and_die("Connection failure: ".  Ace->error );

$outfile = $wormbase->reports . '/orthologs.txt'
    if not defined $outfile;
my $of = IO::File->new($outfile,'w');
$log->write_to("writing to $outfile\n");

my $full_name = $wormbase->long_name;

print $of 
"# $full_name orthologs\n".
"# WormBase version: " . $dbh->version . "\n".
'# Generated:'.get_date."\n".
'# File is in record format with records separated by "=\n"'."\n".
"#      Sample Record\n".
'#      WBGeneID \t PublicName \n'."\n".
'#      Species \t Ortholog \t Public_name \t MethodsUsedToAssignOrtholog \n'."\n".
'# BEGIN CONTENTS'."\n".

my $i = $dbh->fetch_many(-query=>"find Gene Species=\"$full_name\"");
while (my $gene = $i->next) {
  my $nm = $gene->Public_name;
  my @orthologs = $gene->Ortholog;
  next unless @orthologs;

  print $of "=\n";
  print $of "$gene\t$nm\n";
  
  foreach my $o ($gene->Ortholog) {
    my @support = $o->col(3);
    @support = $o->col(2) unless @support;
    my $support_str = join(";", @support);
    
    print $of join("\t",$o->Species,$o,$o->Public_name,$support_str), "\n";
  }
}

$of->close;
$log->mail;

