#!/usr/bin/env perl
use Getopt::Long;
use IO::File;
use Wormbase;

use feature 'say';

my $date = `date`;
my ($store,$worm,$debug,$test,$out,$wormbase);
chomp $date;

GetOptions(
	'species=s'    => \$worm,    # $worm->species
	'store=s'      => \$store,
        'debug=s'      => \$debug,
        'test'         => \$test,
        'out=s'        => \$out,
)||die(@!);


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             );
}
my $log = Log_files->make_build_log($wormbase);
my $outfh = IO::File->new($out)||die(@!);

$log->write_to("Dumping orthologs for ${\$wormbase->species} to $out\n");

my $db = Ace->connect(-path => $wormbase->autoace) ||die (Ace::Error);

my $version = 'WS'.$wormbase->version;


print $outfh <<HERE;
# $worm orthologs
# WormBase version: $version
# Generated: $date
# File is in record format with records separated by "=\\n"
#      Sample Record
#      WBGeneID \\t PublicName \\n
#      Species \\t Ortholog \\t Public_name \\t MethodsUsedToAssignOrtholog \\n
# BEGIN CONTENTS
HERE

my $geneIt = $db->fetch_many(-query => "Find Gene; Species=\"$worm\"");
while (my $gene = $geneIt->next){
   my @orthologs = ($gene->Ortholog,$gene->Ortholog_other);
   next unless @orthologs;
   say $outfh '=';
   say $outfh "$gene\t${\$gene->Public_name}";
   foreach my $o(@orthologs){
        my @support = $o->col(3);
        @support = $o->col(2) unless @support;
	say $outfh join("\t",($o->Species,$o,$o->class eq 'Gene'?$o->Public_name||'""':'""',join(';',@support)));
   }
}

$log->mail;
