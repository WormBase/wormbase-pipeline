#!/usr/bin/env perl
use Getopt::Long;
use IO::File;
use Wormbase;

use feature 'say';

my $date = `date`;
my ($store,$worm,$debug,$test,$out,$wormbase);
chomp $date;

GetOptions(
	'fullspecies=s'  => \$worm,    # $worm->species
	'store=s'        => \$store,
        'debug=s'        => \$debug,
        'test'           => \$test,
        'out=s'          => \$out,
)||die(@!);


if ( $store ) {
  $wormbase = Storable::retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             );
}
my $log = Log_files->make_build_log($wormbase);
my $outfh = IO::File->new($out,'w')||die(@!);

$worm =~ s/_/ /;
$log->write_to("Dumping orthologs for $worm to $out\n");

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
  my $pn = $gene->Public_name;
  next if not $pn;

  my @orth_strings;

  foreach my $o($gene->Ortholog){
    next if not $o->Species;
    my @support = $o->col(3);
    push @orth_strings, join("\t",$o->Species, $o, $o->Public_name, join(';',@support)));
  }
  
  if (@orth_strings) {
    say $outfh '=';
    say $outfh "$gene\t$pn\n";
    foreach my $str (@orth_strings) {
      say $outfh $str;
    }
  }
}

$db->close;
$log->mail;
exit(0);
