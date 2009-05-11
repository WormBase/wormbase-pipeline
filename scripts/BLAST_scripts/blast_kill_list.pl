#!/usr/bin/env perl

use GDBM_File;
use lib '/software/worm/ensembl/bioperl-live/';
use Bio::SeqIO;
use Getopt::Long;
use strict;

my %id_hash;
my ($infile,$killfile,$outfile,$debug);
GetOptions(
    'infile=s'   => \$infile,
    'killfile=s' => \$killfile,
    'outfile=s'  => \$outfile,
    'debug'      => \$debug,
)||die($!);

&prepare_db($killfile);

my $seqio_object = Bio::SeqIO->new(-file => $infile);
my $seqout = Bio::SeqIO->new(-file => ">$outfile",-format => 'Fasta');

while(my $seq=$seqio_object->next_seq){
    if ($id_hash{$seq->display_id} eq _crc64($seq->seq)){
        printf "killing %s %s\n", $seq->display_id, $id_hash{$seq->display_id} if $debug;
        next;
    }
    $seqout->write_seq($seq);
}

# create flatfile database to speed up lookups
sub prepare_db{
    my ($file)=@_;
    tie %id_hash,'GDBM_File', '/tmp/id_hash',&GDBM_WRCREAT, 0666 or die "cannot open /tmp/id_hash DBM file\n";
    open INF, "<$file" || die "cannot open $file\n";
    my $counter=0;
    while (<INF>){
        my @a=split;
        $id_hash{$a[0]}=$a[-1];
        print '.' if (($counter++ % 1000 == 0)&&$debug);
    }
    print "\n" if $debug;
    close INF;
}

# lifted from Renee Baecker
sub _crc64 {
  my ($text) = @_;
  use constant EXP => 0xd8000000;
  my @highCrcTable = 256;
  my @lowCrcTable  = 256;
  my $initialized  = ();
  my $low          = 0;
  my $high         = 0;

  unless($initialized) {
    $initialized = 1;
    for my $i(0..255) {
      my $low_part  = $i;
      my $high_part = 0;
      for my $j(0..7) {
        my $flag = $low_part & 1; # rflag ist f<C3><BC>r alle ungeraden zahlen 1
        $low_part >>= 1;# um ein bit nach rechts verschieben
        $low_part |= (1 << 31) if $high_part & 1; # bitweises oder mit 2147483648 (), wenn $parth ungerade
        $high_part >>= 1; # um ein bit nach rechtsverschieben
        $high_part ^= EXP if $flag;
      }
      $highCrcTable[$i] = $high_part;
      $lowCrcTable[$i]  = $low_part;
    }
  }

  foreach (split '', $text) {
    my $shr = ($high & 0xFF) << 24;
    my $tmph = $high >> 8;
    my $tmpl = ($low >> 8) | $shr;
    my $index = ($low ^ (unpack "C", $_)) & 0xFF;
    $high = $tmph ^ $highCrcTable[$index];
    $low  = $tmpl ^ $lowCrcTable[$index];
  }
  return sprintf("%08X%08X", $high, $low);
}

# cleanup the tempfile
END{
    untie %id_hash;
    unlink '/tmp/id_hash';
}