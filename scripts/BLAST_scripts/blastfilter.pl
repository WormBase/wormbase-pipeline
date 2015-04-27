#!/usr/bin/env perl
#
# takes a fasta file (trembl in our case) and searches the 7 mysql
# blast databases for hits.
# if it can't find any, it will print out ID SPECIES CRC64
#
# needs: to run on the farm
# needs: mysql databases on farmdb1 called worm_ensembl_$SPECIES
# needs: /lustre/scratch101/ensembl/wormpipe/swall_data/trembl2org

use DB_File;
use File::Copy;
use DBI;
use lib $ENV{'WORM_PACKAGES'} . '/ensembl/bioperl-live/';
use Bio::SeqIO;
use strict;

my %ORG;
prepare_gdb();

my %handles = %{&prepare_mysql()};

my $seqio_object = Bio::SeqIO->new(-file => shift);

while (my $seq = $seqio_object->next_seq){
    my $count=check_mysql($seq->display_id);
    next if $count > 0;
    printf "%s\t\"%s\"\t%s\n",$seq->display_id,$ORG{$seq->display_id},_crc64($seq->seq);
}

teardown_gdb();

#########################
# functions
#

DESTROY{
    &teardown_gdb();
}

# GDBM fuffing around
sub prepare_gdb{
    copy "$ENV{'PIPELINE'}/swall_data/trembl2org",'/tmp/trembl2org';
    tie (%ORG, 'DB_File', "/tmp/trembl2org", O_RDWR|O_CREAT, 0666, $DB_HASH) or die "cannot open /tmp/trembl2org DBM file\n";
}
sub teardown_gdb{
    untie %ORG;
    unlink '/tmp/trembl2org';
}

# generates database connections and prepared statements
sub prepare_mysql{
    my %h;
    my @species = ('brugia','pristionchus','japonica','briggsae','remanei','brenneri','elegans','cangaria','hcontortus','mhapla','mincognita','bxylophilus','csp11','csp5','heterorhabditis','loaloa','sratti','trspiralis'.'ovolvulus','sratti');  
    foreach my $key(@species){
        $h{$key}{dbi}=DBI->connect("dbi:mysql:dbname=worm_ensembl_${key};host=$ENV{'WORM_DBHOST'};port=$ENV{'WORM_DBPORT'}",'wormro');
        $h{$key}{sth1}=$h{$key}{dbi}->prepare('SELECT COUNT(protein_align_feature_id) FROM protein_align_feature WHERE hit_name=?');
        $h{$key}{sth2}=$h{$key}{dbi}->prepare('SELECT COUNT(protein_feature_id) FROM protein_feature WHERE hit_name=?');
    }
    return \%h;
}

#check if id exists
sub check_mysql{
    my($id)=@_;
    my $pcount=0;
    foreach my $key(values %handles){
        $$key{sth1}->execute($id);
        $pcount += ($$key{sth1}->fetchrow_array)[0];
        $$key{sth2}->execute($id);
        $pcount += ($$key{sth2}->fetchrow_array)[0];
    }
    return $pcount;
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
        my $flag = $low_part & 1; # rflag ist für alle ungeraden zahlen 1
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
