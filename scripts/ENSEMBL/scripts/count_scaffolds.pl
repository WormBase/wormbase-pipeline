use strict;                                      
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use IO::File;
use Bio::SeqIO;

my ($dbhost, $dbport, $dbuser, $dbname, $mask, $softmask, $outfile, $ebi_header_prefix, $species_string);

GetOptions (
  'host=s'           => \$dbhost,
  'port=s'           => \$dbport,
  'user=s'           => \$dbuser,
  'dbname=s'         => \$dbname,
);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -user   => $dbuser, 
  -host   => $dbhost,
  -port   => $dbport,
  -dbname => $dbname,
    ) or die "Can't connect to Database $dbname";

my $seqio = (defined $outfile) 
    ? Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile")
    : Bio::SeqIO->new(-format => 'fasta', -fh => \*STDOUT);

my @scaffolds;
foreach my $seq ( sort { $b->length <=> $a->length } @{$db->get_SliceAdaptor->fetch_all('toplevel')}) {
  my $seq_name = $seq->seq_region_name;
  push @scaffolds, $seq_name;
}

print scalar @scaffolds;
