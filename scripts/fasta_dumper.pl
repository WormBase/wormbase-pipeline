#!/software/bin/perl -w
#
# fasta_dumper.pl
#
# Script to dump a fasta file for a specified species/class/method
# Option to add WBGeneIDs to the header of the file.
# >2L52.1 gene=WBGene00007063
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2014-06-17 15:27:47 $
use strict;
use lib $ENV{'CVS_DIR'};
use Getopt::Long;
use Wormbase;
use File::Path;
use Storable;

use Bio::PrimarySeq;
use Bio::SeqIO;

my ($debug, $store, $verbose, $database, $test, $wormbase, $species, 
    @classmethodlabel, $out, @seqs, $pep);

GetOptions ( "debug:s"               => \$debug, 
	     "verbose"               => \$verbose,           #verbose quces a little more info to screen
	     "classmethodlabel=s@"   => \@classmethodlabel,  #specify an acedb class
	     "database:s"            => \$database,   #specify a database to query
	     "output:s"              => \$out,        #specify an output file 
	     "test"                  => \$test,       #invoke test env
	     "store:s"               => \$store,      #supply a storable
	     "species:s"             => \$species,    #needed to work out what species is being processed
             "pep"                   => \$pep         #peptide dump 
	   ) ;



if ( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			     -organism => $species,
                           );
}
# Create log file
my $log = Log_files->make_build_log($wormbase);

# Define database to be queried
$database = $wormbase->autoace if not defined $database;
$out = $wormbase->sequences . "/coding_transcripts.dna" if not defined $out;

my %gene_ids = $wormbase->FetchData('worm_gene2geneID_name');

if (not @classmethodlabel) {
  $log->write_to("!!!!Warning - No class:method specified so defaulting to Transcript:Coding_transcript\n");
  push @classmethodlabel,  "Transcript:Coding_transcript";
}

my $full_name = $wormbase->full_name;

my $connection = Ace->connect(-path => $database) || die (Ace->error);
foreach my $classmethodlabel (@classmethodlabel) {
  my ($class, $method, $label) = split(/:/, $classmethodlabel);

  my $query = "Find $class WHERE Species=\"${full_name}\""; 
  if (defined $method) {
    $query .= " AND Method=\"$method\"";
  }

  my $object_it = $connection->fetch_many(-query => $query);
  while(my $object = $object_it->next){
      my $dna;
      if ($pep) {
          $dna = $object->asPeptide();
      }
      else {
          $dna = $object->asDNA();
      }
    my @dna = split(/\n/, $dna);
    shift @dna;
    $dna = join("", @dna);
    chomp($dna);
    
    my $seq = Bio::PrimarySeq->new(-seq => uc($dna),
                                   -id  => $object->name);
    my @desc;
    if (defined $label) {
      push @desc, "type=$label";
    }
    if (exists $gene_ids{$object->name}) {
      push @desc, sprintf("gene=%s", $gene_ids{$object->name})
    }
    if (@desc) {
      $seq->description(join(" ", @desc));
    }

    push @seqs, $seq; 
  }
}

my $seqio = Bio::SeqIO->new(-format => 'fasta',
                            -file => ">$out");
foreach my $seq (sort { $a->id cmp $b->id } @seqs) {
  $seqio->write_seq($seq);
}
$seqio->close();
$log->write_to("Wrote " . scalar(@seqs) . " seqs for @classmethodlabel\n");

$connection->close();

$log->mail();
exit(0);

__END__
