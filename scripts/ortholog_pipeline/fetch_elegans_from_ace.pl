#!/usr/bin/perl

# Retrieves all elegans genes from the database.
# This may include additional superlinked genes and alternative splices

# This works for everything that I need except that it does not
# Correctly report the sequences in clone coordinates

use Ace;
use Ace::Sequence;
use DBI;
use strict;

my $ELEGANS_GENES = shift;
chomp $ELEGANS_GENES;

my $path = shift || '/data_2/acedb/elegans';
my $counter;

# NOTE: THESE MAY SHIFT!
my %chrom_lengths = (
		     CHROMOSOME_I   => '15080475',
		     CHROMOSOME_II  => '15174006',
		     CHROMOSOME_III => '13855084',
		     CHROMOSOME_IV  => '17493985',
		     CHROMOSOME_V   => '20916337',
		     CHROMOSOME_X   => '17749735'
		    );

# A running counter of the genome position.
my $genome_position = 0;

my $date = `date`;
print "#",join("\t",qw/gene protein clone chrom strand pep_length spliced_length start stop chrom_start chrom_stop 
	       genome_start genome_stop gmap source version method/),"\n";


# Fetch a list of all <clean> elegans genes
my $ELEGANS = fetch_clean();

# connect to database
my $db = Ace->connect(-path=>$path) || die "Couldn't open database $!";

# find all genomic sequences that are linked to a chromosome
my @sequences = $db->fetch('Sequence' => 'CHROMOSOME*');
die "Couldn't get any genome sequences" unless @sequences;

# recursively follow down to the gene level
for my $s (@sequences) {
  print STDERR "processing genes on $s...\n";
  my @genes = dump_seqs($s,"$s",1);

  $genome_position += $chrom_lengths{$s};
}

sub dump_seqs {
  my ($seq,$chromosome,$offset,$len,%in) = @_;
  # If we're at a terminal sequence (i.e. a predicted gene), then dump it.
  my @subseq = sort {$a->right <=> $b->right} $seq->get('Subsequence');
  if ( !@subseq ) {
    dump_gene($seq,$chromosome,$offset,$offset+$len-1,%in);
    return;
  }
  
  # otherwise recursively fetch information about segments
  foreach (@subseq) {
    my %info;
    my $s = $_->fetch;
    
    # convert start/end coordinates into cumulative offset and length
    my $start  = $offset + $_->right-1;    # offset to beginning of sequence
    my $length = $_->right(2) - $_->right(1) + 1;
    
    # This section accumulates information that should be propagated
    # from links and superlinks down to the gene level.
    $info{gbk}  = $s->Database(3) if $s->Database;
    if ($s->Interpolated_gMap) {
      $info{gmap} = join "\t",$s->Interpolated_gMap(1)->row;
    }
    
    $info{cosmid} = $s unless $s =~ /\.\w+$/;
    if (my $d = $s->DNA) {
      if (my $dna = $s->asDNA) {
	$dna =~ s/^>.+//;
	$dna =~ s/\n//g;
	($info{clength},$info{cgc}) = get_gc_content($dna);
	$info{dna}     = $dna;
      }
    }
    
    # remember where the current sequence begins in its parent
    $info{start} = $_->right(1);
    $info{end}   = $_->right(2);
    
    # dump subsequences
    dump_seqs($s,$chromosome,$start,$length,%in,%info);
  }
}


sub fetch_clean {
  open IN,$ELEGANS_GENES;
  my %h;
  while (<IN>) {
    chomp;
    # Let's only look at the gene, not at alternative splices...
    # my $clean = ($_ =~ /(.*\.\d+).*/) ? $1 : $_;
    $h{uc $_}++;
  }
  return \%h;
}



# Dump out information about a particular gene and its corresponding locus
sub dump_gene {
  my ($gene,$chromosome,$start,$end,%info) = @_;
  return unless $gene->Coding;

  # Is this gene one of the true clean genes from ws77?
  return unless ($ELEGANS->{uc $gene});
  my $clean_gene = ($gene =~ /(.*\.\d+).*/) ? $1 : $gene;
  my $id      = $gene->Brief_identification;
  my $locus   = $gene->Locus_genomic_seq;
  my ($gpos,$gconf);
  if ($locus) {
    $gpos = join "\t",$locus->Map,$locus->Map(3);
    $gconf = $locus->Map(5);
  }
  
  my $seq     = Ace::Sequence->new($gene);
  my ($dna);
  if ($info{dna}) {
    my ($start,$end) = @info{'start','end'};
    ($start,$end) = ($end,$start) if $end < $start;
    $dna     = substr($info{dna},$start-1,$end-$start+1);
  } else {
    $dna = $seq->dna;
  }
  my $strand = $seq->strand;
  my ($gc_content,$length) = get_gc_content($dna);

  $chromosome =~ s/^CHROMOSOME_//;
  $info{gmap} = ($info{gmap} =~ /[I|II|III|IV|V|X]\s(.*)/) 
    ? $1 : $info{gmap};
		 
  # $info{start} and $info{end} contain the start and end 
  # positions of the gene in clone coords 
  # (basically, but they are slightly off - works for this analysis, though)
  
  # Fetch the spliced sequence
  my @exons = $seq->features('exon');
  my $spliced;
  #  my $exonct = 1;
  for (sort { $a->start <=> $b->start } @exons) {
    #    print $_->id(), "\n";
    #    print "Exon: $exonct\n",
    #      to_fasta($_->dna, $exonct++);
    $spliced .= $_->dna;
  }
  
  # Fetch and format the peptide sequence
  my $pep     = $gene->asPeptide();
  my @lines = split ("\n",$pep);
  my $pep;
  foreach my $l (@lines) {
    next if $l =~ /^>.*/;
    $pep .= $l;
  }
  
  # Calculate the running genomic position for this gene...
  my $genome_start = $start + $genome_position;
  my $genome_stop  = $end  + $genome_position;
  print join("\t",$clean_gene,uc $gene,$info{cosmid},$chromosome,$strand,length $pep,
	     length $spliced,
	     $info{start},$info{end},
	     $start,$end,
	     $genome_start,$genome_stop,
	     $info{gmap},'elegans','-','predicted'),"\n";
}


sub to_fasta {
  my ($seq,$id) = @_;
  my $fasta;
  for (my $i=0;$i<=length $seq;$i+=80) {
    my $chunk = substr($seq,$i,60);
    $fasta .= $chunk;
    $fasta .= "\n" if ($chunk ne '');
  }
  return ($fasta);
}


sub get_gc_content {
  my $dna = shift;
  my $gc      = $dna =~ tr/gcGC/gcGC/;
  my $length     = length $dna || '';
  my $gc_content = $length > 0 ? sprintf("%.2f",$gc/$length): '';
  return ($gc_content,$length);
}

