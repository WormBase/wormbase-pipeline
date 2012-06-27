=pod 

=head1 NAME

Sequence_extract

=head1 SYNOPSIS

 my $seq_obj      = Sequence_extract->invoke($database, $refresh, $wormbase);
 my $seq          = "AH6";
 my $sub_sequence = $seq_obj->Sub_sequence($seq,50,100);


=head1 DESCRIPTION

This object is used to get DNA subsequence of any S_Mapped sequence object in the database

Inherits from Coords_converter

=head1 CONTACT

Anthony  ar2@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut




package Sequence_extract;

use lib $ENV{'CVS_DIR'};
use lib '/software/worm/lib/bioperl-live';
use Carp;
use Wormbase;
use Coords_converter;
use Bio::SeqIO;    
use strict;

our @ISA;
@ISA = ('Coords_converter');


=head2 invoke

  Title   :   invoke
  Usage   :   my $extractor = Sequence_extract->invoke($database,1);
  Function:   Creates the object and loads in the data (generates fresh if requested)
  Returns :   ref to self
  Args    :   Database  (optional) - which database to use. Default is current_DB
              refresh   default is NULL - connect to the database and update coordinates
              wormbase - wormbase object
 Requires:    Database must have CHROMOSOMES directory with dna of each chromosome in CHROMOSOME_*.dna (FASTA format)

=cut


sub invoke {
  my ($class,$database,$refresh,$wormbase) = @_;
  
  # inherit from Coords_converter to get all the coord info
  my $self = {};
  bless $self, $class;
  
  my $cc = Coords_converter->invoke($database, $refresh, $wormbase);
  $self->_coords_converter($cc);

  $database = $wormbase->database('current') unless $database; #defaults to current_DB in CC if not set.
  # get the chromosomal sequences
  my $tace = $wormbase->tace; # <= hmpf
  
  my @chromosome = $wormbase->get_chromosome_names(-mito => 1);
  $self->{chromprefix} = $wormbase->chromosome_prefix;
  if ($wormbase->assembly_type eq 'contig'){
    #contig assemblies	
    my $genome_seq = $wormbase->genome_seq;
    my $seqs = Bio::SeqIO->new(-file => $genome_seq, -format => "fasta");
    while(my $seq = $seqs->next_seq) {
      $self->_sequence_hash->{$seq->id} = $seq->seq;
    }
  }
  else {
    # iterate chromosomes
    my $chr_dir = "$database/CHROMOSOMES";

    foreach my $chromosome (@chromosome){
      #chromosome based assemblies
      my $seqname = $self->{chromprefix} . "$chromosome";
      
      # dump the chromosome if it doesn't exist
      unless( -e "$chr_dir/$seqname.dna" ) {
        open (ACE, "| $tace $database") or croak "cant connect to $database :$!\n";
        
        print "writing DNA seq for $seqname\n";
        print ACE <<EOF;
clear
find sequence $seqname
dna -f $chr_dir/$seqname.dna
EOF
        close ACE;
        $wormbase->remove_blank_lines("$chr_dir/$seqname.dna");
      }
      
      # read the file/sequence into $self
      $/ = "";
      open (SEQ, "$chr_dir/$seqname.dna") or croak "cant open the dna file for $seqname:$!\n";
      my $seq = <SEQ>;
      close SEQ;
      $/ = "\n";
      $seq =~ s/>[\w\-_]+//;
      $seq =~ s/\n//g;
      $self->_sequence_hash->{$seqname} = $seq;
    }
  }
  return $self;
}

=head2 Sub_sequence

  Title   :   Sub_sequence
  Usage   :   $seq_obj->Sub_sequence($seq)          - whole DNA sequence of that object
              $seq_obj->Sub_sequence($seq,50,100)   - 100 bases of $seq starting at index 50 (base 51)
              $seq_obj_>Sub_sequence($seq,50,+100)  - DNA sequence of that object from index 50 (base 51) to 100 bases past the end
              $seq_obj_>Sub_sequence($seq,-50)      - whole DNA sequence of that object with additional 50 bases at the start
              $seq_obj_>Sub_sequence($seq,-50, 500) - 500 bases of $seq starting 50 bases before clone ie -50 to 450

  Function:   Extract a DNA subsequence
  Returns :   DNA sequence as string
  Args    :   Sequence object  (req)
              start coord  - means upstream of seq obj start
              start coord  + means downstream of seq obj end

=cut

sub Sub_sequence {
  my ($self, $seq, $offset, $length) = @_;

  my $subseq = "";
  my $extend = 0;

  my ($chrom_of_clone, $chrom_start_of_clone, $chrom_end_of_clone) = $self->_coords_converter->LocateSpanUp($seq);
  
  my $seq_str = $self->_sequence_hash->{$chrom_of_clone};

  # to extend past end of a seq obj use "+" and no of bases.
  if (defined $offset) {
    $chrom_start_of_clone += $offset;
  }

  if (defined $length) {
    if ($length =~ /^\+(\d+)/) {
      $length = ($chrom_end_of_clone - $chrom_start_of_clone + 1) + $1;
    }
  } else {
    $length = $chrom_end_of_clone - $chrom_start_of_clone + 1;
  }

  my $end = $chrom_start_of_clone + $length - 1;
  if ($end > length($seq_str)) {
    $end = length($seq_str);

    $length = $end - $chrom_start_of_clone + 1;
  }
  
  my $subsequence = substr($seq_str, $chrom_start_of_clone - 1, $length);

  return $subsequence;
}

=head2

  Title   :   DNA_revcomp
  Usage   :   my $revcomp = $seq_obj->($seq)
  Function:   revcomp DNA seq
  Returns :   DNA sequence as string
  Args    :   DNA sequence as string

=cut


sub DNA_revcomp
  {
    my $self = shift;
    my $revseq = reverse shift;
    $revseq =~ tr/a/x/;
    $revseq =~ tr/t/a/;
    $revseq =~ tr/x/t/;
    $revseq =~ tr/g/x/;
    $revseq =~ tr/c/g/;
    $revseq =~ tr/x/c/;

    $revseq =~ tr/A/X/;
    $revseq =~ tr/T/A/;
    $revseq =~ tr/X/T/;
    $revseq =~ tr/G/X/;
    $revseq =~ tr/C/G/;
    $revseq =~ tr/X/C/;
    return ($revseq);
  }

=head2

  Title   :   DNA_rev
  Usage   :   my $rev = $seq_obj->($seq)
  Function:   reverse DNA seq
  Returns :   DNA sequence as string
  Args    :   DNA sequence as string

=cut


sub DNA_rev
  {
    my $self = shift;
    my $revseq = reverse shift;
    return ($revseq);
  }

=head2

  Title   :   DNA_comp
  Usage   :   my $comp = $seq_obj->DNA_comp($seq)
  Function:   complement DNA seq
  Returns :   DNA sequence as string
  Args    :   DNA sequence as string

=cut


sub DNA_comp
  {
    my $self = shift;
    my $seq = shift;
    $seq =~ tr/a/x/;
    $seq =~ tr/t/a/;
    $seq =~ tr/x/t/;
    $seq =~ tr/g/x/;
    $seq =~ tr/c/g/;
    $seq =~ tr/x/c/;

    $seq =~ tr/A/X/;
    $seq =~ tr/T/A/;
    $seq =~ tr/X/T/;
    $seq =~ tr/G/X/;
    $seq =~ tr/C/G/;
    $seq =~ tr/X/C/;
    return ($seq);
  }
  
=head2

  Title   :   flanking_sequences
  Usage   :   my $flanks = $seq_obj->flanking_sequences($seq, 1000, 20000, 30)
  Function:   give flanking sequence of specified location
  Returns :   ref to array contianing two DNA sequences as string
  Args    :   reference seq, start coord, end coord, flank length required.

=cut

sub flanking_sequences {
	my $self = shift;
	my($seq, $start, $end, $length) = @_;
	my @flanks;

	$flanks[0] = $self->Sub_sequence($seq,($start-$length-1),$length);
	$flanks[1] = $self->Sub_sequence($seq,$end,$length);
	
	if( $end < $start ){
		my $tmp = $self->DNA_revcomp($flanks[0]);
		$flanks[0] = $self->DNA_revcomp($flanks[1]);
		$flanks[1] = $tmp;
	}
	
	return \@flanks;
}


#####################
sub _sequence_hash {
  my ($self, $val) = @_;

  if (not exists $self->{_sequence_hash}) {
    $self->{_sequence_hash} = {};
  }
  if (defined $val) {
    $self->{_sequence_hash} = $val;
  }

  return $self->{_sequence_hash};
}

######################
sub _coords_converter {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_coords_converter} = $val;
  }

  return $self->{_coords_converter};
}


1;
