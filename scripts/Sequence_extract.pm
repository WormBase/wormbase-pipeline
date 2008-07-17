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


sub invoke 
  {
    my ($class,$database,$refresh,$wormbase) = @_;

    # inherit from Coords_converter to get all the coord info
    my $self = Coords_converter->invoke($database, $refresh, $wormbase);
    bless $self, $class;

    $database = $wormbase->database('current') unless $database; #defaults to current_DB in CC if not set.
    # get the chromosomal sequences
    my $tace = $wormbase->tace; # <= hmpf
    
    my @chromosome = $wormbase->get_chromosome_names(-mito => 1);
    $self->{chromprefix} = $wormbase->chromosome_prefix;
    if (scalar @chromosome > 100){
    #contig assemblies	
    	my $supercontig_seq = $wormbase->chromosomes."/supercontigs.fa";
    	my $seqs = Bio::SeqIO->new(-file => $supercontig_seq, -format => "fasta");
    	while(my $seq = $seqs->next_seq) {
    		$self->{'SEQUENCE'}->{$seq->id} = $seq->seq;
    	}
    }
    else {
	 # iterate chromosomes
	 foreach my $chromosome (@chromosome){
         #chromosome based assemblies
	  	  my $seqname = $self->{chromprefix} . "$chromosome";

		  # dump the chromosome if it doesn't exist
    	          unless( -e "$database/CHROMOSOMES/$seqname.dna" ) {
      	            open (ACE, "| $tace $database") or croak "cant connect to $database :$!\n";

	            print "writing DNA seq for $seqname\n";
	            print ACE <<EOF;
clear
find sequence $seqname
dna -f $database/CHROMOSOMES/$seqname.dna
EOF
	            close ACE;
                    $wormbase->remove_blank_lines("$database/CHROMOSOMES/$seqname.dna", $log);
	          }

		  # read the file/sequence into $self
	          $/ = "";
	          open (SEQ, "$database/CHROMOSOMES/$seqname.dna") or croak "cant open the dna file for $seqname:$!\n";
	          my $seq = <SEQ>;
	          close SEQ;
	          $/ = "\n";
	          $seq =~ s/>[\w\-_]+//;
	          $seq =~ s/\n//g;
	          $self->{'SEQUENCE'}->{$seqname} = $seq;
	    }

         }
	 return $self;
  }


=head2 Sub_sequence

  Title   :   Sub_sequence
  Usage   :   $seq_obj->Sub_sequence($seq)          - whole DNA sequence of that object
              $seq_obj->Sub_sequence($seq,50,100)   - 100 bases of $seq starting at base 50
              $seq_obj_>Sub_sequence($seq,50,+100)  - DNA sequence of that object from base 50 to 100 bases past the end
              $seq_obj_>Sub_sequence($seq,-50)      - whole DNA sequence of that object with 50 bases at the start ie -50 to end
              $seq_obj_>Sub_sequence($seq,-50, 500) - 500 bases of $seq starting 50 bases before clone ie -50 to 450

  Function:   Extract a DNA subsequence
  Returns :   DNA sequence as string
  Args    :   Sequence object  (req)
              start coord  - means upstream of seq obj start
              start coord  + means downstream of seq obj end

=cut

sub Sub_sequence
  {
    my $self = shift;
    my $seq = shift;
    my $start = shift;
    my $length = shift;
    my $subseq = "";
    my $extend = 0;
    my $chrom;

#    # to extend past end of a seq obj use "+" and no of bases.
    if( $length and $length =~ /\+/) {
      $extend = substr($length,1);
      undef $length;
    }

    #carp "no length given for sequence to extract - using full length of $seq\n" unless (defined $length);

    $seq = $self->{'single_chrom'}->{"$seq"} if $self->{'single_chrom'}->{"$seq"};

    my $prefix = $self->{chromprefix};
    if( $seq =~ /$prefix/ ) {
      $chrom = $seq;

    } elsif( $seq =~ /SUPERLINK/ ) { #passed seq is a SUPERLINK, only elegans uses 'SUPERLINK'
      my $sl = $seq;
      $chrom = $self->_getChromFromSlink("$seq");
      $seq = $chrom; # gets processed as chromosome below

      # modify the starting coordinate
      $start += $self->{"$chrom"}->{SUPERLINK}->{"$sl"}->[0];
      $length = $self->{"$chrom"}->{SUPERLINK}->{"$sl"}->[1] - $self->{"$chrom"}->{SUPERLINK}->{"$sl"}->[0] unless $length;

    } else {
      # This is when the passed seq is a clone
      unless ($self->{'CLONE2CHROM'}->{"$seq"} ) {
	carp "$seq is not a valid sequence\n";
	return;
      }
      SLINKS:
	foreach my $slink (keys %{$self->{'SUPERLINK'}} ) {

	  foreach my $clone (keys %{$self->{SUPERLINK}->{$slink}} ) {

	    if ( "$clone" eq "$seq" ) {
	      $chrom = $self->_getChromFromSlink("$slink");

	      # modify the starting coordinate
	      $start += $self->{"$chrom"}->{SUPERLINK}->{$slink}->[0] - 1; # superlink base coords 
	      $start += $self->{SUPERLINK}->{"$slink"}->{"$clone"}->[0] - 1; # clone base coords

	      # length is entire obj length unless specified
	      $length = 1 + $self->{SUPERLINK}->{"$slink"}->{"$clone"}->[1] - $self->{SUPERLINK}->{"$slink"}->{"$clone"}->[0] unless $length;
	      last SLINKS;
	    } 
	  }
	}
    }

    if( $chrom ) {
      my $seqlength = $self->{LENGTH}->{$chrom};
      if ($start >= $seqlength) {
	carp "subsequence invalid, past end of sequence : seq = $seq\tsequence length = $seqlength\tstart = $start\tlength = $length\textend = $extend\n";
	return;
      }
      $length = $seqlength unless $length; #full sequence of object.
      return $self->{SEQUENCE}->{"$chrom"} unless ($start and $length);#return whole seq if no coords passed
      $subseq = substr( ($self->{SEQUENCE}->{"$chrom"} ),$start, $length+(2*$extend) ); #extend either end
    } else {
      carp "couldn't work out chromosome that sequence $seq derives from\n";
    }

    if ($subseq ) {
      return $subseq;
    } else {
      carp "subsequence invalid : seq = $seq\tstart = $start\tlength = $length\textend = $extend\n";
    }
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
  Usage   :   my $comp = $seq_obj->($seq)
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
    return ($seq);
  }



1;
