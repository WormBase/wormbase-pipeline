=pod 

=head1 NAME

Sequence_extract

=head1 SYNOPSIS

 my $seq_obj      = Sequence_extract->invoke($database, $refresh);
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

use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"
  : glob("~ar2/wormbase/scripts");

use Carp;
use Wormbase;
use Coords_converter;

@ISA = ('Coords_converter');


=head2 invoke

  Title   :   invoke
  Usage   :   my $extractor = Sequence_extract->invoke($database,1);
  Function:   Creates the object and loads in the data (generates fresh if requested)
  Returns :   ref to self
  Args    :   Database  (optional) - which database to use. Default is current_DB
              refresh   default is NULL - connect to the database and update coordinates
 Requires:    Database must have CHROMOSOMES directory with dna of each chromosome in CHROMOSOME_*.dna (FASTA format)

=cut


sub invoke 
  {
    my $class = shift;
    my $database = shift;
    my $refresh = shift;

    # inherit from Coords_converter to get all the coord info
    my $self = Coords_converter->invoke($database,$refresh);

    bless $self, $class;

    # get the chromosomal sequences
    my $tace = &tace;
    my @chromosome = qw( I II III IV V X MtDNA);
    my $seq_file = "$database/CHROMOSOMES/CHROMOSOME_I.dna";
    unless( -e "$seq_file" ) {
      open (ACE, "| $tace $database") or croak "cant connect to $database :$!\n";

      foreach my $chrom ( @chromosome ) {
	print "writing DNA seq for chromosome_$chrom\n";
	      my $command = <<EOF;
	clear
	  find sequence CHROMOSOME_$chrom
	    dna -f $database/CHROMOSOMES/CHROMOSOME_$chrom.dna
EOF
        print ACE $command;
      }
      close ACE;
    }
    foreach (@chromosome) {
      # read seq into $self
      $/ = "";
      open (SEQ, "$database/CHROMOSOMES/CHROMOSOME_$_.dna") or croak "cant open the dna file for $_:$!\n";
      my $seq = <SEQ>;
      close SEQ;
      $/ = "\n";
      $seq =~ s/>CHROMOSOME_\w+//;
      $seq =~ s/\n//g;
      $self->{'SEQUENCE'}->{"CHROMOSOME_$_"} = $seq;
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

    #passed seq is a SUPERLINK
    if( $seq =~ /SUPERLINK/ ) {
      my $sl = $seq;
      $chrom = $self->_getChromFromSlink("$seq");
      $seq = $chrom; # gets processed as chromosome below

      # modify the starting coordinate
      $start += $self->{"$chrom"}->{SUPERLINK}->{"$sl"}->[0];
      $length = $self->{"$chrom"}->{SUPERLINK}->{"$sl"}->[1] - $self->{"$chrom"}->{SUPERLINK}->{"$sl"}->[0] unless $length;
    }

    unless( $seq =~ /CHROMOSOME/ ) {
      # This is when the passed seq is a clone
    SLINKS:
      foreach my $slink (keys %{$self->{'SUPERLINK'}} ) {

	foreach my $clone (keys %{$self->{SUPERLINK}->{$slink}} ) {

	  if( "$clone" eq "$seq" ) {
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
    else {
      $chrom = $seq;
      $start--; # chromosome coords start at 1, substr assumes 0 for 1st char.
    }

    $length = length($self->{SEQUENCE}->{"$chrom"}) unless $length; #full sequence of object.
    $subseq = substr( ($self->{SEQUENCE}->{"$chrom"} ),$start, $length+(2*$extend) ); #extend either end

    if ($subseq ) {
      return $subseq;
    }
    else {
      carp "subsequence invalid : seq = $seq\tstart = $start\n\tlength = $length\n\textend = $extend\n";
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

1;
