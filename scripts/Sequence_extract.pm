package Sequence_extract;

use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"
  : glob("~ar2/wormbase/scripts");

use Carp;
use Wormbase;
use Coords_converter;

@ISA = ('Coords_converter');

sub invoke 
  {
    my $class = shift;
    my $database = shift;
    my $refresh = shift;

    # inherit from Coords_converter to get all the coord info
    my $self = Coords_converter->invoke($database,$refresh);

    bless $self, $class;
    print "now doing $class (Seq_ex) constructor\n\n\n";

    # get the chromosomal sequences
    my $tace = &tace;
    my @chromosome = qw( I II III IV V X );
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
      print "reading chromosome $_\n";
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

    carp "no length given for sequence to extract - using full length of $seq\n" unless (defined $length);

    #passed seq is a SUPERLINK
    if( $seq =~ /SUPERLINK/ ) {
      my $sl = $seq;
      $chrom = $self->_getChromFromSlink("$seq");

      # modify the starting coordinate
      $start += $self->{"$chrom"}->{SUPERLINK}->{"$sl"}->[0] - 1;
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
	    $start += $self->{"$chrom"}->{SUPERLINK}->{$slink}->[0] -1 ; # superlink base coords
	    $start += $self->{SUPERLINK}->{"$slink"}->{"$clone"}->[0] - 1; # clone base coords
	    
	    # length is entire obj length unless specified
	    $length = $self->{SUPERLINK}->{"$slink"}->{"$clone"}->[1] - $self->{SUPERLINK}->{"$slink"}->{"$clone"}->[0] unless $length;
	    last SLINKS;
	  } 
	}
      }
    }
    else {
      $start--; # chromosome coords start at 1, substr assumes 0 for 1st char.
    }

    $length = length($self->{SEQUENCE}->{"$chrom"}) unless $length; #full sequence of object.
    $subseq = substr( ($self->{SEQUENCE}->{"$chrom"} ),$start, $length+$extend );
  }

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
