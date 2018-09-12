
# Take a genome fasta file, and creates a fasta file and agp
# that is suitable for loading into Ensembl
# Used to be a script, $WORM_CODE/scripts/ENSEMBL/scripts/split_genome_for_ensembl.pl
# Usage:
# perl -MCoreCreation::Fasta -e "print CoreCreation::Fasta->new($path)->mito"
# my $fasta = CoreCreation::Fasta->new($path);
# my $needs_split = $fasta -> needs_contig_structure();
# $fasta -> split("$path.split", $agp_path); 
use strict;

package CoreCreation::Fasta;
use Bio::SeqIO;

our $MAX_CONTIG_LEN = 10000000;

sub new {
    my ( $class, $genome_fa ) = @_;
    die "Required: <path to fasta with genome>" unless -f $genome_fa;
    my $self = {};
    $self->{genome_fa} = $genome_fa;
    bless $self, $class;
    return $self;
}
sub _parsed {
   my $self = shift;
   return defined $self->{need_agp};
}
sub _parse {
    my $self      = shift;
    my $genome_fa = $self->{genome_fa};
    my $seqio     = Bio::SeqIO->new(
        -format => 'fasta',
        -file   => ( $genome_fa =~ /\.gz$/ )
        ? "gunzip -c $genome_fa |"
        : $genome_fa
    );
    my @toplevels;
    my $need_agp = 0;
    while ( my $seq = $seqio->next_seq ) {
        if ( $seq->length > $MAX_CONTIG_LEN ) {
            $need_agp = 1;
        }
        push @toplevels, $seq;
    }
    Carp::croak("Fasta appears empty: $self->{genome_fa}") unless @toplevels;
    $self->{toplevels} = \@toplevels;
    $self->{need_agp}  = $need_agp;
}

sub needs_contig_structure {
    my $self = shift;
    $self->_parse() unless $self->_parsed;
    return $self->{need_agp};
}

sub split {
    my ( $self, %h) = @_;
    $self->_parse() unless $self->_parsed;
    my $contig_fa_file = $h{fasta};
    my $agp_file = $h{agp};
    my @toplevels = @{ $self->{toplevels} };

    my $count = 0;
    my $outseqio = Bio::SeqIO->new(
        -format => 'fasta',
        -file   => ">$contig_fa_file"
    );
    open( my $agp_fh, ">$agp_file" )
      or die "Could not open $agp_file for writing\n";

    foreach my $seq (@toplevels) {
        if ( $seq->length < $MAX_CONTIG_LEN ) {
            $outseqio->write_seq($seq);
            printf( $agp_fh "%s\t%d\t%d\t%d\tW\t%s\t%d\t%d\t+\n",
                $seq->id, 1, $seq->length, ++$count, $seq->id, 1,
                $seq->length );
        }
        else {
            my $seg_count = 1;
            for ( my $i = 0 ; $i < $seq->length ; $i += $MAX_CONTIG_LEN ) {
                my $start = $i + 1;
                my $end   = $start + $MAX_CONTIG_LEN - 1;
                $end = $seq->length if $end > $seq->length;

                my $new_id = sprintf( "%s.%d", $seq->id, $seg_count++ );

                my $seq_seg = Bio::PrimarySeq->new(
                    -id  => $new_id,
                    -seq => substr( $seq->seq, $start - 1, $end - $start + 1 )
                );
                $outseqio->write_seq($seq_seg);

                printf( $agp_fh "%s\t%d\t%d\t%d\tW\t%s\t%d\t%d\t+\n",
                    $seq->id, $start, $end, ++$count, $seq_seg->id, 1,
                    $seq_seg->length );
            }
        }
    }

}

sub mito {
    my $self = shift;
    $self->_parse() unless $self->_parsed;
    my @mito;
    for my $top_level (@{$self->{toplevels}}) {
       my $id = $top_level->display_id;
       push @mito , $id if $id =~ /mito/i or $id =~ /mtDNA/i;
    } 
    Carp::croak("Multiple mitochondrial scaffolds?: @mito") if @mito > 1;
    my $result = pop @mito;
    Carp::croak("Not actually right- comma at the end? $result") if $result =~ /,$/;
    return $result;
}
1;

