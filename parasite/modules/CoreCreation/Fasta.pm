
# Take a genome fasta file, and creates a fasta file and agp
# that is suitable for loading into Ensembl
# Used to be a script, $WORM_CODE/scripts/ENSEMBL/scripts/split_genome_for_ensembl.pl
# Usage:
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
    $self->{toplevels} = \@toplevels;
    $self->{need_agp}  = $need_agp;
}

sub needs_contig_structure {
    my $self = shift;
    if ( not defined $self->{need_agp} ) {
        $self->_parse();
    }
    return $self->{need_agp};
}

sub split {
    my ( $self, $contig_fa_file, $agp_file ) = @_;
    if ( not defined $self->{toplevels} ) {
        $self->_parse();
    }
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

1;

