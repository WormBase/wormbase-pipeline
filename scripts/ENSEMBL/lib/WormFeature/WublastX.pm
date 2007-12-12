#===============================================================================
#         FILE:  Wormy.pm
#  DESCRIPTION:  utility stuff to go into WormBase.pm of EnsEMBL
#        FILES:  Wormy.pm, wormy.t (for the tests)
#         BUGS:  YAML depends on a recent PERL version, but EnsEMBL on an old one, so I axed YAML
#        NOTES:  ---
#       Author:  (michael han), <mh6@sanger.ac.uk>
#      COMPANY:
#      VERSION:  1.0
#      CREATED:  16/05/06 15:09:49 BST
#     REVISION:  $Revision: 1.1 $
#===============================================================================

package WublastX;
use strict;
use IO::File;
use Bio::EnsEMBL::DnaPepAlignFeature;

# store things here
sub new {
    my ( $class, $slice, $file, $analysis,$type ) = @_;
    my $self = {
        slice    => $slice,
        analysis => $analysis,
        file     => $file,
	type     => $type
    };
    bless $self, $class;
    $self->_get_lines;
    return $self;
}

# parse things
sub _get_lines {
    my ($self) = @_;
    $self->{hits} = [];
    my $file     = $self->{file};
    my $gff_file = new IO::File "< $file";
    while ( my $line = <$gff_file> ) {
        next if !( $line =~ /wublastx\sprotein_match/ ); # wublastx specific
        my @col = split /\s/, $line;

        my $strand = $col[6] eq '+' ? 1 : -1;            # convert to ensembl style
	my $id=$col[9];
	$id =~s/Protein:(.+)://;
	my $type=$1;

	next if $type ne $self->{type};

	$id =~s/\"//g;
	my $length=$col[4]-$col[3];
        my ( $hstart, $hstop ) = sort { $a <=> $b } ( $col[10], $col[11] );
	if ($hstop-$hstart != $length){$hstop=$hstart+$length} # le fake

	push @{ $self->{hits} },
          Bio::EnsEMBL::DnaPepAlignFeature->new(
            -start         => $col[3],
            -end           => $col[4],
            -strand        => $strand,
            -slice         => $self->{slice},
            -analysis      => $self->{analysis},
            -score         => $col[5],
	    -p_value       => 1, # fake value because someone bound a double to this
	    -percent_id    => 100, # fake value because someone bound a int to this
	    -hseqname      => $id,
	    -hstart        => $hstart,
	    -hend          => $hstop, # sometimes faked to pass a insanity check
	    -cigar_string  => "${length}MM"
          );
#        push @{ $self->{hits} },
#          Bio::EnsEMBL::SimpleFeature->new(
#            -start         => $col[3],
#            -end           => $col[4],
#            -strand        => $strand,
#            -slice         => $self->{slice},
#            -analysis      => $self->{analysis},
#            -score         => $col[5],
#            -display_label => join( " ", @col[ 9 .. 11 ] )
#          );
    }
    return 1;
}

# lifted from Bronwen
sub save {
    my ( $self, $db ) = @_;
    my @pairs = @{ $self->{hits} };
    eval { print "\n check 1: " . $pairs[0]->display_label };
    eval { print "\n check 2: " . $pairs[0]->start . " - " . $pairs[0]->end . "\n" };

    my $pair_adaptor = $db->get_ProteinAlignFeatureAdaptor;

    eval{ $pair_adaptor->store(@pairs)};
    if ($@) {
        die "couldn't store features , problems: " . $@;
    }
    return 1;
}

1;
