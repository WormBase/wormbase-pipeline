#===============================================================================
#         FILE:  Blat.pm
#  DESCRIPTION:  utility stuff to go into WormBase.pm of EnsEMBL
#        FILES:  
#         BUGS:  
#        NOTES:  ---
#       Author:  (michael han), <mh6@sanger.ac.uk>
#      COMPANY:
#      VERSION:  1.0
#      CREATED:  16/05/06 15:09:49 BST
#     REVISION:  $Revision: 1.1 $
#===============================================================================

package Blat;
use strict;
use IO::File;
use Error;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::PepDnaAlignFeature;

# store things here
sub new {
    my ( $class, $slice, $file, $analysis,$gff) = @_;
    throw  Error::Simple("Bad_Number_of_Arguments") if ! ref($gff);
    my $self = {
        slice    => $slice,
        analysis => $analysis,
        file     => $file,
	feature  => $gff->{-feature},
	source   => $gff->{-source},
	translated     => $gff->{-translated}
    };
    bless $self, $class;
    $self->_get_lines;
    return $self;
}

# parse the GFF and create features
sub _get_lines {
    my ($self) = @_;
    $self->{hits} = [];
    my $file     = $self->{file};
    my $gff_file = new IO::File "< $file";
    while ( my $line = <$gff_file> ) {
	my $feature = $self->{feature};
	my $source = $self->{source};

    	next if !( $line =~ /$feature\s$source/ );
        my @col = split /\s/, $line;

        my $strand = $col[6] eq '+' ? 1 : -1;    #convert to ensembl style
	my $id=$col[9];
	$id =~s/Protein:|Sequence://;
	$id =~s/\"//g;
	my $length=$col[4]-$col[3];
        my ( $hstart, $hstop ) = sort { $a <=> $b } ( $col[10], $col[11] );
	if ($hstop-$hstart != $length){$hstop=$hstart+$length} #le fake
	
	my $score= $col[5] eq '.' ? 100 : $col[5];

	# spot the difference between the Protein and DNA constructors ... ;-)
	# needs maybe some cleanup love
	if ($self->{translated}) {
	push @{ $self->{hits} },
          Bio::EnsEMBL::DnaPepAlignFeature->new(
            -start         => $col[3],
            -end           => $col[4],
            -strand        => $strand,
            -slice         => $self->{slice},
            -analysis      => $self->{analysis},
            -score         => $score,
	    -p_value       => 1, # fake value because someone bound a double to this
	    -percent_id    => 100, # fake value because someone bound a int to this
	    -hseqname      => $id,
	    -hstart        => $hstart,
	    -hend          => $hstop, # sometimes faked to pass a insanity check
	    -cigar_string  => "${length}MM"
          )
  	}
	else {
	push @{ $self->{hits} },
          Bio::EnsEMBL::DnaDnaAlignFeature->new(
            -start         => $col[3],
            -end           => $col[4],
            -strand        => 1,
            -slice         => $self->{slice},
            -analysis      => $self->{analysis},
            -score         => $score,
	    -p_value       => 1, # fake value because someone bound a double to this
	    -percent_id    => 100, # fake value because someone bound a int to this
	    -hseqname      => $id,
	    -hstart        => $hstart,
	    -hend          => $hstop, # sometimes faked to pass a insanity check
	    -hstrand       => $strand,
	    -cigar_string  => "${length}MM"
          )
	}
    }
    return 1;
}

# store features in db
# lifted from Bronwen
sub save {
    my ( $self, $db ) = @_;
    my @pairs = @{ $self->{hits} };
    eval { print "\n check 1: " . $pairs[0]->display_label };
    eval { print "\n check 2: " . $pairs[0]->start . " - " . $pairs[0]->end . "\n" };

    my $pair_adaptor = $self->{translated} ? $db->get_ProteinAlignFeatureAdaptor : $db->get_DnaAlignFeatureAdaptor;

    eval{ $pair_adaptor->store(@pairs)};
    if ($@) {
        die "couldn't store features , problems: " . $@;
    }
    return 1;
}

1;
