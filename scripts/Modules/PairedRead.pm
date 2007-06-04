=pod 

=head1 NAME

 PairedRead

=head1 SYNOPSIS


=head1 DESCRIPTION

 This object represents a pair of EST reads.

Inherits from SequenceObj ( SequenceObj.pm )

=head1 CONTACT

Anthony  ar2@sanger.ac.uk


=head1 METHODS

=cut


package PairedRead;
use lib $ENV{'CVS_DIR'} ;
use Carp;
use strict;

sub new {
	my $class = shift;
	my $self = {};
	my $read1 = shift;
	my $read2 = shift;
	
	unless ($read1 and $read2) {
		carp("missing reads in creation of ReadPair\n");
		return undef;
	}
    $self->{'reads'} = [($read1, $read2)];
	
    bless ( $self, $class );
    return $self;
}


sub reads {
	my $self = shift;
	return $self->{'reads'};
}

sub mapped {
	my $self = shift;
	my $map  = shift;
	$self->{'mapped'} = 1 if $map;
	return $self->{'mapped'};
}

1;