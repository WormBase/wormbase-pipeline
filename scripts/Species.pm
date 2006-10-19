#!/usr/bin/env perl
#
# Demo for a species specific Class
#
# look at Wormbase::initialize for usage

use strict;
use Carp;

use Wormbase;

package Species;

sub flatten_params {
    shift;    # get rid of the class name
    my %param_hash = %{ shift() };
    my @param_array;
    while ( my ( $k, $v ) = each %param_hash ) { push @param_array, $k, $v }
    return @param_array;
}

# use with  (-prefix => 'whatever', -mito => 1)
sub get_chromosome_names {
	my $self= shift;
	my %options = @_;
	my @chromosomes=$self->chromosome_names;
	my $prefix=$options{'-prefix'};
	$prefix=$self->chromosome_prefix if ($prefix && length $prefix < 2);
	push @chromosomes, $self->mt_name if $options{'-mito'} && $self->mt_name;
	return map {$prefix.$_} @chromosomes if $options{'-prefix'};
	return @chromosomes
}

sub _new {
    my $class = shift;
    my %param = %{ shift(@_) };

    # additional parameters go here

    my $self = $class->initialize( $class->flatten_params( \%param ) );

    # add stuff here

    bless $self, $class;
}




# methods to overwrite
sub chromosome_names {return 'none'}
sub mt_name {return undef}
sub chromosome_prefix {'PREFIX_'}


########################################
package Elegans;

our @ISA = qw(Wormbase Species);

sub chromosome_prefix {'CHROMOSOME_'}
sub chromosome_names {qw(I II III IV V X)}
sub mt_name {'MtDNA'}


########################################
package Briggsae;
our @ISA = qw(Wormbase Species);

sub _new {
    my $class = shift;
    my %param = %{ shift(@_) };

    my $basedir='/nfs/disk100/wormpub';

    # additional parameters
    $param{'-autoace'}="$basedir/BUILD/autoace/briggsae";
    $param{'-autoace'}="$basedir/TEST_BUILD/autoace/briggsae" if ($param{'-test'});

    my $self = $class->initialize( $class->flatten_params( \%param ) );

    # stuff post object creation goes here

    bless $self, $class;
}

sub chromosome_prefix {'chr'}
sub chromosome_names {qw(I I_random II II_random III III_random IV IV_random V V_random X X_random Un)} # CB1

#########################################
package Remanei;
use Carp;
our @ISA = qw(Wormbase Species);

sub _new {
    carp('WARNING: Remanei is not yet fully integrated in the build');
	
    my $class = shift;
    my %param = %{ shift(@_) };

    # additional parameters
    $param{'-autoace'}='/nfs/disk100/wormpub/TEST_BUILD/autoace/remanei';

    my $self = $class->initialize( $class->flatten_params( \%param ) );

    # add stuff post object creation goes here

    bless $self, $class;
}

sub chromosome_prefix {'Contig'}

# pulls the names from a fasta file
sub chromosome_names {
	my $self=shift;

	unless ($self->{'chromosome_names'}) {
		use IO::File;

		my $inf= new IO::File '/nfs/disk100/wormpub/DATABASES/remanei/supercontigs.fa','r' || croak(@$); #reference file
		my $prefix=$self->chromosome_prefix;
		my @ids;
		while(<$inf>){ push @ids,"$1" if $_=~/>$prefix(.+)\n/}
		$inf->close;
		$self->{'chromosome_names'}=\@ids;
	}
	return @{$self->{'chromosome_names'}}
}

#######################################################

if ( __FILE__ eq $0 ) {

    package Species;

    sub organism { my $self = shift; print 'I am a C.' . ref($self) . "\n" }

    package Wormbase;

    sub worm { print "I am a worm\n" }

    package Main;
    sub test {
    	my $worm = shift;
	print '-'x40,"\n";
	$worm->organism;
	$worm->worm;
	print 'my compost heap is at :' . $worm->autoace(), "\n";
	print "my chromosomes (custom prefix) are : ",join(' / ',$worm->get_chromosome_names('-prefix' => 'BLEEP_',-mito => 1)),"\n";
	print "my chromosomes (w/o Mt + default prefix) are : ",join(' / ',$worm->get_chromosome_names('-prefix' => 1,)),"\n";
	print "my chromosomes are : ",join(' / ',$worm->get_chromosome_names(-mito => 1)),"\n";
    	print '-'x40,"\n";
    }


    eval{my $a = Wormbase->new(-test => 1);&test($a)};
    print "$@\n" if ($@);

    foreach my $orgs (qw(Elegans Briggsae Remanei TeddyBear)){
	    eval{
		    my $a = Wormbase->new( '-organism' => $orgs,-test => 1);
		    &test($a);
		};
	    print "$@\n" if ($@);
    }

}
1;

