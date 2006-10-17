#!/usr/bin/env perl
#
# Demo for a species specific Class
#
# look at Wormbase::initialize for usage

use Wormbase;

package Elegans;

@ISA=qw(Wormbase);

use strict;

sub _new {
	my $class=shift;
	my %param=%{shift(@_)};
	my @flat_param; # needed for Wormbase constructor
	while(my($k,$v)=each %param){
		push @flat_param,$k,$v;
	}
	my $self=$class->SUPER::new(@flat_param);

	# add stuff here

	bless $self,$class;
}


if (__FILE__ eq $0){
	package Elegans;

	sub organism {my $self=shift;print 'I am a C.'.ref($self)."\n"}

	package Wormbase;

	sub worm {print "I am a worm\n"}

	my $a=Wormbase->initialize('-organism' => 'Elegans');
	$a->worm;
	$a->organism;
	print 'my compost heap is at :'.$a->autoace() ,"\n";
	my $b=Wormbase->initialize('-organism' => 'TeddyBear');
}
1;


