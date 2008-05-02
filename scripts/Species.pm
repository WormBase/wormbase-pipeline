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
# pulls the names from a fasta file
sub chromosome_names {
	my $self=shift;

	unless ($self->{'chromosome_names'}) {
		use IO::File;
		my $species = lc(ref($self));
		my $inf= new IO::File $self->wormpub."/DATABASES/$species/CHROMOSOMES/supercontigs.fa",'r' || croak(@$); #reference file
		my $prefix=$self->chromosome_prefix;
		my @ids;
		while(<$inf>){ push @ids,"$1" if $_=~/^>$prefix(\S+)/}
		$inf->close;
		$self->{'chromosome_names'}=\@ids;
	}
	return @{$self->{'chromosome_names'}}
}

sub mt_name {return undef}
sub chromosome_prefix {'PREFIX_'}
sub pep_prefix {return undef}
sub wormpep_prefix {return undef}

########################################
#
# Core species
#
#########################################

package Elegans;

our @ISA = qw(Wormbase Species);

sub chromosome_prefix {'CHROMOSOME_'}
sub chromosome_names {qw(I II III IV V X)}
sub mt_name {'MtDNA'}
sub pep_prefix {'CE'}
sub pepdir_prefix{'worm'};
sub cds_regex{qr/^[A-Z0-9_cel]+\.[1-9]\d?\d?[A-Za-z]?$/};  #the cel is for telomeric clone CDSs cTel54X.1
sub ncbi_tax_id {'6239'};
sub full_name {
	my $self = shift;
	my %param = @_ ;
	if($param{'-short'}){
		return 'C. elegans';
	}
	else { return'Caenorhabditis elegans'
	};
}
sub wormpep_prefix{'WP'}

########################################
package Briggsae;
our @ISA = qw(Wormbase Species);

sub _new {
    my $class = shift;
    my %param = %{ shift(@_) };
    my $self = $class->initialize( $class->flatten_params( \%param ) );

    # stuff post object creation goes here

    # overriding wormpep directory with brigpep
    $self->{'wormpep'}  = $self->brigpep;

    bless $self, $class;
}

sub chromosome_prefix {'chr'}
#sub chromosome_names {qw(I I_random II II_random III III_random IV IV_random V V_random X X_random Un)} # CB1
sub chromosome_names {qw(I I_random II II_random III III_random IV IV_random V V_random X Un)} # CB3

sub pep_prefix {'CBP'}
sub pepdir_prefix{'brig'};
sub cds_regex{qr/^CBG\d{5}$/};
sub ncbi_tax_id {'6238'};

sub full_name {
	my $self = shift;
	my %param = @_ ;
	if($param{'-short'}){
		return 'C. briggsae';
	}
	else { return'Caenorhabditis briggsae'
	};
}
sub wormpep_prefix{'BP'}

sub wormpep_files {
  my $self = shift;
  return ( "brigpep", "brigpep.accession", "brigpep.dna", "brigpep.history", "brigpep.fasta", "brigpep.table",
	   "brigpep.diff" );
}


#########################################
package Remanei;
use Carp;
our @ISA = qw(Wormbase Species);

sub _new {
    my $class = shift;
    my %param = %{ shift(@_) };
    my $self = $class->initialize( $class->flatten_params( \%param ) );

    # add stuff post object creation goes here

    bless $self, $class;
}
sub full_name {
	my $self = shift;
	my %param = @_ ;
	if($param{'-short'}){
		return 'C. remanei';
	}
	else { return'Caenorhabditis remanei'
	};
}

sub chromosome_prefix {'Crem_Contig'}
sub chromosome_names {
	my @supercontigs;
	my $i = 0;
	while($i < 9660){
		push (@supercontigs, $i) unless ($i == 5540);
		$i++;
	}
	return @supercontigs;
}

sub wormpep_prefix {'RP'}
sub pep_prefix {'RP'}
sub pepdir_prefix{'rema'};
sub ncbi_tax_id {'31234'};
sub cds_regex{qr/cr01\.sctg\d+\.wum\.\d+\.\d+/};

#######################################################
package Brenneri;
use Carp;
our @ISA = qw(Wormbase Species);

sub _new {	
    my $class = shift;
    my %param = %{ shift(@_) };

    # additional parameters
    #$param{'-autoace'}='/nfs/disk100/wormpub/TEST_BUILD/autoace/brenneri';

    my $self = $class->initialize( $class->flatten_params( \%param ) );

    # add stuff post object creation goes here

    bless $self, $class;
}

sub chromosome_prefix {'Contig'}
sub pep_prefix {'CN'}
sub pepdir_prefix{'bre'};
sub ncbi_tax_id {'135651'};
sub full_name {
	my $self = shift;
	my %param = @_ ;
	if($param{'-short'}){
		return 'C. brenneri';
	}
	else { return'Caenorhabditis brenneri'
	};
}

#######################################################

package Japonica;
use Carp;
our @ISA = qw(Wormbase Species);

sub _new {	
    my $class = shift;
    my %param = %{ shift(@_) };

    my $self = $class->initialize( $class->flatten_params( \%param ) );

    # add stuff post object creation goes here

    bless $self, $class;
}

sub chromosome_prefix {'Contig'}
sub pep_prefix {'JA'}
sub pepdir_prefix{'jap'};
sub ncbi_tax_id {'281687'};
sub full_name {
	my $self = shift;
	my %param = @_ ;
	if($param{'-short'}){
		return 'C. japonica';
	}
	else { return'Caenorhabditis japonica'
	};
}

#######################################################

package Pristionchus;
use Carp;
our @ISA = qw(Wormbase Species);

sub _new {
    my $class = shift;
    my %param = %{ shift(@_) };
    my $self = $class->initialize( $class->flatten_params( \%param ) );

    # add stuff post object creation goes here

    bless $self, $class;
}
sub full_name {
	my $self = shift;
	my %param = @_ ;
	if($param{'-short'}){
		return 'P. pacificus';
	}
	else { return'Pristionchus pacificus'
	};
}

sub chromosome_prefix {'Supercontig'}
sub chromosome_names {
	my @contigs;
	my $i = 0;
	while($i < 5060){
		push(@contigs, "$i");
		$i++;
	}
	return @contigs;
}

sub pep_prefix {'PP'}
sub pepdir_prefix{'ppa'};
sub ncbi_tax_id {'54126'};

#######################################################

package Heterorhabditis;
use Carp;
our @ISA = qw(Wormbase Species);

sub _new {
	my $class = shift;
    my %param = %{ shift(@_) };
    my $self = $class->initialize( $class->flatten_params( \%param ) );

    # add stuff post object creation goes here

    bless $self, $class;
}
sub full_name {
	my $self = shift;
	my %param = @_ ;
	if($param{'-short'}){
		return 'H. bacteriophora';
	}
	else { return'Heterorhabditis bacteriophora'
	};
}
sub ncbi_tax_id {'37862'};
sub chromosome_prefix {'Contig'}
sub pep_prefix {'HB'}
sub pepdir_prefix{'het'};

######################################################
#
# Tier3 species
#
######################################################

package Brugia;
use Carp;
our @ISA = qw(Wormbase Species);

sub _new {	
    my $class = shift;
    my %param = %{ shift(@_) };

    my $self = $class->initialize( $class->flatten_params( \%param ) );

    # add stuff post object creation goes here

    bless $self, $class;
}
sub full_name {
	my $self = shift;
	my %param = @_ ;
	if($param{'-short'}){
		return 'B. malayi';
	}
	else { return'Brugia malayi'
	};
}
sub chromosome_prefix {'Bmal_supercontig'}
sub pep_prefix {'BM'}
sub pepdir_prefix{'brug'};
sub ncbi_tax_id {'6279'};


#######################################################

#######################################################

if ( __FILE__ eq $0 ) {

    package Species;

    sub organism { my $self = shift; 'I am a ' . ref($self)."\n" }

    package Wormbase;

    sub worm { "I am a worm\n" }

    package Main;
    sub test {
    	my $worm = shift;
	print '-'x40,"\n";
	print $worm->organism();
	print $worm->full_name(),"\n";
	print $worm->worm();
	print 'my compost heap is at :' . $worm->autoace(), "\n";
	print $worm->wormpep_prefix(),"\n";
	print $worm->pep_prefix(),"\n";
	print $worm->ncbi_tax_id(),"\n";
	print $worm->chromosome_prefix(),"\n";
#	print "my chromosomes (custom prefix) are : ",join(' / ',$worm->get_chromosome_names('-prefix' => 'BLEEP_',-mito => 1)),"\n";
#	print "my chromosomes (w/o Mt + default prefix) are : ",join(' / ',$worm->get_chromosome_names('-prefix' => 1,)),"\n";
#	print "my chromosomes are : ",join(' / ',$worm->get_chromosome_names(-mito => 1)),"\n";
   	print '-'x40,"\n";
    }


    eval{my $a = Wormbase->new(-test => 1);&test($a)};
    print "$@\n" if ($@);


    foreach my $orgs (qw(Elegans Briggsae Remanei Brenneri Japonica Heterorhabditis Pristionchus Brugia)){
	    eval{
		    my $a = Wormbase->new( '-organism' => $orgs,-test => 1);
		    &test($a);
		};
	    print "$@\n" if ($@);
    }

}
1;


__END__

=pod

=head1 NAME Species.pm

=head1 DESCRIPTION extends Wormbase.pm

Species.pm provides Elegans, Briggsae and Remanei classes, which are inherited from Wormbase
and Species. Species provides generics for the more specific child classes.

runs tests if not used as module.

=head1 Species

=head2 flatten_params($hashref)

flattens an hash into an array

=head2  get_chromosome_names([-prefix => 1 / PREFIX , -mito => 1])

returns list of chromosomes including mitochondria (-mito => 1),
default prefix (-prefix => 1) or a custom prefix (-prefix => 'blablub_')

uses mock methods if called from within Species

=head2 _new([params,...])

calls the Wormbase constructor and blesses into the new class.

=head2 mock_methods (should be overwritten in the child classes)

=over 3 

=item
chromosome_names

=item 
mt_name

=item 
chromosome_prefix

=back

=head1 Elegans

inherits from Species and Wormbase

just overwrites the mock methods from Species

=head1 Briggsae

inherits from Species and Wormbase

changes basedirectory to: $basedir/autoace/briggsae
does use the CB1 assembly names.

=head1 Remanei

inherits from Species and Wormbase

uses '/nfs/disk100/wormpub/DATABASES/remanei/supercontigs.fa'
to get the sequence names by stripping prefix='Contig' from them.

=cut 
