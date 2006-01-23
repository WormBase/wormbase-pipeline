#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  Physical_map.pl
#
#        USAGE:  ./Physical_map.pl
#
#  DESCRIPTION:
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:   (), <>
#      COMPANY:
#      VERSION:  1.0
#      CREATED:  23/11/05 11:52:00 GMT
#     REVISION:  $Revision: 1.1 $ 
#===============================================================================

package Main;

# specify -database as options if you dont like to use current_DB
sub run {
	my ($options)=@_;
	my $acedb = '~wormpub/DATABASES/current_DB';
#	my $acedb = '/wormsrv2/autoace/';
	$acedb=$options->{'-database'} if $options->{'-database'};
	my $mapper= Physical_mapper->new($acedb,@ARGV); #generate a new mapper based on the files
	foreach my $file(@ARGV) {
		open INF, $file;
		while (<INF>) {
			s/\"//g;
			my @a=split;
			if ($a[1] eq 'Allele' && $a[2] eq 'SNP') {
				my $chr=$a[0];
				my $pos=($a[3]+$a[4])/2;
				my $id=$a[-1];
				print $mapper->to_ace($id,$pos,$chr); # create a SNP Ace file with interpolated positions
			}
		}
		close INF;
	}
}


# class to generate a mapper on the GFF files for gene/gene which have a map position
package Physical_mapper;

use strict;
use warnings;

sub new {
    my ( $class, $acedb, @files) = @_;

    # hash->chromosomes->bp_coordinates->genetic_position
    my %map = %{Map_func::build($acedb,@files)};
    my %sorted_map;

    foreach my $key ( keys %map ) {
        @{ $sorted_map{$key} } = sort { $a <=> $b } keys %{ $map{$key} };
    }
    my $self = {
        'pmap' => \%map,
        'smap' => \%sorted_map
    };
    bless $self, $class;
    return $self;
}

sub map {
    my ( $self, $pos, $chr ) = @_;
    my $last = 0;
    my $next = 0;
    return undef if $pos<@{ $self->{'smap'}->{$chr} }[0];
    for ( my $i = 0 ; $i < scalar @{ $self->{'smap'}->{$chr} } ; $i++ ) {
        my $current  = @{ $self->{'smap'}->{$chr} }[$i];
        my $mlast    = $self->{'pmap'}->{$chr}->{$last};
        my $mcurrent = $self->{'pmap'}->{$chr}->{$current};

        $next = $current;
        if ( $pos <= $current && $pos >= $last ) {
            my $pdiff    = ( $pos - $last ) / ( $current - $last ); # might  be wrong prefix ->better now?
            my $mlast    = $self->{'pmap'}->{$chr}->{$last};
            my $mcurrent = $self->{'pmap'}->{$chr}->{$current};

            my $mpos = ( $mcurrent - $mlast ) * $pdiff + $mlast;
            return ($mpos);
        }
        else {
            $last = $current;
        }
    }


    # glorious correction routine for the ends
    # should be f(x)=dx/dy + x1
    # my $x1 = $last
    # my $y1 = $self->{'pmap'}->{$chr}->{$last};
    # my $x2 = $next
    # my $y2 = $self->{'pmap'}->{$chr}->{$next};
#    my $f = sub { my $x = shift; return ( $x1 - $x2 ) / ( $y1 - $y2 ) + $x1 }

      return undef;
}

sub to_ace {
    my ( $self, $id, $map, $chr ) = @_;
    $chr =~ s/CHROMOSOME_//;
    my $mpos = $self->map( $map, $chr );
    if ($mpos) {
        return "Variation : \"$id\"\nInterpolated_map_position \"$chr\" $mpos \"Physical_map.pm\"\n\n";
    }
    else {
        return undef;
    }
}

# functions moved over from the old script to be used for generating the hashes
package Map_func;

#####################
# utility functions
#
sub extract_quotes {
	my ($line)=@_;
	return $1 if $line=~/\"(\S+)\"/;
	return undef;
}

sub read_hash {
	my ($file) = @_;
	my %hash;
	open INF,$file;
	while (<INF>) {
		my @a=split;
		$hash{$a[0]}=$a[1];
	}
	close INF;
	return \%hash;
}
                     #
######################

# grab from AceDB mapped genes
sub get_phys {
	my ($database)=@_;
	my %map;
	use Ace;
	my @genes;
	my $db = Ace->connect(-path  => $database);
	push @genes,$db->find('find Gene * where Map & SMap');
	foreach my $gene (@genes) {
		if (! defined $gene->Map(3)) {
			print STDERR "cannot find gene $gene ...\n";
			next;
		}
		my $name= "$gene";
		my $pos=  $gene->Map(3)->name;
		$map{$name}=$pos;
	}
		
	return \%map;
}

# build connection map
sub build {
	my ($acedb,@infiles)=@_;
	my %genes; # gene -> phys_map
	my %gen_map =%{get_phys($acedb)};
	foreach my $file(@infiles){
		open IN, $file;
		while (<IN>){
			s/\"//g;
			my @a=split;
			next if !($a[1] eq 'gene' && $a[2] eq 'gene');
			my $gene_id=$a[9];
			next if ! $gen_map{$gene_id};	

			my $chrom=$a[0]; my $map_pos=($a[3]+$a[4])/2; my $gen_pos=$gen_map{$gene_id};
			$chrom=~s/CHROMOSOME_//;
			#should be chromosome->map_pos->gen_pos
			$genes{$chrom}->{$map_pos}=$gen_pos;
		}
		close IN;
	}
	return \%genes;
}


1;
