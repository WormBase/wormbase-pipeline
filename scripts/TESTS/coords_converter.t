#!/usr/bin/env perl

use lib $ENV{CVS_DIR};

use Wormbase;
use Coords_converter;

use Test::More tests => 7;

# create a Test WormBase object
my $wormbase = Wormbase->new(-test => 1);

my $coords=Coords_converter->invoke($wormbase->autoace,1,$wormbase);

# 1 - test offset of the clones in the parent seq
    my $offset = $coords->CloneOffset( 'AH6' );
    is ($offset,9512569,'CloneOffset');

# 2 - test getting the clone from a chromosome position
    my $clone = $coords->GetCloneFromCoord('I',1000000);
    is ($clone,'C54G6','GetCloneFromCoord');

#test getting the span based on chromosome coordinates
    my @span = $coords->LocateSpan('CHROMOSOME_I',5239404,5271341);
    is ($span[0],'B0261','LocateSpan - clone name'); # 3
    is ($span[1],104,'LocateSpan - start');  # 4
    is ($span[2],32041,'LocateSpan - end');    # 5

#convert clone coordinates to chromosome ones
    my ($chrom,$coord) = $coords->Coords_2chrom_coords('AH6',2343);
    is ($chrom,'CHROMOSOME_II','Coords_2chrom_coords - chromosome name');   # 6
    is ($coord,9514911,'Coords_2chrom_coords - chromosome position'); # 7
