#!/usr/bin/env perl
#===============================================================================
#
#         FILE:  Map_Helper.pm
#
#  DESCRIPTION:  
#
#        FILES:  ---
#         BUGS:  ---
#        NOTES:  ---
#      $Author: mh6 $
#      COMPANY:  
#     $Version:  $
#      CREATED:  14/11/05 14:33:43 GMT
#        $Date: 2005-11-14 15:31:28 $
#===============================================================================
package Map_Helper;

use strict;
use warnings;


##########################################################
# gff input file helper moved over from mapping scripts
# will update $genes (hashref of Genes Exons)
# 
sub get_from_gff {
    my ( $file, $type, $unless, $genes ) = @_;
    open( GFF_in, "<$file" )
      || die "Failed to open gff file\n\n";
    while (<GFF_in>) {
        chomp;
        s/\#.*//;
        next unless /$unless/;
        my @f = split /\t/;

        my ($name) = ( $f[8] =~ /\"(\S+)\"/ );
        my $x = new Exon;
        $x->id($name);
        $x->start( $f[3] );
        $x->stop( $f[4] );
        $x->type($type);
        if ( $$genes{ $x->id } ) {
            push @{ $$genes{ $x->id }->exons }, $x;
            $$genes{ $x->id }->start( $x->start )
              if $x->start < $$genes{ $x->id }->start;
            $$genes{ $x->id }->stop( $x->stop )
              if $x->stop < $$genes{ $x->id }->stop;
        }
        else {
            $$genes{ $x->id } = new Gene;
            $$genes{ $x->id }->start( $x->start );
            $$genes{ $x->id }->stop( $x->stop );
        }
        #print "Gene : '$name' [$type]\n"
        #  if $verbose;
    }
    close(GFF_in);

}


###########################################
# helper function to search for overlaps
# will update %$output with Exons
# 
sub map_it {
    my ($output,$pcr,$sorted_genes,$genes)=@_;
    foreach my $testPCR ( keys %$pcr ) {
        my $PCRstart = $$pcr{$testPCR}->[0];
        my $PCRstop  = $$pcr{$testPCR}->[1];

        foreach my $gene (@$sorted_genes) {

            if ( $PCRstart > $$genes{$gene}->stop ) {
                next;
            }

            elsif ( $PCRstop < $$genes{$gene}->start ) {
                last;
            }

            else {
                my %added_exons;
                foreach my $exon ( @{ $$genes{$gene}->exons } ) {
                    if ( $PCRstart > $exon->stop 
                      ||  $PCRstop < $exon->start ) {
                        next
                          if $added_exons{ $exon->type . $exon->id };
                          push @{ $$output{$testPCR} }, $exon;
                          $added_exons{ $exon->type . $exon->id } = 1;
                      }
                }
            }
        }
    }
}


1;

__END__

=pod

=head2 NAME - Map_Helper.pm

=head2 USAGE

 use Map_Helper
 Map_Helper::Function()

=over 4

=head2 DESCRIPTION 

Helper module to be included in mapping scripts

=back

=head3 FUNCTIONS

=over 4

=item map_it($output,$pcr,$sorted_genes,$genes)

maps $pcr s to overlapping $genes.

$output is a hashref containing item -> matched exons

$pcr    is a hashref containing the to be queried arrays id->[start,stop]

$sorted_genes is a arrayref containing genenames sorted by start then stop

$genes is a hashref containing the genes/transcripts/whatever

=item  get_from_gff( $file, $type, $unless, $genes )

reads a gff into a hash of Genes containing Exons (%$genes), key is the Exon->id

struct( Exon => [ start => '$', stop => '$', type => '$', id => '$' ] );
struct( Gene => [ start => '$', stop => '$', exons => '@' ] );


=back

=head3 AUTHOR
$Author: mh6 $

=cut
