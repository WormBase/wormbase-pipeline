#!/usr/bin/env perl

=pod

=head1 Filter check_ncbi.sh output

Filters output of the C<check_ncbi.sh> script.  Accepts options matching the 
major sections (NCBI_ONLY, MATCHING_METADTA, etc.) of the output (not case 
sensitive).

Provide the C<check_ncbi.sh> output on STDIN, or the output file name as command 
line argument.

=head1 SYNOPSIS

  check_ncbi_filter.pl [-ncbi_only] [-parasite_only] [-match_not_intended] [-matching_metadata] [-matching_scaffold_names] [-matching_scaffold_lengths] [-mismatched] [<file>]

=head2 OPTIONS

-ncbi_only
   entries under NCBI_ONLY heading

-parasite_only
   entries under PARASITE_ONLY heading

-match_not_intended, -not_intended
   entries under MATCH_NOT_INTENDED heading

-matching_metadata, -metadata
   entries under MATCHING_METADATA heading

-matching_scaffold_names, -scaffold_names, -names
   entries under MATCHING_SCAFFOLD_NAMES heading

-matching_scaffold_lengths, -scaffold_lengths, -lengths
   entries under MATCHING_SCAFFOLD_LENGTHS heading

-mismatched
   entries under MISMATCHED heading

-v, --invert-match
   invert matching, i.e. display entries I<not> under the specified headings

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use YAML;

my ($ncbi_only, $parasite_only, $match_not_intended, $matching_metadata, $matching_scaffold_names, $matching_scaffold_lengths, $mismatched, $invert, $help);
GetOptions( 'ncbi_only'                   => \$ncbi_only,
            'parasite_only'               => \$parasite_only,
            'match_not_intended'          => \$match_not_intended,
                  'not_intended'          => \$match_not_intended,
            'matching_metadata'           => \$matching_metadata,
                     'metadata'           => \$matching_metadata,
            'matching_scaffold_names'     => \$matching_scaffold_names,
                     'scaffold_names'     => \$matching_scaffold_names,
                              'names'     => \$matching_scaffold_names,
            'matching_scaffold_lengths'   => \$matching_scaffold_lengths,
                     'scaffold_lengths'   => \$matching_scaffold_lengths,
                              'lengths'   => \$matching_scaffold_lengths,
            'mismatched'                  => \$mismatched,
            'v'                           => \$invert,
            'invert-match'                => \$invert,
            'help'                        => \$help
            )
            || pod2usage({-exitval=>1});
$help && pod2usage({-verbose=>2, -exitval=>0});

# hash providing lookup of the tags that are wanted
my %wanted = ( 'NCBI_ONLY'                   => $ncbi_only,
               'PARASITE_ONLY'               => $parasite_only,
               'MATCH_NOT_INTENDED'          => $match_not_intended,
               'MATCHING_METADATA'           => $matching_metadata,
               'MATCHING_SCAFFOLD_NAMES'     => $matching_scaffold_names,
               'MATCHING_SCAFFOLD_LENGTHS'   => $matching_scaffold_lengths,
               'MISMATCHED'                  => $mismatched
               );

my $tag_patterns  = '('.join(')|(',keys %wanted).')';
my $tag_regex     = qr/^($tag_patterns)\s*:\s*$/;

# $tag_regex = qr/^($tag_patterns)/;

my $output_on;
while(my $line = <>) {
   my($tag) = $line =~ $tag_regex;
   if(defined $tag) {
      # this line is a tag, so turn output on/off according to flag in %wanted
      $output_on = (defined $invert && $invert) ? !$wanted{$tag} : $wanted{$tag};
   }
   $output_on && print $line;
}

