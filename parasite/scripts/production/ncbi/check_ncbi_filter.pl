#!/usr/bin/env perl

=pod

=head1 Filter check_ncbi.sh output

Filters output of the C<check_ncbi.sh> script.  Accepts options matching the 
major sections (NCBI_ONLY, MATCHING_METADTA, etc.) of the output (not case 
sensitive).

Provide the C<check_ncbi.sh> output on STDIN, or the output file name as command 
line argument.

=head1 SYNOPSIS

  check_ncbi_filter.pl [options] [<file>]

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

-csv
   print CSV instead of YAML dump

-tsv
   print TSV instead of YAML dump

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use YAML;

use constant ASSEMBLY_KEYS => [qw(assembly_accession provider_name report_url submission_date assembly_url assembly_name)];

my (  $ncbi_only, $parasite_only, $match_not_intended, $matching_metadata, $matching_scaffold_names, $matching_scaffold_lengths, $mismatched,
      $invert, $csv, $tsv,
      $help);
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
            'csv'                         => \$csv,
            'tsv'                         => \$tsv,
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


my $data = YAML::Load(join('',<>));

# remove unwanted headings
for my $this_heading (keys %{$data}) {
   delete $data->{$this_heading} unless ( $invert ? ! $wanted{$this_heading} : $wanted{$this_heading} );
}

if($csv || $tsv) {
   for my $this_heading (sort keys %{$data}) {
      foreach my $this_species (sort keys %{$data->{$this_heading}}) {
         foreach my $this_bioproject (sort keys %{$data->{$this_heading}->{$this_species}}) {
            my $this_assemblies_array = $data->{$this_heading}->{$this_species}->{$this_bioproject};
            unless( ref([]) eq ref($this_assemblies_array) ) {
               require Data::Dumper;
               print STDERR "$this_species $this_bioproject doesn't contain an array (of assemblies)\n";
               print STDERR Dumper $this_assemblies_array;
               die "baled";
            }
            foreach my $this_assembly ( @{$this_assemblies_array} ) {
               unless( ref({}) eq ref($this_assembly) ) {
                  require Data::Dumper;
                  print STDERR "$this_species $this_bioproject has an asssemly which isn't a hash\n";
                  print STDERR Dumper $this_assembly;
                  die "baled";
               }
               print join( ($tsv ? "\t" : ','),
                           $this_heading,
                           $this_species,
                           $this_bioproject,
                           map((defined $this_assembly->{$_} ? $this_assembly->{$_} : ''), @{+ASSEMBLY_KEYS})
                           ),
                     "\n";
            }
         }
      }
   }
} else {
   print YAML::Dump( $data );
}
