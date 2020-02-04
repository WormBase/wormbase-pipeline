#!/usr/bin/env perl

=pod

=head1 SYNOPSIS

  doc_for_assembly.pl [accession0 ... accessionN]
  
=head1 DECRIPTION

Retrieves assembly report from NCBI and dump as YAML.

One or more NCBI assembly accession numbers must be provided either on standard input, or as arguments.

The YAML is the product of importing the NCBI report to a Perl data structure with XML::Simple::XMLin,
and then dumping that structure with YAML::Dump.

=head1 TO DO

=over

=item * Validation and verification

=back

=cut

use strict;
use warnings;
use YAML;
use XML::Simple;
use LWP::UserAgent;

my @assemblies;
@assemblies = split "\n", do {local $/; <>} unless -t STDIN;
@assemblies = @ARGV if @ARGV;
die unless @assemblies;
foreach my $this_assembly (@assemblies){
  # more granularity here to clarify cause of errors
  # using $_ makes code fragile
  # # my $doc = `curl -s "https://www.ncbi.nlm.nih.gov/assembly/$_?report=xml&format=text" `;
  # my $doc = `curl -s "https://www.ncbi.nlm.nih.gov/assembly/$this_assembly?report=xml&format=text" `;
  my $ncbi_url = "https://www.ncbi.nlm.nih.gov/assembly/$this_assembly?report=xml&format=text";
  
  my $ua = LWP::UserAgent->new() || die "failed to create LWP::UserAgent: $!";
  my $response = $ua->get($ncbi_url) || die "Error requesting $ncbi_url: $!";
  $response->is_success() || die "Request for $ncbi_url failed: ".$response->status_line();
  
  my $doc = $response->decoded_content();

  print Dump XMLin(XMLin($doc));
} 
