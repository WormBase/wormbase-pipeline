#!/usr/bin/env perl
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
 
=head1 NAME
 
dump_hal_genome_ids.pl
 
=head1 DESCRIPTION
 
Dumps a mapping of genome IDs from a HAL file.
 
=head1 SYNOPSIS
    perl dump_hal_genome_sequence_names.pl --hal_file <hal_file_path>
 
=head1 OPTIONS
 
=over
 
=item B<[--help]>
 
Prints help message and exits.
 
=item B<[--hal_file path]>
 
Input HAL alignment file.
 
=cut
 
 
use strict;
use warnings;
 
use Getopt::Long;
use JSON;
use Pod::Usage;
 
use Bio::EnsEMBL::Compara::HAL::HALXS::HALAdaptor;
 
 
my ( $help, $hal_file);
GetOptions(
    "help|?"     => \$help,
    "hal_file=s" => \$hal_file,
) or pod2usage(-verbose => 2);
 
# Handle "print usage" scenarios
pod2usage(-exitvalue => 0, -verbose => 1) if $help;
pod2usage(-verbose => 1) if !$hal_file;
 
my $hal_adaptor = Bio::EnsEMBL::Compara::HAL::HALXS::HALAdaptor->new($hal_file);
 
my @hal_genomes = $hal_adaptor->genomes();
my @rel_genomes = grep { $_!~/^Anc[0-9]+$/ } @hal_genomes;
 
my $hal_genome_info;
foreach my $rel_genome (@rel_genomes) {
  print $rel_genome . "\n";
}
