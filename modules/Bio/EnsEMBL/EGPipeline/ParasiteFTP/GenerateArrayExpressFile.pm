=head1 LICENSE

Copyright [2009-2015] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package Bio::EnsEMBL::EGPipeline::ParasiteFTP::GenerateArrayExpressFile;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');
use base qw/Bio::EnsEMBL::Production::Pipeline::Base/;

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
  };
}

sub run {
  my ($self) = @_;
  my $out_dir = $self->param_required('out_dir');

  my $out_file = "$out_dir/assembly_names.txt";
  open(my $fh, '>', $out_file);

  my $division = "EnsemblParasite";

  my $all_dbas = Bio::EnsEMBL::Registry->get_all_DBAdaptors(-GROUP => 'core');

  foreach my $dba (@$all_dbas) {
	my $mc = $dba->get_MetaContainer();

	if ($division eq $mc->get_division()) {

		  my $assembly = $mc->single_value_by_key('assembly.name');
		  my $bp = $mc->single_value_by_key('species.ftp_genome_id');
		  my $taxon = $mc->single_value_by_key('species.taxonomy_id');
		  my $dbname = $mc->dbc->dbname();
  		  my ($sp) = ($dbname =~ /([^_]+_[^_]+)_.*?$/);

		  print $fh "$sp\t$bp\t$taxon\t$assembly\n";
		}
	}

  close $fh;

}

1;
