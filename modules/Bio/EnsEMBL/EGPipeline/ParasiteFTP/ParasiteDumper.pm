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

package Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Hive::Process');
use File::Path qw(make_path remove_tree);
use File::Copy;

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
  };
}

sub run {
  my ($self) = @_;

#get options
  my $script = $self->param_required('script');
  my $params = $self->param('params');
  my $species = $self->param_required('species');
  my $out_dir = $self->param_required('out_dir');
  my $ps_rel = $self->param_required('ps_rel');
  my $suffix = $self->param_required('suffix');

#connect to core database and get info
  my $mc       = $self->core_dba->get_MetaContainer();
  my $host     = $mc->dbc->host();
  my $port     = $mc->dbc->port();
  my $bioproject = $mc->single_value_by_key('species.ftp_genome_id');
  my $dbname = $mc->dbc->dbname();
  my ($sp) = ($dbname =~ /([^_]+_[^_]+)_.*?$/); 
  $mc->dbc->disconnect_if_idle();

#create directory structure
  my $dir = "$out_dir/WBPS$ps_rel/species/$sp/$bioproject";
  make_path($dir);

#define file name
  my $prefix = "$sp.$bioproject.WBPS$ps_rel";
  my $out_file = "$dir/$prefix.$suffix";

  $params .= " --host $host";
  $params .= " --port $port";
  $params .= " --user ensro";
  $params .= " -dbname $dbname";
  $params .= " -outfile $out_file";

# disconnect from Hive - the script might be running for a while
  $self -> dbc && $self -> dbc -> disconnect_if_idle();

#call dump script
  my $command = "perl $script $params";

  unless (system($command) == 0) {
  	$self->throw("Failed to execute script: '$command'.");
  }

#path to file flows into next analysis
  $self->dataflow_output_id({out_file => $out_file}, 4);

}


1;
