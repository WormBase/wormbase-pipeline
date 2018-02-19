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

package Bio::EnsEMBL::EGPipeline::ParasiteFTP::BlastDumper;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');
use File::Path qw(make_path);
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
  my $out_dir = $self->param_required('blast_dir');
  my $suffix = $self->param_required('suffix');

#connect to core database and get info
  my $mc       = $self->core_dba->get_MetaContainer();
  my $host     = $mc->dbc->host();
  my $port     = $mc->dbc->port();
  my $prod_name = ucfirst( $mc->single_value_by_key('species.production_name') );
  my $assembly = $mc->single_value_by_key('assembly.default');
  my $dbname = $mc->dbc->dbname();
  $mc->dbc->disconnect_if_idle(); 

#define file name
  make_path($out_dir);
  my $prefix = "$prod_name.$assembly";
  my $out_file = "$out_dir/$prefix.$suffix";

  
  $params .= " --host $host";
  $params .= " --port $port";
  $params .= " --user ensro";
  $params .= " -dbname $dbname";
  $params .= " -outfile $out_file";

#call dump script
  my $command = "perl $script $params";

  unless (system($command) == 0) {
  	$self->throw("Failed to execute script: '$command'.");
  }


}


1;
