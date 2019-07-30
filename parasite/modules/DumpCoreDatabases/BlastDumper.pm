package DumpCoreDatabases::BlastDumper;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
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
  my $out_dir = $self->param_required('out_dir');
  my $suffix = $self->param_required('suffix');

#connect to core database and get info
  my $mc = Bio::EnsEMBL::Registry->get_adaptor($species, 'Core', 'MetaContainer');
  my $host     = $mc->dbc->host();
  my $port     = $mc->dbc->port();
  my $prod_name = ucfirst( $mc->single_value_by_key('species.production_name') );
  my $assembly = $mc->single_value_by_key('assembly.default');
  my $dbname = $mc->dbc->dbname();
  $mc->dbc->disconnect_if_idle();

#define file name
  make_path($out_dir);
  my $prefix = "$prod_name.$assembly";
  $prefix =~ s/\s+/_/g;
  my $out_file = "$out_dir/$prefix.$suffix";

  
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


}


1;
