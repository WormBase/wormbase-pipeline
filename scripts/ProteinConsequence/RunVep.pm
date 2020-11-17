=head1 NAME

ProteinConsequence::MapVariation - run VEP

=cut

=head1 DESCRIPTION

Runs VEP on WormBase variations.

=cut

package ProteinConsequence::RunVep;

use strict;

use File::Copy;
use File::Path qw(make_path remove_tree);
use Data::Dumper;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

sub run {
    my $self = shift;
    
    $self->param('vep_failure', 1);
    my $batch_id = $self->required_param('batch_id');
    my $input_file = $self->required_param('output_dir') . "/$batch_id/$batch_id";
    $self->param('vep_input_file', $input_file);
    my $output_file = $input_file . '.vep.txt';
    
    my $species = $self->required_param('species');
    my $output_dir = $self->required_param('output_dir');
    my $fasta_file = "$output_dir/$species.fa.gz";
    my $gff_file = "$output_dir/$species.gff.gz";

    my $plugin_str = 'ProtFuncAnnot,mod=WB,pass=' .
	$self->required_param('password');

    # Run VEP
    $self->dbc->disconnect_when_inactive(1);
    my $vep_dir = $self->required_param('vep_dir');
    my $cmd = "$vep_dir/vep -i $input_file -gff $gff_file -fasta $fasta_file --hgvs --hgvsg -shift_hgvs=0 --symbol --distance 0 --output_file $output_file --numbers --force_overwrite --plugin $plugin_str";	
    
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
    
    $self->dbc->disconnect_when_inactive(0);

    if ($exit_code != 0) {
	$self->warning("ERROR: $flat_cmd - exit code: $exit_code: $stderr") if $exit_code != 0;
    }
    else {
	$self->remove_header() unless $input_file =~ /\w0+1$/;
	$self->param('vep_failure', 0);
	system("rm $input_file");
    }
    system("rm ${output_file}_summary.html");
}


sub write_output {
    my $self = shift;

    if ($self->param('vep_failure')) {
	$self->dataflow_output_id({batch_id => $self->param('batch_id')}, 4);
    }
    else {
	$self->dataflow_output_id({batch_id => $self->param('batch_id')}, 2);
    }
}
    
    
1;
