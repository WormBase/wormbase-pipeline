=head1 NAME

ProteinConsequence::BaseProteinConsequence - base class for ProteinConsequence

=cut

=head1 DESCRIPTION

Base class for WormBase build VEP pipeline.

=cut

package ProteinConsequence::BaseProteinConsequence;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::AnalysisJob;

use base qw(Bio::EnsEMBL::Hive::Process);


=head2 param
    Title:    param
    Function: getter/setter for pipeline parameters
=cut

sub param {
    my $self = shift;
    
    unless ($self->input_job) {
        # if we don't have an input job, add a dummy one (used when we're not 
        # running as part of a pipeline proper)
        $self->input_job(Bio::EnsEMBL::Hive::AnalysisJob->new);
    }

    return $self->SUPER::param(@_);
}


=head2 required_param
 
    Title:    required_param
    Function: retrieves a pipelne parameter or dies if parameter does not exist
    Args:     parameter name
    Returns:  parameter value

=cut

sub required_param {
    my $self        = shift;
    my $param_name  = shift;
    
    my $param_value = $self->param($param_name, @_);
    
    die "$param_name is a required parameter" unless defined $param_value;
    
    return $param_value;
}


=head2 remove_header

    Title:    remove_header
    Function: remove header lines from VEP output files
    Args:     none
    Return:   n/a

=cut

sub remove_header {
    my $self = shift;

    my $file = $self->required_param('vep_input_file') . '.vep.txt';
    my @cmds = ("grep -v '^#' $file > $file.tmp", "mv $file.tmp $file");
    for my $cmd (@cmds) {
	my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
	die "Couldn't remove header from $file: $exit_code: $stderr" unless $exit_code == 0;
    }

    return;
}


###
### Logging methods to replace standard WormBase logging methods
###

=head2 write_to
    
    Title:    write_to
    Function: writes message to eHive database log table
    Args:     message string
    Returns:  n/a

=cut

sub write_to {
    my ($self, $msg) = @_;

    $self->warning($msg);

    return;
}


=head2 write_to
    
    Title:    log_and_dir
    Function: writes message to eHive database log table and kills pipeline
    Args:     message string
    Returns:  n/a

=cut

sub log_and_die {
    my ($self, $msg) = @_;

    die $msg;
}


=head2 mail

    Title:    mail
    Function: dummy method to prevent errors when code expecting Log_file object
    Args:     none
    Returns:  n/a

=cut

sub mail {
    return;
}


1;
