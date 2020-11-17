package ProteinConsequence::BaseProteinConsequence;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::AnalysisJob;

use base qw(Bio::EnsEMBL::Hive::Process);

sub param {
    my $self = shift;
    
    unless ($self->input_job) {
        # if we don't have an input job, add a dummy one (used when we're not 
        # running as part of a pipeline proper)
        $self->input_job(Bio::EnsEMBL::Hive::AnalysisJob->new);
    }

    return $self->SUPER::param(@_);
}


sub required_param {
    my $self        = shift;
    my $param_name  = shift;
    
    my $param_value = $self->param($param_name, @_);
    
    die "$param_name is a required parameter" unless defined $param_value;
    
    return $param_value;
}


### Logging methods to overwrite standard WormBase logging methods
sub write_to {
    my ($self, $msg) = @_;

    $self->warning($msg);

    return;
}


sub log_and_die {
    my ($self, $msg) = @_;

    die $msg;
}


sub mail {

}


1;
