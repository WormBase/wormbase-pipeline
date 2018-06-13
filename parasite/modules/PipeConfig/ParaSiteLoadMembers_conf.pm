package PipeConfig::ParaSiteLoadMembers_conf;

use strict;
use warnings;

use base('Bio::EnsEMBL::Compara::PipeConfig::EBI::EG::LoadMembers_conf');

sub default_options {
    my ($self) = @_;
    return {
        %{$self->SUPER::default_options},
        exclude_gene_analysis => {  'macrostomum_lignano_prjna284736' =>  ['mlignano_schatz_gene_bad']  },
        curr_core_sources_locs => [],
        prev_core_sources_locs => [],
        reuse_member_db => '',
    } ;
}
1;
