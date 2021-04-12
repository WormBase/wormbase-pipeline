=head1 NAME

ProteinConsequence::ProteinConsequence_conf - eHive config file

=cut

=head1 DESCRIPTION

eHive config file for WormBase Build variant mapping and VEP annotation pipeline.

=cut

package ProteinConsequence::ProteinConsequence_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub default_options {
    my ($self) = @_;

     return {

        # pipeline wide settings
        
	hive_root_dir                  => $ENV{'WBVEP_DIR'} . '/ensembl-hive',
        hive_force_init                => 1,
        hive_use_param_stack           => 0,
        hive_use_triggers              => 0,
        hive_auto_rebalance_semaphores => 0, 
        hive_no_init                   => 0,
        hive_default_max_retry_count   => 0,
        hive_debug_init                => 1,
	
	standard_max_workers           => 250,
	highmem_max_workers            => 100,
	hive_max_workers               => 350,
	
	# pipeline details
        pipeline_name           => 'wb_vep_' . lc($self->o('species')),
	pipeline_database       => $self->o('pipeline_database'),

	# folder locations
	vep_dir                 => $ENV{'WBVEP_DIR'} . '/ensembl-vep',
	pipeline_base_dir       => $ENV{'WBVEP_WORKING_DIR'},
        output_dir              => $self->o('pipeline_base_dir') . '/' . $self->o('pipeline_name'),  
        
        
        # connection details for the hive's own database
        pipeline_db => {
            -host   => $ENV{'WORM_DBHOST'},
            -port   => $ENV{'WORM_DBPORT'},
            -user   => $ENV{'WORM_DBUSER'},
            -pass   => $self->o('password'),            
            -dbname => $self->o('pipeline_database'),
            -driver => 'mysql',
        },

        # configuration for the various resource options used in the pipeline
        
	lsf_queue             => $ENV{'LSF_DEFAULT_QUEUE'},
        default_lsf_options   => '-q' . $self->o('lsf_queue') . ' -R"select[mem>2000] rusage[mem=2000]" -M2000',    
        highmem_lsf_options   => '-q' . $self->o('lsf_queue') . ' -R"select[mem>5000] rusage[mem=5000]" -M5000',     
        supermem_lsf_options  => '-q' . $self->o('lsf_queue') . ' -R"select[mem>16000] rusage[mem=16000]" -M16000',  

	# batch size
	batch_size              => 10000,

	# VEP parameters
	max_deletion_size       => 10000,

    };
}


sub resource_classes {
    my ($self) = @_;
    return {
	'default'  => { 'LSF' => $self->o('default_lsf_options')  },
	'highmem'  => { 'LSF' => $self->o('highmem_lsf_options')  },
	'supermem'  => { 'LSF' => $self->o('supermem_lsf_options')  },
    };
}


sub pipeline_analyses {
    my ($self) = @_;
    
    my @common_params = (
        debug      => $self->o('debug'),
	test       => $self->o('test'),
	output_dir => $self->o('output_dir'),
	species    => $self->o('species'),
    );

    return [	
        {   -logic_name => 'init_jobs',
            -module     => 'ProteinConsequence::InitJobs',
            -parameters => {
		@common_params,
		batch_size => $self->o('batch_size'),
		database   => $self->o('database'),		
		fasta      => $self->o('fasta'),
		gff_dir    => $self->o('gff_dir'),
            },
	    -input_ids => [{}],
            -rc_name    => 'highmem',
            -max_retry_count => 0,
            -flow_into  => {
                '2->A' => [ 'map_variation' ],
		'A->1' => [ 'collate_data' ],
		
            },
        },

        {   -logic_name     => 'map_variation',
            -module         => 'ProteinConsequence::MapVariation',
            -parameters     => {
                @common_params,
		database          => $self->o('database'),
		fasta             => $self->o('fasta'),
		max_deletion_size => $self->o('max_deletion_size'),
            },
            -failed_job_tolerance => 100,
            -max_retry_count => 0,
            -analysis_capacity  => $self->o('standard_max_workers'),
	    -hive_capacity      => $self->o('hive_max_workers'),
	    -rc_name        => 'highmem',
            -flow_into      => {
		3 => ['map_variation_highmem'],
		2 => ['run_vep'],
		-1 => ['map_variation_highmem'],
            },
        },

        {   -logic_name     => 'map_variation_highmem',
            -module         => 'ProteinConsequence::MapVariation',
            -parameters     => {
                @common_params,
		database          => $self->o('database'),
		fasta             => $self->o('fasta'),
		max_deletion_size => $self->o('max_deletion_size'),
            },
	    -analysis_capacity  => $self->o('highmem_max_workers'),
	    -hive_capacity      => $self->o('hive_max_workers'),
	    -rc_name        => 'supermem',
	    -flow_into      => {
		2 => ['run_vep'],
	    },
            -failed_job_tolerance => 0,
        },

	{   -logic_name     => 'run_vep',
            -module         => 'ProteinConsequence::RunVep',
            -parameters     => {
                @common_params,
		password        => $self->o('password'),
		vep_dir         => $self->o('vep_dir'),
            },
            -failed_job_tolerance => 100,
            -max_retry_count => 1,
            -analysis_capacity  => $self->o('standard_max_workers'),
	    -hive_capacity      => $self->o('hive_max_workers'),
	    -rc_name        => 'default',
            -flow_into      => {
		2 => ['create_ace'],
		4 => ['run_partial_vep'],
		-1 => ['run_partial_vep'],
            },
        },

	{   -logic_name     => 'run_partial_vep',
	    -module         => 'ProteinConsequence::RunPartialVep',
	    -parameters     => {
		password        => $self->o('password'),
		vep_dir         => $self->o('vep_dir'),
                @common_params,
	    },
	    -failed_job_tolerance => 100,
	    -max_retry_count => 0,
	    -analysis_capacity => $self->o('highmem_max_workers'),
	    -hive_capacity      => $self->o('hive_max_workers'),
	    -rc_name => 'highmem',
	    -flow_into => {
		2 => ['create_ace'],
		4 => ['run_partial_vep_supermem'],
		-1 => ['run_partial_vep_supermem'],
	    },
	},

	{   -logic_name     => 'run_partial_vep_supermem',
	    -module         => 'ProteinConsequence::RunPartialVep',
	    -parameters     => {
		password        => $self->o('password'),
		vep_dir         => $self->o('vep_dir'),
                @common_params,
	    },
	    -failed_job_tolerance => 0,
	    -max_retry_count => 3,
	    -analysis_capacity => $self->o('highmem_max_workers'),
	    -hive_capacity     => $self->o('hive_max_workers'),
	    -rc_name => 'supermem',
	    -flow_into => {
		2 => ['create_ace'],
	    }
	},

	{   -logic_name    => 'create_ace',
	    -module        => 'ProteinConsequence::CreateAce',
	    -parameters    => {
		@common_params,
		database          => $self->o('database'),
		max_deletion_size => $self->o('max_deletion_size'),
	    },
	    -failed_job_tolerance => 0,
	    -max_retry_count => 1,
	    -analysis_capacity => $self->o('standard_max_workers'),
	    -hive_capacity     => $self->o('hive_max_workers'),
	    -rc_name => 'highmem',
	},

	{   -logic_name    => 'collate_data',
	    -module        => 'ProteinConsequence::CollateData',
	    -parameters    => {
		@common_params,
		database          => $self->o('database'),
		pipeline_database => $self->o('pipeline_database'),
		pipeline_host     => $ENV{'WORM_DBHOST'},
		pipeline_port     => $ENV{'WORM_DBPORT'},
		pipeline_user     => $ENV{'WORM_DBUSER'},
		password          => $self->o('password'),
		log_dir           => $self->o('log_dir'),
		ace_dir           => $self->o('ace_dir'),
	    },
	    -failed_job_tolerance => 0,
	    -max_retry_count => 1,
	    -rc_name => 'default',
	},

		
    ];
}

1;

