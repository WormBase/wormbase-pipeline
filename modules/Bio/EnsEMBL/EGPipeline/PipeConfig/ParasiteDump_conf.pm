
=pod

=head1 NAME
Bio::EnsEMBL::EGPipeline::PipeConfig::ParasiteDump_conf

=head1 DESCRIPTION

Configuration for running the ParaSite FTP dumping pipeline. Please see documentation on Confluence https://www.ebi.ac.uk/seqdb/confluence/display/WORMBASE/Parasite+Production+Pipelines

=head1 Author

Myriam Shafie

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::ParasiteDump_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::Version 2.3;
use base 'Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf';

sub default_options {
  my ($self) = @_;

  return {

    %{$self->SUPER::default_options},
    max_hive_capacity => 200,

    division => ['parasite'],
    dump_genome     => $ENV{WORM_CODE}."/scripts/ENSEMBL/scripts/dump_genome.pl",
    dump_transcript => $ENV{WORM_CODE}."/scripts/ENSEMBL/scripts/dump_transcripts.pl",
    dump_gff3       => $ENV{WORM_CODE}."/scripts/ENSEMBL/scripts/dump_gff3.pl",
    dump_gtf        => $ENV{WORM_CODE}."/scripts/ENSEMBL/scripts/dump_gtf_from_ensembl.pl",
    wbps_version    => $ENV{PARASITE_VERSION},
    

    dump_ftp        => 1,
    dump_blast      => 1,
    dump_vep        => 0,

    pipeline_db => {
      -host   => $self->o('host'),
      -port   => $self->o('port'),
      -user   => $self->o('user'),
      -pass   => $self->o('pass'),
      -dbname => $ENV{'USER'} . "_parasitedump_" .$self->o('wbps_version'),
      -driver => 'mysql',
    },


    #
    # ftp dumping options
    #
    ftp_species => [],
    # put core species in here
    ftp_antispecies => [],

    ftp_outdir => "/nfs/ftp/pub/databases/wormbase/staging/parasite/releases",
    ftp_genome	             => 1,
    ftp_genome_masked	     => 1,
    ftp_genome_softmasked    => 1,
    ftp_cds                  => 1,
    ftp_mrna                 => 1,
    ftp_protein              => 1,
    ftp_gff3                 => 1,
    ftp_gtf                  => 1,
    ftp_create_core_symlinks => 1,


    #
    # BLAST dumping options
    #
    blast_outdir => "/nfs/nobackup/ensemblgenomes/wormbase/parasite/ftp_dumps/vep",
    blast_species => [],
    blast_antispecies => [],

    blast_genome         => 1,
    blast_genome_masked  => 1,
    blast_mrna           => 1,
    blast_protein        => 1,
  };
}

# Force an automatic loading of the registry in all workers.
sub beekeeper_extra_cmdline_options {
  my $self = shift;

  return "-reg_conf ".$self->o("registry");
}

sub pipeline_analyses {
  my ($self) = @_;

  my $ftp_analysis = [];

  if ($self->o('ftp_genome')) {
    push @$ftp_analysis, 'DumpGenome';
  }
  if ($self->o('ftp_genome_masked')) {
    push @$ftp_analysis, 'DumpGenomeMasked';
  }
  if ($self->o('ftp_genome_softmasked')) {
    push @$ftp_analysis, 'DumpGenomeSoftMasked';
  }
  if ($self->o('ftp_cds')) {
    push @$ftp_analysis, 'DumpTranscriptsCDS';
  }
  if ($self->o('ftp_mrna')) {
    push @$ftp_analysis, 'DumpTranscriptsmRNA';
  }
  if ($self->o('ftp_protein')) {
    push @$ftp_analysis, 'DumpProteins';
  }
  if ($self->o('ftp_gff3')) {
    push @$ftp_analysis, 'DumpGFF3';
  }
  if ($self->o('ftp_gtf')) {
    push @$ftp_analysis, 'DumpGTF';
  }

  my $blast_analysis = [];

  if ($self->o('blast_genome')) {
    push @$blast_analysis, 'DumpGenomeBlast';
  }
  if ($self->o('blast_genome_masked')) {
    push @$blast_analysis, 'DumpGenomeMaskedBlast';
  }
  if ($self->o('blast_mrna')) {
    push @$blast_analysis, 'DumpTranscriptsBlast';
  }
  if ($self->o('blast_protein')) {
    push @$blast_analysis, 'DumpProteinsBlast';
  }

  my $ini = [];
  if ($self->o('dump_ftp')) {
    push @$ini, 'SpeciesFactory';
  }
  if ($self->o('dump_blast')) {
    push @$ini, 'SpeciesFactoryBlast';
  }

  return [
    {
      -logic_name      => 'InitialisePipeline',
      -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -max_retry_count => 0,
      -input_ids       => [{}],
      -flow_into       => {
        1 => $ini,
      },
      -meadow_type     => 'LOCAL',
    },
    #
    # BLAST
    # 
    {
      -logic_name        => 'SpeciesFactoryBlast',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -max_retry_count   => 1,
      -parameters        => {
                              division        => $self->o('division'),
                              species         => $self->o('blast_species'),
                              antispecies     => $self->o('blast_antispecies'),
                              meta_filters    => {},
                              chromosome_flow => 0,
                              variation_flow  => 0,
			      regulation_flow => 0,
			      core_flow       => 4, 
                            },
      -flow_into         => {
				4 => $blast_analysis,
			    },
      -meadow_type	 => 'LOCAL',
    },
    {
      -logic_name        => 'DumpGenomeBlast',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::BlastDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script    => $self->o('dump_genome'),
                              blast_dir => $self->o('blast_outdir'),
                              suffix    => "dna.toplevel.fa",
                              params    => "-ebiblastheader WBPS",
                            },
      -rc_name           => '2Gb_job',
    },
    {
      -logic_name        => 'DumpGenomeMaskedBlast',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::BlastDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script    => $self->o('dump_genome'),
                              blast_dir => $self->o('blast_outdir'),
                              suffix    => "dna_rm.toplevel.fa",
                              params    => "-mask -ebiblastheader WBPS",
                            },
      -rc_name           => '2Gb_job',
    },
    {
      -logic_name        => 'DumpTranscriptsBlast',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::BlastDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script    => $self->o('dump_transcript'),
                              blast_dir => $self->o('blast_outdir'),
                              suffix    => "cdna.all.fa",
                              params    => "-mrna -ebiblastheader WBPS",
                            },
      -rc_name           => '2Gb_job',
    },
    {
      -logic_name        => 'DumpProteinsBlast',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::BlastDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_transcript'),
                              blast_dir => $self->o('blast_outdir'),
                              suffix    => "pep.all.fa",
                              params    => "-pep -ebiblastheader WBPS",
                            },
      -rc_name           => '2Gb_job',
    },
    #
    # FTP
    #
    {
      -logic_name        => 'SpeciesFactory',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -max_retry_count   => 1,
      -parameters        => {
                              species         => $self->o('ftp_species'),
                              antispecies     => $self->o('ftp_antispecies'),
                              division        => $self->o('division'),
                              meta_filters    => {},
                              chromosome_flow => 0,
                              variation_flow  => 0,
			      regulation_flow => 0,
			      core_flow => 4, 
                            },
      -flow_into         => {
				4 => $ftp_analysis,
			    },
      -meadow_type	 => 'LOCAL',
    },
    {
      -logic_name        => 'DumpGenome',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_genome'),
                              out_dir    => $self->o('ftp_outdir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
			      suffix    => "genomic.fa",
                            },
      -flow_into         => { 
				4 => ['Gzip'],
				},
      -rc_name           => '2Gb_job',
    },
    {
      -logic_name        => 'DumpGenomeMasked',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_genome'),
                              out_dir    => $self->o('ftp_outdir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
                              suffix    => "genomic_masked.fa",
			      params    => "-mask",
                            },
      -flow_into         => { 
                                4 => ['Gzip'],
                                },
      -rc_name           => '2Gb_job',
    },
    {
      -logic_name        => 'DumpGenomeSoftMasked',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_genome'),
                              out_dir    => $self->o('ftp_outdir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
                              suffix    => "genomic_softmasked.fa",
                              params    => "-softmask",
                            },
      -flow_into         => { 
                                4 => ['Gzip'],
                                },
      -rc_name           => '2Gb_job',
    },
    {
      -logic_name        => 'DumpTranscriptsCDS',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_transcript'),
                              out_dir    => $self->o('ftp_outdir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
                              suffix    => "CDS_transcripts.fa",
                              params    => "-cds",
                            },
          -flow_into         => { 
                                4 => ['Gzip'],
                                },
      -rc_name           => '2Gb_job',
    },
    {
      -logic_name        => 'DumpTranscriptsmRNA',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_transcript'),
                              out_dir    => $self->o('ftp_outdir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
                              suffix    => "mRNA_transcripts.fa",
                              params    => "-mrna",
                            },
      -flow_into         => { 
                                4 => ['Gzip'],
                                },
      -rc_name           => '2Gb_job',
    },
    {
      -logic_name        => 'DumpProteins',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_transcript'),
                              out_dir    => $self->o('ftp_outdir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
                              suffix    => "protein.fa",
                              params    => "-pep",
                            },
      -flow_into         => {
                             	4 => ['Gzip'],
                                },
      -rc_name           => '2Gb_job',
    },
    {
      -logic_name        => 'DumpGFF3',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_gff3'),
                              out_dir    => $self->o('ftp_outdir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
                              suffix    => "annotations.gff3",
                            },
      -flow_into         => { 
                                4 => ['Gzip'],
                            },
      -rc_name           => '4Gb_job',
    },
    {
      -logic_name        => 'DumpGTF',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_gtf'),
                              out_dir    => $self->o('ftp_outdir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
                              suffix    => "canonical_geneset.gtf",
                            },
      -flow_into         => { 
                                4 => ['Gzip'],
                            },
      -rc_name           => '4Gb_job',
    },
    {
      -logic_name        => 'Gzip',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -batch_size        => 10,
      -max_retry_count   => 0,
      -parameters        => {
                              cmd => 'gzip -n -f #out_file#',
                            },
      -rc_name           => '2Gb_job',
    },
  ];
}


sub resource_classes {
    return {
         'default'      => {'LSF' => '-q production-rh7 -C0 -M1000   -R"select[mem>1000]   rusage[mem=1000]"' },
         '2Gb_job'      => {'LSF' => '-q production-rh7 -C0 -M2000   -R"select[mem>2000]   rusage[mem=2000]"' },
         '4Gb_job'      => {'LSF' => '-q production-rh7 -C0 -M4000   -R"select[mem>4000]   rusage[mem=4000]"' },
    };
}

1;
