
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
    max_hive_capacity => 150,

    division => [],
    dump_genome     => $ENV{WORM_CODE}."/scripts/ENSEMBL/scripts/dump_genome.pl",
    dump_transcript => $ENV{WORM_CODE}."/scripts/ENSEMBL/scripts/dump_transcripts.pl",
    dump_gff3       => $ENV{WORM_CODE}."/scripts/ENSEMBL/scripts/dump_gff3.pl",
    dump_gtf        => $ENV{WORM_CODE}."/scripts/ENSEMBL/scripts/dump_gtf_from_ensembl.pl",
    wbps_version    => $ENV{PARASITE_VERSION},
    
    pipeline_db => {
      -host   => $self->o('host'),
      -port   => $self->o('port'),
      -user   => $self->o('user'),
      -pass   => $self->o('pass'),
      -dbname => $ENV{'USER'} . "_parasitedump_" .$self->o('wbps_version'),
      -driver => 'mysql',
    },

    species => [],
    antispecies => [],
    ftp_outdir => '',
    blast_outdir => '',
    dump_ftp => join("," , 
      qw/DumpGenome
      DumpGenomeMasked
      DumpGenomeSoftMasked
      DumpTranscriptsCDS
      DumpTranscriptsmRNA
      DumpProteins
      DumpGFF3
      DumpGTF/),
    dump_blast =>join( ",", 
     qw/DumpGenomeBlast
      DumpGenomeMaskedBlast
      DumpTranscriptsBlast
      DumpProteinsBlast/),
  };
}

# Force an automatic loading of the registry in all workers.
sub beekeeper_extra_cmdline_options {
  my $self = shift;

  return "-reg_conf ".$self->o("registry");
}

sub pipeline_analyses {
  my ($self) = @_;

  my $ini = [];
  
  if (!$self->o('ftp_outdir') && ! $self->o('blast_outdir')) {
    die("Required options: at least one of ftp_outdir, blast_outdir. Nothing do to otherwise");
  }
  if ($self->o('ftp_outdir')) {
    push @$ini, 'SpeciesFactoryFTP';
  }
  if ($self->o('blast_outdir')) {
    push @$ini, 'SpeciesFactoryBlast';
  }

  my @ftp_analysis = $self->o('ftp_outdir') ? split "," , $self->o('dump_ftp') : ();
  my @blast_analysis = $self->o('blast_outdir') ? split "," , $self->o('dump_blast') : ();
  
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
                              species         => $self->o('species'),
                              antispecies     => $self->o('antispecies'),
                              division        => $self->o('division'),
                              meta_filters    => {},
                              chromosome_flow => 0,
                              variation_flow  => 0,
			      regulation_flow => 0,
			      core_flow       => 4, 
                            },
      -flow_into         => {
				4 => \@blast_analysis,
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
      -logic_name        => 'SpeciesFactoryFTP',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -max_retry_count   => 1,
      -parameters        => {
                              species         => $self->o('species'),
                              antispecies     => $self->o('antispecies'),
                              division        => $self->o('division'),
                              meta_filters    => {},
                              chromosome_flow => 0,
                              variation_flow  => 0,
			      regulation_flow => 0,
			      core_flow => 4, 
                            },
      -flow_into         => {
				4 => \@ftp_analysis,
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
                              params => "",
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
                              params => "",
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
