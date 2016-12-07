=cut


=pod

=head1 NAME
Bio::EnsEMBL::EGPipeline::PipeConfig::ParasiteFTP_conf

=head1 DESCRIPTION

Configuration for running the ParaSite FTP dumping pipeline.

=head1 Author

Myriam Shafie

=cut

package Bio::EnsEMBL::EGPipeline::PipeConfig::ParasiteFTP_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::Version 2.3;
use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');
use File::Spec::Functions qw(catdir);

sub default_options {
  my ($self) = @_;
  return {

    %{$self->SUPER::default_options},
    max_hive_capacity => 100,
    pipeline_name => 'ftp_dump_'.$self->o('ensembl_release'),

    #species factory options
    species => [],
    antispecies => [],
    division => [],
    run_all => 0,
    meta_filters => {},

    #vep options
    species_flags => {},
    vep_command  => '--build all',
    hc_random => 0.05,
    pipeline_dir => "/nfs/nobackup/ensemblgenomes/wormbase/parasite/hive_test/vep",
#    ensembl_cvs_root_dir => $ENV{ENSEMBL_CVS_ROOT_DIR},
    ensembl_cvs_root_dir => "/homes/ms41/Hive/",
    perl_command => $^X,
    refseq => 0,
    merged => 0,
    convert => 0,
    debug => 0,
    eg => 1,

    out_dir => "/nfs/nobackup/ensemblgenomes/wormbase/parasite/hive_test/wormbase/staging/parasite/releases",
    blast_dir => "/nfs/nobackup/ensemblgenomes/wormbase/parasite/hive_test/blast",

    sym_script => $ENV{WORM_CODE}."/parasite/scripts/production/make_parasite_FTP_site_symlinks.pl",
    dump_genome => $ENV{WORM_CODE}."/scripts/ENSEMBL/scripts/dump_genome.pl",
    dump_transcript => $ENV{WORM_CODE}."/scripts/ENSEMBL/scripts/dump_transcripts.pl",
    dump_gff3 => $ENV{WORM_CODE}."/scripts/ENSEMBL/scripts/dump_gff3.pl",
    dump_gtf => $ENV{WORM_CODE}."/scripts/ENSEMBL/scripts/dump_gtf_from_ensembl.pl",


    no_genome	   => 0,
    no_genome_masked	  => 0,
    no_genome_softmasked        => 0,
    no_cds => 0,
    no_mrna => 0,
    no_protein     => 0,
    no_gff3        => 0,
    no_gtf         => 0,
    analysis    => [],
    create_core_symlinks => 1,
    no_genome_blast => 0,
    no_genome_masked_blast => 0,
    no_mrna_blast => 0,
    no_protein_blast => 0,
    no_vep => 0,

};
}
# Force an automatic loading of the registry in all workers.
sub beekeeper_extra_cmdline_options {
  my $self = shift;
  return "-reg_conf ".$self->o("registry");
}

sub pipeline_analyses {
  my ($self) = @_;

  my @common_params = map {$_ => $self->o($_) || undef} qw(
    ensembl_release 
    ensembl_cvs_root_dir
    pipeline_dir
    perl_command
    refseq
    merged
    convert
    debug
  );

  my $analysis = [];
  if (!$self->o('no_genome')) {
    push @$analysis, 'DumpGenome';
  }
  if (!$self->o('no_genome_blast')) {
    push @$analysis, 'DumpGenomeBlast';
  }
  if (!$self->o('no_genome_masked')) {
    push @$analysis, 'DumpGenomeMasked';
  }
  if (!$self->o('no_genome_softmasked')) {
    push @$analysis, 'DumpGenomeSoftMasked';
  }
  if (!$self->o('no_cds')) {
    push @$analysis, 'DumpTranscriptsCDS';
  }
  if (!$self->o('no_mrna')) {
    push @$analysis, 'DumpTranscriptsmRNA';
  }
  if (!$self->o('no_protein')) {
    push @$analysis, 'DumpProteins';
  }
  if (!$self->o('no_gff3')) {
    push @$analysis, 'DumpGFF3';
  }
  if (!$self->o('no_gtf')) {
    push @$analysis, 'DumpGTF';
  }
  if (!$self->o('no_genome_masked_blast')) {
    push @$analysis, 'DumpGenomeMaskedBlast';
  }
  if (!$self->o('no_mrna_blast')) {
    push @$analysis, 'DumpTranscriptsBlast';
  }
  if (!$self->o('no_protein_blast')) {
    push @$analysis, 'DumpProteinsBlast';
  }
  if (!$self->o('no_vep')) {
    push @$analysis, 'DumpVEP';
  }

  my $flow;
  if ($self->o('create_core_symlinks')) {
     $flow={'1->A' => ['SpeciesFactory'],'A->1' => ['MakeCoreSymlinks']};	
  }
  else {
     $flow={1 => ['SpeciesFactory']};
  }

  return [
    {
      -logic_name      => 'InitialisePipeline',
      -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -max_retry_count => 0,
      -input_ids       => [{}],
      -flow_into       => $flow,
      -meadow_type     => 'LOCAL',
    },
    {
      -logic_name        => 'SpeciesFactory',
      -module            => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -max_retry_count   => 1,
      -parameters        => {
                              species         => $self->o('species'),
                              antispecies     => $self->o('antispecies'),
                              division        => $self->o('division'),
                              run_all         => $self->o('run_all'),
                              meta_filters    => $self->o('meta_filters'),
                              chromosome_flow => 0,
                              variation_flow  => 0,
                            },
      -flow_into         => {
				2 => $analysis,
			    },
      -meadow_type	 => 'LOCAL',
    },
    {
      -logic_name    => 'DumpVEP',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module        => 'Bio::EnsEMBL::Variation::Pipeline::DumpVEP::DumpVEP',
      -parameters    => {
        species_flags  => $self->o('species_flags'),
        vep_command    => $self->o('vep_command'),
        hc_random      => $self->o('hc_random'),
        @common_params,
	eg	       => $self->o('eg'),
        eg_version     => $ENV{PARASITE_VERSION},
      },
      -rc_name       => 'default',
    },
    {
      -logic_name        => 'DumpGenome',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_genome'),
                              out_dir    => $self->o('out_dir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
			      suffix    => "genomic.fa",
                            },
      -flow_into         => { 
				2 => ['Gzip'],
				},
    },
    {
      -logic_name        => 'DumpGenomeMasked',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_genome'),
                              out_dir    => $self->o('out_dir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
                              suffix    => "genomic_masked.fa",
			      params    => "-mask",
                            },
      -flow_into         => { 
                                2 => ['Gzip'],
                                },
    },
    {
      -logic_name        => 'DumpGenomeSoftMasked',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_genome'),
                              out_dir    => $self->o('out_dir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
                              suffix    => "genomic_softmasked.fa",
                              params    => "-softmask",
                            },
      -flow_into         => { 
                                2 => ['Gzip'],
                                },
    },
    {
      -logic_name        => 'DumpTranscriptsCDS',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_transcript'),
                              out_dir    => $self->o('out_dir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
                              suffix    => "CDS_transcripts.fa",
                              params    => "-cds",
                            },
          -flow_into         => { 
                                2 => ['Gzip'],
                                },
    },
    {
      -logic_name        => 'DumpTranscriptsmRNA',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_transcript'),
                              out_dir    => $self->o('out_dir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
                              suffix    => "mRNA_transcripts.fa",
                              params    => "-mrna",
                            },
      -flow_into         => { 
                                2 => ['Gzip'],
                                },
    },
    {
      -logic_name        => 'DumpProteins',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_transcript'),
                              out_dir    => $self->o('out_dir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
                              suffix    => "protein.fa",
                              params    => "-pep",
                            },
      -flow_into         => {
                             	2 => ['Gzip'],
                                },
    },
    {
      -logic_name        => 'DumpGFF3',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_gff3'),
                              out_dir    => $self->o('out_dir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
                              suffix    => "annotations.gff3",
                            },
      -flow_into         => { 
                                2 => ['Gzip'],
                                },
    },
    {
      -logic_name        => 'DumpGTF',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::ParasiteDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_gtf'),
                              out_dir    => $self->o('out_dir'),
                              ps_rel    => $ENV{PARASITE_VERSION},
                              suffix    => "canonical_geneset.gtf",
                            },
      -flow_into         => { 
                                2 => ['Gzip'],
                                },
    },
    {
      -logic_name        => 'DumpGenomeBlast',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::BlastDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_genome'),
                              blast_dir => $self->o('blast_dir'),
                              suffix    => "dna.toplevel.fa",
                              params    => "-ebiblastheader WBPS",
                            },
    },
    {
      -logic_name        => 'DumpGenomeMaskedBlast',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::BlastDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_genome'),
                              blast_dir => $self->o('blast_dir'),
                              suffix    => "dna_rm.toplevel.fa",
                              params    => "-mask -ebiblastheader WBPS",
                            },
    },
    {
      -logic_name        => 'DumpTranscriptsBlast',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::BlastDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_transcript'),
                              blast_dir => $self->o('blast_dir'),
                              suffix    => "cdna.all.fa",
                              params    => "-mrna -ebiblastheader WBPS",
                            },
    },
    {
      -logic_name        => 'DumpProteinsBlast',
      -hive_capacity     => $self->o('max_hive_capacity'),
      -module            => 'Bio::EnsEMBL::EGPipeline::ParasiteFTP::BlastDumper',
      -max_retry_count   => 1,
      -parameters        => {
                              script     => $self->o('dump_transcript'),
                              blast_dir => $self->o('blast_dir'),
                              suffix    => "pep.all.fa",
                              params    => "-pep -ebiblastheader WBPS",
                            },
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
      -rc_name           => 'normal',
    },
    {
      -logic_name        => 'MakeCoreSymlinks',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -max_retry_count   => 1,
      -parameters        => {
                              out_dir   => $self->o('out_dir'),
                              script    => $self->o('sym_script'),
                              ps_rel    => $ENV{PARASITE_VERSION},
                              wb_rel    => $ENV{WORMBASE_VERSION},
                              cmd	=> "perl #script# -relnum #ps_rel# -wbrelnum #wb_rel# -staging -checksum -stagingdir #out_dir#",
                            }, 
    },

    ];
}
1;
