#!/usr/bin/env perl
#
# do the collected per-species dumps for the build
# these go into the REPORTS directory
#
#  it runs scripts that generate one output file per species
#  in the respective SPECIES/REPORTS/ directory
#
#  it will by default overwrite existing output files and warn


use lib $ENV{CVS_DIR};
use lib "$ENV{CVS_DIR}/Modules/Third_party";

use Wormbase;
use Log_files;
use Modules::WormSlurm;

use Getopt::Long;

use strict;

my ($store,$debug,$test,$database,@species);
GetOptions(
     'store=s'    => \$store,
     'debug=s'    => \$debug,
     'test'       => \$test,
     'species=s@' => \@species,
)||die(@!);


# WormBase build stanza
my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n");
}else{
  $wormbase = Wormbase->new( -debug => $debug, 
                             -test => $test);
}

my @script_conf = (
  { script => 'dump_species_functional_descriptions.pl', output => 'functional_descriptions.txt', all => 1,  mem => '12g' },
  { script => 'dump_protein_domains.pl',                 output => 'protein_domains.csv', all => 1 ,  mem => '12g' },
  { script => 'dump_species_orthologs.pl',               output => 'orthologs.txt', all => 1,  mem => '12g'        },
  { script => 'dump_confirmed_genes.pl',                 output => 'confirmed_genes.fa', all => 1 ,  mem => '12g'  },
  { script => 'dump_gpi.pl',                             output => 'gene_product_info.gpi',  all => 1, mem => '12g' },
  { script => 'dump_gpiv2.pl',                           output => 'gene_product_info.gpi2', all => 1, mem => '12g' },
  { script => 'dump_species_gene_interactions.pl',       output => 'interactions.txt',  mem => '12g'               },
  { script => 'dump_interpolated.pl',                    output => 'interpolated_clones.txt',  mem => '12g'        },
  { script => 'dump_swissprot.pl',                       output => 'swissprot.txt',  mem => '12g'                  },
  { script => 'dump_ko.pl',                              output => 'knockout_consortium_alleles.xml',  mem => '12g'},
  { script => 'dump_alaska_ids.pl',                      output => 'alaska_ids.tsv', all => 1,  mem => '12g'       },
  { script => 'dump_cdna2orf.pl',                        output => 'cdna2orf.txt',  mem => '12g'                   },
  { script => 'dump_pcr_list.pl',                        output => 'pcr_product2gene.txt',  mem => '12g'           },
  { script => 'dump_geneid_list.pl',                     output => 'geneIDs.txt', all => 1,  mem => '12g'          },
  { script => 'dump_molecules.pl',                       output => 'molecules.ace',  mem => '12g'                  },
  { script => 'dump_geneid_list.pl',                     output => 'geneOtherIDs.txt', options => '-other', all=>1,  mem => '12g' },
  { script => 'uniprotxrefs.pl',                         output => 'uniprot_papers.txt', all=>1,  mem => '12g'     },
);


# global aspects
my $log = Log_files->make_build_log($wormbase);
my $defaultMem = '4g';
my $slurm_out = $wormbase->build_lsfout;

my %files_to_check;
my %core_species = $wormbase->species_accessors;

my %slurm_jobs;
foreach my $wb ($wormbase, values %core_species ) {
  my $spe = $wb->species;

  next if @species and not grep { $_ eq $spe } @species;

  my $report_dir = $wb->reports;
  
  foreach my $item (@script_conf) {
    my $script = $item->{script};
    my $output = $item->{output};
    # note that all scripts need to be run against autoace, because only that has all
    # of the necessary objects filled in
    my $options = $item->{options} . " -database " . $wormbase->autoace;
    my $mem = (exists $item->{mem}) ? $item->{mem} : $defaultMem;
    next unless $item->{all} or $spe eq 'elegans';
    
    next unless &check_script($script);

    my $outfile = "$report_dir/${spe}.${output}";
    &clean_previous_output($outfile);

    my $cmd = $wb->build_cmd("preFTPdump/perSpecies/$script -outfile $outfile $options");
    my $job_id = WormSlurm::submit_job_with_name($cmd, 'production', $mem, '4:00:00', "${slurm_out}/${script}.${spe}.slurmout", "${slurm_out}/${script}.${spe}.slurmerr", "perSpeciesDumps_${spe}");
    $slurm_jobs{$job_id} = $cmd;
    push @{$files_to_check{$script}}, $outfile;
  }  
}

$log->write_to("Waiting for Slurm jobs to finish.\n");
WormSlurm::wait_for_jobs(keys %slurm_jobs);
for my $job_id (keys %slurm_jobs) {
    my $exit_code = WormSlurm::get_exit_code($job_id);
    if ($exit_code != 0) {
	$log->error("ERROR: Slurm job $job_id (" . $slurm_jobs{$job_id} . ") exited non zero: " . $exit_code . "\n");
    }
}


foreach my $script (keys %files_to_check) {
  foreach my $file (@{$files_to_check{$script}}) {
    if (-e $file){
      my $size = -s $file;
      if (-s $file){
        $log->write_to("OK: $script => $file ($size bytes)\n");
      } else{
        $log->error("ERROR: $script created an empty $file\n");
      }
    } else{
      $log->error("ERROR: $script didn't create $file\n");
    }
  }
}

$log->mail;

# basically try to delete $file
sub clean_previous_output{
  my ($file) = @_;
  if (-e $file){
    unlink $file;
    $log->write_to("cleaning up $file from a previous run\n");
    $log->write_to("can't clean up $file, check permissions\n") if -e $file;
  }
}

# assuming CVS_DIR/scripts/preFTPdumps
sub check_script{
  my ($script) = @_;
  my $path = $ENV{CVS_DIR}."/preFTPdumps/perSpecies/$script";
  if (-e $path){
    return 1;
  } else {
    $log->error("ERROR: can't find $path\n");
    return undef;
  }
}


