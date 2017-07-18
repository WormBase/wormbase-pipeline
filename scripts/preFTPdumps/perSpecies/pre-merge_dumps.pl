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
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;

use Getopt::Long;

use strict;

my ($species,$store,$debug,$test,$database);
GetOptions(
     'species=s'  => \$species,
     'store=s'    => \$store,
     'debug=s'    => \$debug,
     'test'       => \$test,
     'database=s' => \$database,
)||die(@!);

# WormBase build stanza
my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")
}else{
  $wormbase = Wormbase->new( -debug => $debug, -test => $test,-autoace=> $database,-organism => $species)
}

# global aspects
my $log = Log_files->make_build_log($wormbase);
my $defaultMem = 4000;
my $lsf = LSF::JobManager->new(-M => $defaultMem, 
                               -R => "select[mem>$defaultMem] rusage[mem=$defaultMem]",
                               -J => 'perSpeciesDumps', 
                               -e => '/dev/null',
                               -o => '/dev/null'
                              );

# that is the root of the output file name based on the FTP spec
my $report_dir = $wormbase->reports;

# list of script => (output => (),options => "")
my %script2files = (
  'dump_species_functional_descriptions.pl' => { output => "$report_dir/functional_descriptions.txt" },
  'dump_protein_domains.pl'                 => { output => "$report_dir/protein_domains.csv" },
  'dump_species_orthologs.pl'               => { output => "$report_dir/orthologs.txt" },
  'dump_species_gene_interactions.pl'       => { output => "$report_dir/interactions.txt" },
  'dump_interpolated.pl'                    => { output => "$report_dir/interpolated_clones.txt" },
  'dump_promoters.pl'                       => { output => "$report_dir/potential_promotors.fa" },
  'dump_resource_gene_ids.pl'               => { output => "$report_dir/resource_gene_ids.txt" },
  'dump_swissprot.pl'                       => { output => "$report_dir/swissprot.txt" },
  'dump_ko.pl'                              => { output => "$report_dir/knockout_consortium_alleles.xml" },
  'dump_confirmed_genes.pl'                 => { output => "$report_dir/confirmed_genes.fa" },
  'dump_cdna2orf.pl'                        => { output => "$report_dir/cdna2orf.txt" },
  'dump_geneid_list.pl'                     => { output => "$report_dir/geneIDs.txt" },
  'dump_geneid_list.pl'                     => { output => "$report_dir/geneOtherIDs.txt", -options => "-other" },
  'dump_pcr_list.pl'                        => { output => "$report_dir/pcr_product2gene.txt" },
);

# a.) run
foreach my $script (keys %script2files){
   next unless &check_script($script);
   &clean_previous_output($script2files{$script}->{output});
   &queue_script($script, $script2files{$script});
}

# finish
$log->write_to("Waiting for LSF jobs to finish.\n");
$lsf->wait_all_children( history => 1 );
for my $job ( $lsf->jobs ) {
    if ($job->history->exit_status != 0) {
      $log->write_to("WARNING: Job $job (" . $job->history->command . ") exited non zero: " . $job->history->exit_status . "\n");
    }
}
$lsf->clear;


# check
while (my($script,$options) = each %script2files){
    my $file = $options->{output};
    if (-e $file){
       my $size = -s $file;
       if (-s $file){
         $log->write_to("OK: $script => $file ($size bytes)\n");
       }else{
         $log->write_to("WARNING: $script created an empty $file\n");
       }
    }else{
       $log->write_to("WARNING: $script didn't create $file\n");
    }
}

$log->mail;

# LSF submit $script
sub queue_script{
   my ($script, $data) = @_;
   
   my $cmd = "preFTPdumps/perSpecies/$script";
   $cmd .= " -output $data->{outfile}";
   if (exists $data->{options}) {
     $cmd .= " $data->{options}";
   }
   $cmd = $wormbase->build_cmd($cmd);
   $lsf->submit($cmd);

}

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
     return 1
   }else{
     $log->write_to("ERROR: can't find $path\n");
     return undef
   }
}


