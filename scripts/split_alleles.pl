#!/usr/bin/env perl 
#===============================================================================
#
#         FILE:  split_alleles.pl
#
#        USAGE:  ./split_alleles.pl 
#
#  DESCRIPTION: creates bins of alleles for mapping and submits them to LSF
#
#       AUTHOR:   (Michael Paulini), <mh6@sanger.ac.uk>
#      COMPANY:  
#      CREATED:  13/05/10 12:14:18 BST
#===============================================================================

use Ace;
use IO::File;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin";
use Wormbase;
use Modules::map_Alleles;

use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;

use strict;

sub print_usage{
print  <<USAGE;
split_alleles.pl options:
	-debug USER_NAME             sets email address and debug mode
	-store FILE_NAME             use a Storable wormbase configuration file
	-outdir DIR_NAME             print allele_mapping_VERSION.ace to DIR_NAME
	-database DATABASE_DIRECTORY use a different AceDB
	-help                        print this message
	-test                        use the test database
	-species SPECIES_NAME        specify a non-elegans species
USAGE

exit 1;	
}

my ( $debug, $store,$database,$help,$test,$species,$wb,$noload,$outdir,$binsize,
     $read_ids, $write_ids, $map, $load_collated);

GetOptions(
  'species=s'=> \$species,
  'debug=s'  => \$debug,
  'store=s'  => \$store,
  'outdir=s' => \$outdir,
  'database=s'  => \$database,
  'help'        => \$help,
  'test'        => \$test,
  'binsize=s'   => \$binsize,
  'readids=s'     => \$read_ids,
  'writeids=s'    => \$write_ids,
  'map'           => \$map,
  'loadcollated' => \$load_collated,
) or &print_usage();

&print_usage if $help;

if ($store) {
  $wb = Storable::retrieve($store) 
      or croak("cannot restore wormbase from $store");
}
else { 
  $wb = Wormbase->new( -debug => $debug, 
                       -test => $test, 
                       -organism => $species, 
                       -autoace => $database ); 
}

my $log = Log_files->make_build_log($wb);

if ($debug) {
    print "DEBUG \"$debug\"\n\n";
}

$outdir = $wb->autoace."/TMP" if not defined $outdir;
mkdir $outdir if not -d $outdir;

$database = $wb->autoace() if not defined $database;
$species = $wb->species if not defined $species;

my $collated_out_file = $wb->acefiles . "/map_alleles_output.ace";

my $variations;
if ($read_ids or $write_ids) {
  if ($read_ids and $write_ids) {
    $log->log_and_die("You cannot supply both -readids and -writeids. Choose one\n");
  }

  if ($read_ids) {
    #
    # read_ids
    #
    $variations = [];

    open(my $idfh, $read_ids) or $log->log_and_die("Could not open $read_ids for reading\n");
    while(<$idfh>) {
      /^\#/ and next;
      /^(\S+)/ and push @$variations, $1;
    }
  } else {
    #
    # write ids
    #
    open(my $idfh, ">$write_ids") or $log->log_and_die("Could not open $write_ids for writing\n");
    $variations = &get_all_allele_ids_table_maker();
    foreach my $var_id (@$variations) {
      print $idfh "$var_id\n";
    }
    close($idfh);
    $log->mail();
    exit(0);
  }
} else {
  $variations = &get_all_allele_ids_table_maker();
}

if ($map) {
  $binsize = 50000 if not defined $binsize;
  
  my (@bins, @id_files, @out_files);
  
  while (my $a = shift @$variations){
    if (not @bins or scalar(@{$bins[-1]}) == $binsize) {
      push @bins, [];
    }
    
    push @{$bins[-1]}, $a;
  }
  
  my (@jobs, %lsf_jobs);
  
  
  for(my $bn=1; $bn <= @bins; $bn++) {
    my $list = $bins[$bn-1];
    
    my $id_file = "$outdir/map_alleles.$$.${species}.${bn}.txt";
    my $out_file = "$outdir/allele_mapping.$$.${species}.${bn}.ace";
    
    open(my $fh, ">$id_file") or 
        $log->log_and_die("Could not open IDfile $id_file for writing\n");  
    foreach my $nm (@$list) {
      print $fh "$nm\n";
    }
    close($fh);
    
    push @jobs, {
      bin => $bn,
      idfile => $id_file,
      outfile => $out_file,
      success => 0,
    };
    
  }
  
  foreach my $mem_gb (5, 10) {
    my @outstanding_jobs = grep { not $_->{success} } @jobs;
    
    if (@outstanding_jobs) {
      my $lsf = LSF::JobManager->new();
      
      foreach my $job (@outstanding_jobs) {
        my $idfile = $job->{idfile};
        my $outfile = $job->{outfile};
        my $batch = $job->{bin};
        
        my $cmd =  "perl $Bin/map_Alleles.pl -idfile $idfile -outfile $outfile -noload -species $species";
        $cmd .= ' -noremap' unless $wb->species eq 'elegans'; # bit crude, but it works
        $cmd .= " -debug $debug" if $debug;
        $cmd .= " -test"  if $test;
        
        my $job_name = "${species}_splitalleles_${batch}_${mem_gb}GB";
        my $lsf_output = "${outdir}/${job_name}.lsf_out";
        my @bsub_options =(-o => $lsf_output,
                           -M => $mem_gb * 1000, 
                           -R => sprintf("select[mem>%d] rusage[mem=%d]", $mem_gb * 1000, $mem_gb * 1000),
                           -J => $job_name);
        $lsf->submit(@bsub_options, $cmd);
        $job->{lsfout} = $lsf_output;
      }
      
      $lsf->wait_all_children();
      # Pause here to make sure that LSF has fully registered completion of the jobs
      sleep(10);
      $lsf->clear;
      
      foreach my $job (@outstanding_jobs) {
        my $status = &check_exit_status($job->{lsfout});
        if ($status == 0) {
          $job->{success} = 1;
        }
      }
    }
  }
  
  my (@success_out, @success_id, @success_lsfout);
  
  my $errors = 0;
  foreach my $job (@jobs) {
    if ($job->{success}) {
      push @success_out, $job->{outfile};
      push @success_id, $job->{idfile};
      push @success_lsfout, $job->{lsfout};
    } else {
      $log->error("The following command failed, so retaining files for investigation: " . $job->{cmd} . "\n");
      $errors++;
    }
  }
  unlink @success_lsfout;

  if (not $errors) {
    # Finally, generated collated file
    $wb->run_command("cat @success_out > $collated_out_file", $log);
    unlink @success_out;
  }
}  

if ($load_collated) {
  if (not -e $collated_out_file) {
    $log->log_and_die("Could not locate $collated_out_file for loading\n");
  }
  $wb->load_to_database($wb->autoace, $collated_out_file, 'map_alleles.pl',$log);
}

$log->mail();
exit(0);

#########################
sub get_all_allele_ids_table_maker {
  
  my $species = $wb->full_name;
  my $query = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Mandatory 
Visible 
Class 
Class Variation 
From 1 
Condition Flanking_sequences AND Live AND Species = "$species"
 
EOF

  my $def_file = "/tmp/all_allele_ids.def";
  open(my $tmdeffh, ">$def_file") or $log->log_and_die("Could not open $def_file for writing\n");
  print $tmdeffh $query;
  close($tmdeffh);

  my @list;

  my $tmfh = $wb->table_maker_query($wb->autoace, $def_file);
  while(<$tmfh>) {
    /^\"(\S+)\"$/ and do {
      push @list, $1;
    }
  }
  close($tmfh);

  return \@list;
}

#########################
sub check_exit_status {
  my ($file) = @_;

  my $status = 999;

  open(my $fh, $file) or $log->log_and_die("Could not open $file for reading\n");
  while(<$fh>) {
    if (/Exited with exit code (\d+)/) {
      $status = $1;
    } elsif (/^Exited/) {
      $status = 999;
    } elsif (/Successfully completed/) {
      $status = 0;
    }
  }

  return $status;
}




