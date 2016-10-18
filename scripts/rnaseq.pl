#!/software/bin/perl -w

# script to test making a feature_data track for RNASeq data
# this script submits one copy of itself for each chromosome/contig run under LSF

# Gary Williams
# 28 Feb 2012


use strict;
use lib $ENV{'CVS_DIR'};
use Getopt::Long;
use Carp;

use Ace;
use Sequence_extract;
use Coords_converter;

use Wormbase;
use Log_files;
use Storable;
use Modules::RNASeq;

use List::Util qw(max);

use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;

# for NameDB_handler
#use lib '/nfs/WWWdev/SANGER_docs/lib/Projects/C_elegans/';
#use NameDB;


my ($help, $test, $verbose, $store, $wormbase, $species, $debug, $outdir);
my ($quick, $method, $chunk_id, $chunk_total, $outfname);

GetOptions ("help"         => \$help,
            "test"         => \$test,
            "verbose"      => \$verbose,
	    "debug:s"      => \$debug,
	    "store:s"      => \$store,
	    "species:s"    => \$species,
	    "outdir:s"     => \$outdir,      # the default is the BUILD/species/acefiles directory

	    "chunkid:s"    => \$chunk_id,    # set when running a sub-job of this script to find the hits in the hits file for a set of chromosomes
	    "chunktotal:s" => \$chunk_total, # set when running a sub-job of this script to find the hits in the hits file for a set of chromosomes
	    "acefname:s"   => \$outfname,    # set when running a sub-job of this script to find the hits in the hits file for a set of chromosomes

            );


#$test = 1;
#$debug = "gw3";


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
                             );
}

if (!defined $species) {$species = $wormbase->species}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

my $acefiles = $wormbase->acefiles;
my $database = $wormbase->autoace;

print "Using database $database\n";

my %F_clones; # forward reads
my %R_clones; # reverse reads
my %F_stranded_clones; # forward stranded reads
my %R_stranded_clones; # reverse stranded reads
  
my $chromosome_prefix = $wormbase->chromosome_prefix;

##########################
# read in the data
##########################

my $RNASeq = RNASeq->new($wormbase, $log, 0, 0);

my $Software = "/nfs/panda/ensemblgenomes/wormbase/software/packages";

my $FEATURE = "RNASeq";

my $RNASeqSRADir    = $RNASeq->{RNASeqSRADir};
my $RNASeqGenomeDir = $RNASeq->{RNASeqGenomeDir};
chdir $RNASeqSRADir;

$outdir = "$acefiles/rnaseq" if (!defined $outdir);

$log->write_to("Initialising data\n");

my $seq_obj = Sequence_extract->invoke($database, 0, $wormbase);

my $coords = Coords_converter->invoke($database, 0, $wormbase);



if ($chunk_id) { # getting the alignments for a set of chromosomes
  ###########################
  # run the sub-job under LSF
  ###########################

  &do_chromosome($chunk_id, $chunk_total, $outfname);
  
} else {
  
  $log->write_to("Parsing reads\n");
  
  
  my $library_count = 0;
  
  #######################
  # create the hits files
  #######################

  my $script = "rnaseq.pl";

  my $lsf;
  my $scratch_dir = $wormbase->logs;
  my $job_name;
  my @bsub_options;
  my $store_file = $wormbase->build_store; # get the store file to use in all commands

  ##########################################################
  # now run the sub-job for each chromosome/contig under LSF
  ##########################################################
  
  system("rm -rf $outdir");
  if (!-e $outdir) {
    mkdir $outdir, 0777;
  }

  # create the LSF Group "/RNASeq/$species" which is now limited to running 15 jobs at a time
#  my $status = system("bgadd -L 15 /RNASeq/$species"); 

  my @chromosomes = $wormbase->get_chromosome_names(-mito => 1, -prefix => 1);

  my $chunk_total = 200; # maximum number of chunked jobs to run

  my $memory = 6000;
  if ($species eq 'elegans') {$memory = 10000}

  $chunk_total = scalar(@chromosomes) if $chunk_total > scalar(@chromosomes);
  $log->write_to("bsub commands . . . . \n\n");
  $lsf = LSF::JobManager->new();
  foreach my $chunk_id (1..$chunk_total) {
    my $batchname = "batch_${chunk_id}";
    my $outfname = "rnaseq_${batchname}.ace";
    my $err = "$scratch_dir/rnaseq.$batchname.err.$$";
    my $out = "$scratch_dir/rnaseq.$batchname.err.$$";

    my $cmd = "$script";
    $cmd .= " -chunkid $chunk_id -chunktotal $chunk_total -acefname $outfname";
    $cmd .= " -outdir $outdir" if ($outdir);
    $cmd .= " -species $species" if ($species);
    $cmd .= " -test" if $test;
    $cmd .= " -debug $debug" if $debug;
    $cmd = $wormbase->build_cmd_line($cmd, $store_file);
    $log->write_to("$cmd\n");
    $log->write_to("logs files:\n $err\n $out\n");
    print "cmd to be executed: $cmd\n";
    @bsub_options = ();
    push @bsub_options, (-M => "$memory", # in EBI both -M and -R are in Gb
			 -R => "select[mem>$memory] rusage[mem=$memory]", 
			 -J => 'rnaseq_alignments',
			 -e => $err,
			 -o => $out,
			 -g => "/RNASeq/$species", # add this job to the LSF group '/RNASeq/$species'
		      );
    $lsf->submit(@bsub_options, $cmd);
  }  


  print "Waiting for LSF jobs to finish.\n";
  $lsf->wait_all_children( history => 1 );
  for my $job ( $lsf->jobs ) {
    if ($job->history->exit_status != 0) {
      $log->write_to("Job $job (" . $job->history->command . ") exited non zero: " . $job->history->exit_status . "\n");
    }
  }
  $lsf->clear;




  ##########################################################
  # concatenate the files to make misc_RNASeq_hits_${species}.ace
  ##########################################################


  $log->write_to("Concatenating the resulting ace files to make misc_RNASeq_hits_${species}.ace\n");
  my $final_file = $wormbase->misc_dynamic."/misc_RNASeq_hits_${species}.ace";
  system("rm -f $final_file");
  system("cat $outdir/RNASeq_*.ace $outdir/virtual_objects_RNASeq_*.ace > $final_file");

}


$log->write_to("Finished.\n");
$log->mail();
exit(0);

##############################################################
# SUBROUTINES
##############################################################

# routine to be executed in sub-job of this script in a LSF queue
sub do_chromosome {

  my ($chunk_id, $chunk_total, $outfname) = @_;

  $log->write_to("Get experiments from config\n");
  my $experiments = $RNASeq->get_transcribed_long_experiments();

  $outdir = "$acefiles/rnaseq" if (!defined $outdir);
  mkdir $outdir, 0777;
  my $output = "$outdir/RNASeq_${outfname}";
  my $Foutput = "$outdir/RNASeq_F_${outfname}";
  my $Routput = "$outdir/RNASeq_R_${outfname}";
  my $Plusoutput = "$outdir/RNASeq_P_${outfname}";
  my $Minusoutput = "$outdir/RNASeq_M_${outfname}";
  my $vfile = "$outdir/virtual_objects_RNASeq_${outfname}";
  open(VIRT, ">$vfile") or $log->log_and_die("Could not open $vfile for writing\n");
  open (ACE, ">$output") || $log->log_and_die("Can't open the file $output\n");
  open (FACE, ">$Foutput") || $log->log_and_die("Can't open the file $Foutput\n");
  open (RACE, ">$Routput") || $log->log_and_die("Can't open the file $Routput\n");
  open (PACE, ">$Plusoutput") || $log->log_and_die("Can't open the file $Plusoutput\n");
  open (MACE, ">$Minusoutput") || $log->log_and_die("Can't open the file $Minusoutput\n");


  my @chromosomes = $wormbase->get_chunked_chroms(-prefix => 0,
                                                  -mito => 1,
						  -chunk_total => $chunk_total,
						  -chunk_id    => $chunk_id);


  foreach my $chromosome (@chromosomes) {
    print "Making the ace files for chromsome $chromosome\n";

    %F_clones = ();
    %R_clones = ();
    %F_stranded_clones = ();
    %R_stranded_clones = ();
    my $library_count = 0;
    my $library_stranded_count = 0;

    # only get the known forward-sense stranded experiments - be very conservative: don't use reverse_stranded, unknown, unstranded, those that havn't been checked yet, etc.
    my $stranded_experiments_1 = $RNASeq->filter_experiments($experiments, 'strandedness', 'stranded');
    my $stranded_experiments = $RNASeq->exclude_experiments($stranded_experiments_1, 'library_layout', 'PAIRED'); # to be on the safe side, exclude all paired-end libraries as their orientations can be exotic

    foreach my $experiment_accession (keys %{$experiments}) {
      if (readhits($experiment_accession, $chromosome)) {
	$library_count++;
      } # count the libraries successfully read
      if (exists $stranded_experiments->{$experiment_accession} ) {
	if (read_stranded_hits($experiment_accession, $chromosome)) {
	  $library_stranded_count++;
	} # count the libraries successfully read
      }
    }
    

    # write the results to an ace file 
    print "writing ace file\n";
    writeace($library_count, $library_stranded_count, $chromosome, $outfname);
  }

  close(ACE);
  close(FACE);
  close(RACE);
  close(PACE);
  close(MACE);
  close(VIRT);

}

##############################################################
# read the aligned reads of the bam file into the clone sequence arrays
# returns '1' if it succeeded in reading the data
sub readhits {
  my ($experiment_accession, $chromosome) = @_;

  my $hitsfile = "$experiment_accession/".$RNASeq->{'alignmentDir'}."/hits.tmp";
  
  if (! -f $hitsfile) {return 0}
  
  # read the hits
  print "\treading hits\n";
  open(BED, "< $hitsfile") || $log->log_and_die("Can't open the file $hitsfile\n");
  
  my $sense;
  my ($clone, $clone_start, $clone_end);
  my ($indel, $change);
  my $clone_last = 0;
  my $in_chromosome = 0;
  my $offset;
  
  my $line_count = 0;
  
  while (my $line = <BED>) {
    if ($line =~ /^track/) {next}
    if ($line =~ /^\s*$/) {next}
    
    my @cols = split /\s+/, $line;
    
    my $chrom = $cols[2];
    my $start = $cols[3] + 1;
    my $end = $cols[4];
    my $reads = $cols[1];
    my $sense = $cols[5];
    
    if ($chromosome ne $chrom) {
      if ($in_chromosome) {
        last;
      } else {
        next;
      } 
    }
      
    $in_chromosome = 1;
    
    if ($sense eq '+') {
      map { $F_clones{$chrom}->[$_] += $reads } ($start..$end); # efficient way of adding the reads value to a range of elements in an array
    } else {
      map { $R_clones{$chrom}->[$_] += $reads } ($start..$end); # efficient way of adding the reads value to a range of elements in an array
    }
    
    $line_count++;
  }    
  
  return 1;
}
#############################################################
# read the aligned stranded reads of the bam file into the clone sequence arrays
# returns '1' if it succeeded in reading the data
sub read_stranded_hits {
  my ($experiment_accession, $chromosome) = @_;

  my $hitsfile = "$experiment_accession/".$RNASeq->{'alignmentDir'}."/hits.stranded";
  
  if (! -f $hitsfile) {return 0}
  
  # read the hits
  print "\treading stranded hits\n";
  open(BED, "< $hitsfile") || $log->log_and_die("Can't open the file $hitsfile\n");
  
  my $sense;
  my ($clone, $clone_start, $clone_end);
  my ($indel, $change);
  my $clone_last = 0;
  my $in_chromosome = 0;
  my $offset;
  
  my $line_count = 0;
  
  while (my $line = <BED>) {
    if ($line =~ /^track/) {next}
    if ($line =~ /^\s*$/) {next}
    
    my @cols = split /\s+/, $line;
    
    my $chrom = $cols[2];
    my $start = $cols[3] + 1;
    my $end = $cols[4];
    my $reads = $cols[1];
    my $sense = $cols[5];
    
    if ($chromosome ne $chrom) {
      if ($in_chromosome) {
        last;
      } else {
        next;
      } 
    }
      
    $in_chromosome = 1;
    
    if ($sense eq '+') {
      map { $F_stranded_clones{$chrom}->[$_] += $reads } ($start..$end); # efficient way of adding the reads value to a range of elements in an array
    } else {
      map { $R_stranded_clones{$chrom}->[$_] += $reads } ($start..$end); # efficient way of adding the reads value to a range of elements in an array
    }
    
    $line_count++;
  }    
  
  return $line_count; # hits.stranded files that are empty count as not existing
}

##############################################################
# write the results to ace files  

sub writeace {
  my ($library_count, $library_stranded_count, $chromosome, $outfname) = @_;

  my (@tiles);

  my $chr_len = $coords->Superlink_length($chromosome_prefix . $chromosome);

  for(my $i=0; $i < $chr_len; $i += 300000) {
    my $chr_start = $i + 1;
    my $chr_end = $chr_start + 300000 - 1;
    $chr_end = $chr_len if $chr_end > $chr_len;
    push @tiles, {
		  start => $chr_start, 
		  end   => $chr_end,
		  segs  => [],
		  F_segs => [],
		  R_segs => [],
		  P_segs => [], # stranded on plus-strand
		  M_segs => [], # stranded on minus-strand
    };
  }

  my @F_counts;
  my @R_counts;
  my @F_stranded_counts;
  my @R_stranded_counts;
  
  @F_counts = (exists $F_clones{$chromosome}) ? @{$F_clones{$chromosome}} : (); # counts in each base in this clone
  @R_counts = (exists $R_clones{$chromosome}) ? @{$R_clones{$chromosome}} : (); # counts in each base in this clone
  @F_stranded_counts = (exists $F_stranded_clones{$chromosome}) ? @{$F_stranded_clones{$chromosome}} : (); # counts in each base in this clone
  @R_stranded_counts = (exists $R_stranded_clones{$chromosome}) ? @{$R_stranded_clones{$chromosome}} : (); # counts in each base in this clone
    
  for( my $tile_idx = 1; $tile_idx <= @tiles; $tile_idx++) {
    my $tile = $tiles[$tile_idx-1];
  
    my $level=0;         # level of output block for average number of reads
    my $clone_start = 0; # start position of growing block - initialise to 0 to show we don't have a block yet

    # and the same variables for the forward and reverse asymmetry output
    my $F_level=0;
    my $F_clone_start = 0; # start position of growing block - initialise to 0 to show we don't have a block yet
    my $R_level=0;
    my $R_clone_start = 0; # start position of growing block - initialise to 0 to show we don't have a block yet

    # and the same variables for the plus and minus strand output
    my $P_level=0;
    my $P_clone_start = 0; # start position of growing block - initialise to 0 to show we don't have a block yet
    my $M_level=0;
    my $M_clone_start = 0; # start position of growing block - initialise to 0 to show we don't have a block yet

    for(my $pos = $tile->{start}; $pos <= $tile->{end}; $pos++) {
      my $counts;
      my $F_asymmetry;
      my $R_asymmetry;
      my $P_counts;
      my $M_counts;

      # first do the sum of the forward and reverse counts (averaged over all libraries)
      if (!defined  $F_counts[$pos]) {$F_counts[$pos] = 0}
      if (!defined  $R_counts[$pos]) {$R_counts[$pos] = 0}
      if (!defined  $F_stranded_counts[$pos]) {$F_stranded_counts[$pos] = 0}
      if (!defined  $R_stranded_counts[$pos]) {$R_stranded_counts[$pos] = 0}

      if ($library_count) {
	$counts = ($F_counts[$pos] + $R_counts[$pos]) / $library_count; # get the average read hits over all libraries
	$counts = int ($counts + 0.5); # this removes regions with < half an average hit
	
	if ($counts > 1.5 * $level || $counts < 0.66 * $level) { # have a substantial change in the amount of reads
	  if ($clone_start != 0 && $level != 0) { # have an existing region that needs to be closed and written
	    push @{$tile->{segs}}, [$clone_start - $tile->{start} + 1, $pos - $tile->{start}, $level];                         
	  }
	  $level = $counts; # set new level
	  
	  if ($level) {
	    $clone_start = $pos;
	  } else {
	    $clone_start = 0; # level is zero, so not growing a region
	  }
	}
      }

      if ($library_stranded_count) {
	$P_counts = $F_stranded_counts[$pos] / $library_stranded_count; # get the average read hits over all libraries
	$P_counts = int ($P_counts + 0.5); # this removes regions with < half an average hit
	
	if ($P_counts > 1.5 * $P_level || $P_counts < 0.66 * $P_level) { # have a substantial change in the amount of reads
	  if ($P_clone_start != 0 && $P_level != 0) { # have an existing region that needs to be closed and written
	    push @{$tile->{P_segs}}, [$P_clone_start - $tile->{start} + 1, $pos - $tile->{start}, $P_level];                         
	  }
	  $P_level = $P_counts; # set new level
	  
	  if ($P_level) {
	    $P_clone_start = $pos;
	  } else {
	    $P_clone_start = 0; # level is zero, so not growing a region
	  }
	}
      }

      if ($library_stranded_count) {
	$M_counts = $R_stranded_counts[$pos] / $library_stranded_count; # get the average read hits over all libraries
	$M_counts = int ($M_counts + 0.5); # this removes regions with < half an average hit
	
	if ($M_counts > 1.5 * $M_level || $M_counts < 0.66 * $M_level) { # have a substantial change in the amount of reads
	  if ($M_clone_start != 0 && $M_level != 0) { # have an existing region that needs to be closed and written
	    push @{$tile->{M_segs}}, [$M_clone_start - $tile->{start} + 1, $pos - $tile->{start}, $M_level];                         
	  }
	  $M_level = $M_counts; # set new level
	  
	  if ($M_level) {
	    $M_clone_start = $pos;
	  } else {
	    $M_clone_start = 0; # level is zero, so not growing a region
	  }
	}
      }

      # now get the asymmetry where forward or reverse is more than $multiple times the other
      my $multiple = 3;
      # increment the counts when zero to stop small values (<$multiple) being displayed when the other sense has a value of zero
      $F_counts[$pos]++ if !$F_counts[$pos];
      $R_counts[$pos]++ if !$R_counts[$pos];
      if ($F_counts[$pos] > $R_counts[$pos] * $multiple) {
	$R_counts[$pos] = 0;
      } elsif ($R_counts[$pos] > $F_counts[$pos] * $multiple) {
	$F_counts[$pos] = 0;
      } else {
	$F_counts[$pos] = 0;
	$R_counts[$pos] = 0;
      }


      if ($F_counts[$pos] > 1.5 * $F_level || $F_counts[$pos] < 0.66 * $F_level) { # have a substantial change in the amount of reads
	if ($F_clone_start != 0 && $F_level != 0) { # have an existing region that needs to be closed and written
          push @{$tile->{F_segs}}, [$F_clone_start - $tile->{start} + 1, $pos - $tile->{start}, $F_level];
	}
	$F_level = $F_counts[$pos]; # set new level
	
	if ($F_level) {
	  $F_clone_start = $pos;
	} else {
	  $F_clone_start = 0; # level is zero, so not growing a region
	}
      }


      if ($R_counts[$pos] > 1.5 * $R_level || $R_counts[$pos] < 0.66 * $R_level) { # have a substantial change in the amount of reads
	if ($R_clone_start != 0 && $R_level != 0) { # have an existing region that needs to be closed and written
          push @{$tile->{R_segs}}, [$R_clone_start - $tile->{start} + 1, $pos - $tile->{start}, $R_level];
	}
	$R_level = $R_counts[$pos]; # set new level
	
	if ($R_level) {
	  $R_clone_start = $pos;
	} else {
	  $R_clone_start = 0; # level is zero, so not growing a region
	}
      }
    }

    if ($clone_start != 0 && $level != 0) { # have an final region that needs to be closed and written
      push @{$tile->{segs}}, [$clone_start - $tile->{start} + 1, $tile->{end} - $tile->{start} + 1, $level];
    }
    if ($P_clone_start != 0 && $P_level != 0) { # have an final region that needs to be closed and written
      push @{$tile->{P_segs}}, [$P_clone_start - $tile->{start} + 1, $tile->{end} - $tile->{start} + 1, $P_level];
    }
    if ($M_clone_start != 0 && $M_level != 0) { # have an final region that needs to be closed and written
      push @{$tile->{M_segs}}, [$M_clone_start - $tile->{start} + 1, $tile->{end} - $tile->{start} + 1, $M_level];
    }
    if ($F_clone_start != 0 && $F_level != 0) { # have an final region that needs to be closed and written
      push @{$tile->{F_segs}}, [$F_clone_start - $tile->{start} + 1, $tile->{end} - $tile->{start} + 1, $F_level];
    }
    if ($R_clone_start != 0 && $R_level != 0) { # have an final region that needs to be closed and written
      push @{$tile->{R_segs}}, [$R_clone_start - $tile->{start} + 1, $tile->{end} - $tile->{start} + 1, $R_level];
    }
  }


  print VIRT "\nSequence : \"${chromosome_prefix}${chromosome}\"\n";
  
  for(my $tile_idx = 1; $tile_idx <= @tiles; $tile_idx++) {
    my $tile = $tiles[$tile_idx-1];

    my $vseq = "${chromosome_prefix}${chromosome}:RNASeq:$tile_idx";
    
    if (@{$tile->{segs}}) {
      printf VIRT "S_Child Feature_data %s %d %d\n", $vseq, $tile->{start}, $tile->{end};
      
      print ACE "\nFeature_data : \"$vseq\"\n";
      foreach my $seg (@{$tile->{segs}}) {
	print ACE "Feature RNASeq @$seg \"Region of RNASeq reads\"\n";
      }
    }
  }


  for(my $tile_idx = 1; $tile_idx <= @tiles; $tile_idx++) {
    my $tile = $tiles[$tile_idx-1];

    my $vseq = "${chromosome_prefix}${chromosome}:RNASeq_stranded_reads:$tile_idx";

    if (@{$tile->{P_segs}}) {
      printf VIRT "S_Child Feature_data %s %d %d\n", $vseq, $tile->{start}, $tile->{end};

      print PACE "\nFeature_data : \"$vseq\"\n";
      foreach my $seg (@{$tile->{P_segs}}) {
	print PACE "Feature RNASeq_stranded @$seg \"Region of forward sense RNASeq reads\"\n";
      }
    }
  }

  for(my $tile_idx = 1; $tile_idx <= @tiles; $tile_idx++) {
    my $tile = $tiles[$tile_idx-1];

    my $vseq = "${chromosome_prefix}${chromosome}:RNASeq_stranded_reads:$tile_idx";

    if (@{$tile->{M_segs}}) {
      printf VIRT "S_Child Feature_data %s %d %d\n", $vseq, $tile->{start}, $tile->{end};
      
      print MACE "\nFeature_data : \"$vseq\"\n";
      foreach my $seg (@{$tile->{M_segs}}) {
	my ($st, $en, $sc) = @$seg;
	print MACE "Feature RNASeq_stranded $en $st $sc \"Region of reverse sense RNASeq reads\"\n"; # reverse the start and end positions for reverse strand
      }
    }
  }

  for(my $tile_idx = 1; $tile_idx <= @tiles; $tile_idx++) {
    my $tile = $tiles[$tile_idx-1];

    my $vseq = "${chromosome_prefix}${chromosome}:RNASeq_forward_reads:$tile_idx";

    if (@{$tile->{F_segs}}) {
      printf VIRT "S_Child Feature_data %s %d %d\n", $vseq, $tile->{start}, $tile->{end};
      
      print FACE "\nFeature_data : \"$vseq\"\n";
      foreach my $seg (@{$tile->{F_segs}}) {
	print FACE "Feature RNASeq_F_asymmetry @$seg \"Region of forward asymmetric RNASeq reads\"\n";
      }
    }
  }

  for(my $tile_idx = 1; $tile_idx <= @tiles; $tile_idx++) {
    my $tile = $tiles[$tile_idx-1];

    my $vseq = "${chromosome_prefix}${chromosome}:RNASeq_reverse_reads:$tile_idx";

    if (@{$tile->{R_segs}}) {
      printf VIRT "S_Child Feature_data %s %d %d\n", $vseq, $tile->{start}, $tile->{end};
      
      print RACE "\nFeature_data : \"$vseq\"\n";
      foreach my $seg (@{$tile->{R_segs}}) {
	print RACE "Feature RNASeq_R_asymmetry @$seg \"Region of reverse asymmetric RNASeq reads\"\n";
      }
    }
  }

}

##########################################
# ($clone, $offset, $clone_last) = get_chrom_details($start, $chrom)



sub get_chrom_details {

  my ($start, $chrom) = @_;

  my $clone = $coords->GetCloneFromCoord($chrom, $start);

  my ($chromosome, $offset) = $coords->CloneOffset($clone);

  $offset--;

  my $clone_last = $coords->Superlink_length($clone);
  $clone_last += $offset;

  return ($clone, $offset, $clone_last);
}


##########################################
sub usage {
  my $error = shift;
  
  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

__END__

=pod

=head2 NAME - rnaseq.pl

=head1 USAGE

=over 4

=item 

=back

This script reads in the RNASeq hits and makes ace files showing the regions with hits (average hits per library) and the regions where there are twice as many forward sense reads as reverese sense reads and vice versa.

script_template.pl MANDATORY arguments:

=over 4

=item none

=back

=head1 REQUIREMENTS

=over 4

=item 

=back

=head1 AUTHOR

=over 4

=item Gary Williams

=back

=cut
