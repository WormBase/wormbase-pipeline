#!/usr/bin/env perl
#
# test_RNASeq.pl
# 
# by Gary Williams                        
#
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2014-05-01 10:29:23 $      

use strict;
use lib $ENV{'CVS_DIR'};
use Carp;
use Modules::RNASeq;
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species, $new_genome, $check, $runlocally, $notbuild);
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
            "verbose"    => \$verbose,
            "store:s"    => \$store,
            "species:s"  => \$species, # the default is elegans
	    "new_genome" => \$new_genome, # redo all of the alignments against a fresh copy of the genome (this can take 2 weeks for elegans)
	    "check"      => \$check, # test to see if any cufflinks etc. results are missing
	    "runlocally" => \$runlocally, # use the current machine to run the jobs instead of LSF (for getting the last few memory-hungry jobs done)
	    "notbuild"   => \$notbuild, # don't try to make GTF or run cufflinks (for when it is run before doing the Build) 
	   );


$species = 'elegans' unless $species;

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
			   );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


######################################
# variables and command-line options # 
######################################

my $status;
my $data;

#my $database = $wormbase->database('current');
my $database = $wormbase->autoace;

my $RNASeq = RNASeq->new($wormbase, $log, $new_genome, $check);

my $masked = 0;
$log->write_to("Set up genome\n");
$RNASeq->setup_genome_for_star();

if (!$notbuild) {
  $log->write_to("Make GTF file\n");
  $RNASeq->run_make_gtf_transcript($database);
}

$log->write_to("Get experiments from config\n");
my $data = $RNASeq->get_transcribed_long_experiments();

# if -new_genome is set remove all the existing data so it must be aligned again
# otherwise, unless -check is set, remove only the cufflinks data
$log->write_to("Remove old experiment files\n");
$RNASeq->remove_old_experiment_files($data);


# create the LSF Group "/RNASeq/$species" which is now limited to running 15 jobs at a time
$status = system("bgadd -L 15 /RNASeq/$species"); 

my $lsf = LSF::JobManager->new();

my $total = keys %{$data};
my $count_done = 0;

foreach my $experiment_accession (keys %{$data}) {
  
  $count_done++;
  
  # only fire off a job when needed
  if ($check && $RNASeq->check_all_done($experiment_accession, $notbuild)) {
    print "($count_done of $total) Already finished $experiment_accession - not repeating this.\n";
    next;    
  }
  
  # for elegans and briggsae, request 6 Gb memory as 'samtools sort' can take at least 4.5Gb and probably more sometimes
  # for the others, ask for 6 Gb of memory
  my $memory = "6000";
  if ($species ne 'elegans' && $species ne 'briggsae') {
    $memory = "6000";
  }
  
  my $job_name = "worm_".$wormbase->species."_RNASeq";
  my $scratch_dir = $wormbase->logs;
  my $err = "$scratch_dir/RNASeq_align.pl.lsf.${experiment_accession}.err";
  my $out = "$scratch_dir/RNASeq_align.pl.lsf.${experiment_accession}.out";
  my @bsub_options = (-e => "$err", -o => "$out");
  push @bsub_options, (
		       -M => $memory, # in EBI both -M and -R are in Mb
		       -R => "select[mem>$memory] rusage[mem=$memory]",
		       -J => $job_name,
		       -g => "/RNASeq/$species", # add this job to the LSF group '/RNASeq/$species'
		      );
  my $cmd = "RNASeq_align_job.pl -expt $experiment_accession";
  if ($check) {$cmd .= " -check";}
  if ($new_genome) {$cmd .= " -new_genome";}
  if ($notbuild) {$cmd .= " -notbuild";}
  if ($database) {$cmd .= " -database $database";}
  if ($species) {$cmd .= " -species $species";}
  $log->write_to("$cmd\n");
  print "($count_done of $total) Running: $cmd\n";
  if ($runlocally) {
    $cmd .= " -threads 4";
    $wormbase->run_script($cmd, $log);
  } else {
    $cmd = $wormbase->build_cmd($cmd);
    $lsf->submit(@bsub_options, $cmd);
    sleep(10); # wait for a while between submissions so that the ENA FTP server is not overwhelmed by simultaneous requests
  }
}

if (!$runlocally) {
  $lsf->wait_all_children( history => 1 );
  $log->write_to("This set of jobs have completed!\n");
  for my $job ( $lsf->jobs ) {
    if ($job->history->exit_status != 0) {
      $log->write_to("Job $job (" . $job->history->command . ") exited non zero: " . $job->history->exit_status . "\n");
    }
  }
}
$lsf->clear;

####################################################################################
# all analyses have now run - munge results and put them where they are expected
####################################################################################

# now sleep for 1 minute to give the file-system a chance to
# sort out the files we have just written - was having intermittant
# problems with one or two of the intron.ace files as if it hadn't
# finished flushing the file buffers before the file was read in the
# next sections.
sleep 60;

if (!$notbuild) { # in the Build and have GTF and cufflinks done
  make_fpkm();
}

# make the ace file of RNASeq spanned introns to load into acedb
my $splice_file = $wormbase->misc_dynamic."/RNASeq_splice_${species}.ace";
chdir $RNASeq->{RNASeqSRADir};
$status = $wormbase->run_command("rm -f $splice_file", $log);
$status = $wormbase->run_command("cat */Introns/virtual_objects.${species}.RNASeq.ace > $splice_file", $log);
$status = $wormbase->run_script("acezip.pl -file $splice_file", $log);
$status = $wormbase->run_command("cat */Introns/Intron.ace > ${splice_file}.tmp", $log);
$status = $wormbase->run_script("acezip.pl -file ${splice_file}.tmp", $log);
# flatten the results of all libraries at a position into one entry
open (FEAT, "< ${splice_file}.tmp") || $log->log_and_die("Can't open file ${splice_file}.tmp\n");
open (FLAT, ">> $splice_file") || $log->log_and_die("Can't open file $splice_file\n");
my %splice;
while (my $line = <FEAT>) {
  if ($line =~ /^Feature_data/ || $line =~ /^\s*$/) { # new clone
    foreach my $start (keys %splice) {
      foreach my $end (keys %{$splice{$start}}) {
	my $total= 0;
	my $string = "";
	foreach my $library (keys %{$splice{$start}{$end}}) {
	  my $value = $splice{$start}{$end}{$library};
	  $total += $value;
	  $string .= "$library $value "
	}
	# filter out any spurious introns with only 1 or 2 reads
	if ($total > 2) {print FLAT "Feature RNASeq_splice $start $end $total \"$string\"\n";}
      }
    }
    # reset things for the new clone
    %splice = ();
    print FLAT $line;
  } else {
    my @feat = split /\s+/, $line;
    $splice{$feat[2]}{$feat[3]}{$feat[5]} = $feat[4];
  }
}
# and do the last clone
foreach my $start (keys %splice) {
  foreach my $end (keys %{$splice{$start}}) {
    my $total= 0;
    my $string = "";
    foreach my $library (keys %{$splice{$start}{$end}}) {
      my $value = $splice{$start}{$end}{$library};
      $total += $value;
      $string .= "$library $value "
    }
    # filter out any spurious introns with only 1 or 2 reads
    if ($total > 2) {print FLAT "Feature RNASeq_splice $start $end $total \"$string\"\n";}
  }
}
close(FLAT);
close(FEAT);
$status = $wormbase->run_command("rm -f ${splice_file}.tmp", $log);





$log->mail();
print "Finished.\n" if ($verbose);
exit(0);

############################################################################
sub make_fpkm {
  # now write out a table of what has worked
  # and make the expresssion tarball for Wen to put into SPELL
  # and write out .ace files for the FPKM expression levels of genes, transcripts, pseudogenes and CDSs
  
  $log->write_to("\nResults\n");
  $log->write_to("-------\n\n");
  
  chdir $RNASeq->{RNASeqSRADir};
  
  # get the a hashref of the Life_stages etc. of the Condition objects
  my ($life_stage, $condition_species, $sex, $strain, $condition_reference, $tissue) = get_condition_details($database);
  my %life_stages = %{$life_stage};
  my ($sra_study, $sra_experiment, $condition, $analysis_reference) = get_analysis_details($database);
  my %condition_of_analysis = %{$condition};
  my %papers = %{$analysis_reference};
  
  # get the valid CDS and Pseudogenes IDs
  my %CDSs = get_CDSs();
  my %Pseudogenes = get_Pseudogenes();
  
  # open a .ace file to hold the FPKM expression levels of genes, transcripts and CDSs
  my $misc_dynamic = $wormbase->misc_dynamic;
  open (EXPRACE, "> $misc_dynamic/RNASeq_expression_levels_${species}.ace") || $log->log_and_die("Can't open $misc_dynamic/RNASeq_expression_levels_${species}.ace\n");
  
  # write a manifest file for the SPELL data for Wen and Raymond
  open (MANIFEST, "> SPELL_manifest.dat") || $log->log_and_die("Can't open SPELL_manifest.dat\n");
  
  foreach my $experiment_accession (keys %{$data}) {
    
    $log->write_to("$experiment_accession");
    if (-e "$experiment_accession/star_out/accepted_hits.bam.done") {$log->write_to("\tSTAR OK");} else {{$log->write_to("\tSTAR ERROR");}}
    if (-e "$experiment_accession/cufflinks/genes.fpkm_tracking.done") {$log->write_to("\tcufflinks OK");} else {$log->write_to("\tcufflinks ERROR");}
    
    # we have had a problem with incomplete intron.ace files
    if (-e "$experiment_accession/Introns/Intron.ace.done" && ! -z "$experiment_accession/Introns/Intron.ace") {
      my $last_line = `tail -1 "$experiment_accession/Introns/Intron.ace"`;
      my $last_char = chop $last_line;
      if ($last_char ne "\n") {
	$wormbase->run_command("rm -rf $experiment_accession/Introns/Intron.ace", $log);
	$wormbase->run_command("rm -rf $experiment_accession/Introns/Intron.ace.done", $log);
	$log->write_to("\tincomplete intron file removed");
      }
    }
    
    if (-e "$experiment_accession/Introns/Intron.ace.done") {$log->write_to("\tIntrons OK");} else {$log->write_to("\tintron ERROR");}
    $log->write_to("\n");
    
    if (!exists $data->{$experiment_accession}{analysis}) {
      $log->write_to("$experiment_accession lacks an acedb Analysis object - this will not included in the SPELL expression data\n");
      next;
    }
    
    my $analysis = $data->{$experiment_accession}{analysis};
    my $condition_of_analysis = $condition_of_analysis{$analysis};
    my $life_stage = $life_stages{$condition_of_analysis};
    if (!defined $life_stage) {$log->write_to("WARNING: no Condition object found for Analysis object $analysis\n");}
    my $paper = $papers{$analysis};
    if (!defined $paper) {$log->write_to("WARNING: no Reference tag set for Analysis object: $analysis\n");}
    my %Gene_SRX; # used to get the 'best' gene value when there are non-overlapping transcripts causing cufflinks to output two values for one gene
    my %CDS_SRX;
    my %Pseudogene_SRX;
    
    # Wen says: "Gary, the file is good, I just wonder what should we do
    # with the "0" values. SPELL data are all log2 transformed. That is
    # how they are stored in mysql and on the website there are also
    # options for users to disply linear or log2 transformed data.  I will
    # not be able to handle the "0" in SPELL. Can we replace "0" with a
    # very small value, such as 0.0000000001?  If you have a lot tables,
    # you tarzip them together and place them on a ftp site for me to
    # download."
    
    if (-e "$experiment_accession/cufflinks/genes.fpkm_tracking.done") {
      if (!defined $life_stage) {$life_stage=""}
      if (!defined $paper) {$paper=""}
      print MANIFEST "$experiment_accession\t$analysis\t$life_stage\t$paper\n";
      open (EXPR, "<$experiment_accession/cufflinks/genes.fpkm_tracking") || $log->log_and_die("Can't open $experiment_accession/cufflinks/genes.fpkm_tracking\n");
      open (EXPROUT, ">$experiment_accession.out") || $log->log_and_die("Can't open $experiment_accession.out\n");
      while (my $line = <EXPR>) {
	# 0 tracking_id	1 class_code	2 nearest_ref_id	3 gene_id	4 gene_short_name	5 tss_id	6 locus	7 length	8 coverage	9 FPKM	10 FPKM_conf_lo	11 FPKM_conf_hi	12 FPKM_status
	my @f = split /\s+/, $line;
	if ($f[0] eq 'tracking_id') {next;}
	if ($f[0] =~ /CUFF/) {$log->log_and_die("Cufflinks for $experiment_accession has failed to put the Gene IDs in the output file - problem with the GTF file?\n");}
	my $gene_expr = $f[9];
	if ($gene_expr == 0) {$gene_expr = 0.0000000001}
	if ($f[12] ne "FAIL") {
	  if (!exists $Gene_SRX{$f[0]} || $Gene_SRX{$f[0]} < $gene_expr) { # if the gene is repeated, store the highest expression value
	    $Gene_SRX{$f[0]} = $gene_expr;
	  }
	  
	}
      }
      close (EXPR);
      
      # now dump the hash 
      foreach my $gene_id (keys %Gene_SRX) {
	my $gene_expr = $Gene_SRX{$gene_id};
	# print the file for Wen's SPELL data
	print EXPROUT "$gene_id\t$gene_expr\n";
	# and print to the Gene model ace file
	if (defined $life_stage) {
	  print EXPRACE "\nGene : \"$gene_id\"\n";
	  print EXPRACE "RNASeq_FPKM  \"$life_stage\"  \"$gene_expr\"  From_analysis \"$analysis\"\n";
	}
      }
      
      close(EXPROUT);
    }
    
    my $number_of_complaints = 5; #complain only 5 times about Sequences names not in the normal format, otherwise we are swamped by warnings in remanei
    # now get the isoform (coding transcript and CDS) expression values and write out to the ace file
    open (EXPR, "<$experiment_accession/cufflinks/isoforms.fpkm_tracking") || $log->log_and_die("Can't open $experiment_accession/cufflinks/isoforms.fpkm_tracking\n");
    while (my $line = <EXPR>) {
      my @f = split /\s+/, $line;
      if ($f[0] eq 'tracking_id') {next;}
      if ($f[0] =~ /CUFF/) {$log->log_and_die("Cufflinks for $experiment_accession has failed to put the Transcript IDs in the output file - problem with the GTF file?\n");}
      if ($f[12] ne "FAIL") {
	if (defined $life_stage && $life_stage ne "") {
	  my ($sequence_name) = ($f[0] =~ /(^\S+?\.\d+[a-z]?)/);
	  if (!defined $sequence_name) {
	    # complain up to 5 times if it is not even like 'T27B1.t1' - a tRNA gene
	    if ($f[0] !~ /(^\S+?\.t\d+)/ && $number_of_complaints-- > 0 && $species eq 'elegans') {print "Sequence name $f[0] is not in the normal format\n";} 
	    $sequence_name = $f[0]; # use the name as given
	  }
	  
	  # Pseudogenes don't have Transcript objects in acedb, only do CDSs
	  if (exists $CDSs{$sequence_name}) {
	    print EXPRACE "\nTranscript : \"$f[0]\"\n";
	    print EXPRACE "RNASeq_FPKM  \"$life_stage\"  \"$f[9]\"  From_analysis \"$analysis\"\n";
	  }
	  
	  # sum the Pseudogene and CDS scores from their isoforms
	  if (exists $CDSs{$sequence_name}) {
	    $CDS_SRX{$sequence_name} += $f[9]; # sum the FPKMs for this CDS
	  } elsif (exists $Pseudogenes{$sequence_name}) {
	    $Pseudogene_SRX{$sequence_name} += $f[9]; # sum the FPKMs for this Pseudogene
	  }
	}
      }
    }
    close (EXPR);
    
    # now get the CDSs and Pseudogenes that were seen and output their FPKMs to the ace file
    foreach my $CDS (keys %CDS_SRX) {
      if (defined $life_stage && $life_stage ne "") {
	print EXPRACE "\nCDS : \"$CDS\"\n";
	print EXPRACE "RNASeq_FPKM  \"$life_stage\"  \"$CDS_SRX{$CDS}\"  From_analysis \"$analysis\"\n";
      }
    }
    foreach my $Pseudogene (keys %Pseudogene_SRX) {
      if (defined $life_stage && $life_stage ne "") {
	print EXPRACE "\nPseudogene : \"$Pseudogene\"\n";
	print EXPRACE "RNASeq_FPKM  \"$life_stage\"  \"$Pseudogene_SRX{$Pseudogene}\"  From_analysis \"$analysis\"\n";
      }
  }
    
  } # foreach SRX
  $log->write_to("\n");
  
  close(MANIFEST);
  close(EXPRACE);
  
  # make a tarball
  my $out_file = "expr.rnaseq.tar.gz";
  unlink $out_file if -e $out_file;
  $status = $wormbase->run_command("tar zcf $out_file *.out SPELL_manifest.dat", $log);
  my $outdir = $wormbase->spell;
  $status = $wormbase->run_command("cp $out_file $outdir", $log); # this will probably be changed to the autoace/OUTPUT directory soon
}

############################################################################
# my ($sra_study, $sra_experiment, $condition, $analysis_reference) = get_analysis_details($database);
  sub get_analysis_details {

  my ($database) = @_;

  my %analysis_reference;
  my %sra_study;
  my %sra_experiment;
  my %condition;

  my $table_def = &write_paper_def;
  my $table_query = $wormbase->table_maker_query($database, $table_def);
  while(<$table_query>) {
    chomp;
    s/\"//g;  #remove "
    next if (/acedb/ or /\/\//);
    my @data = split(/\t/,$_);
    my ($analysis, $analysis_reference, $sra_or_study, $sraid, $condition) = @data;

    if (!defined $analysis) {next;}

    $analysis_reference{$analysis} = $analysis_reference;
    if (defined $sra_or_study) {
      if ($sra_or_study eq 'Study') {
	$sra_study{$analysis} = $sraid;
      } elsif ($sra_or_study eq 'SRA') {
	$sra_experiment{$analysis} = $sraid;
      }
    }
    $condition{$analysis} = $condition;
  }

  return (\%sra_study, \%sra_experiment, \%condition, \%analysis_reference);
}
############################################################################
#  my ($life_stage, $condition_species, $sex, $strain, $condition_reference) = get_condition_details($database);
sub get_condition_details {

  my ($database) = @_;

  my %life_stage;
  my %condition_species;
  my %sex;
  my %strain;
  my %condition_reference;
  my %tissue;

  my $table_def = &write_life_stage_def;
  my $table_query = $wormbase->table_maker_query($database, $table_def);
  while(<$table_query>) {
    chomp;
    s/\"//g;  #remove "
    next if (/acedb/ or /\/\//);
    my @data = split(/\t/,$_);
    my ($condition, $life_stage, $condition_species, $strain, $condition_reference, $sex, $tissue) = @data;

    if (!defined $condition) {next;}

    $life_stage{$condition} = $life_stage; # multiple life-stages get overwritten leaving just the last one
    $condition_species{$condition} = $condition_species;
    $sex{$condition} = $sex;
    $strain{$condition} = $strain; if (defined $strain{$condition} && $strain{$condition} eq '') {$strain{$condition} = undef}
    $condition_reference{$condition} = $condition_reference;
    $tissue{$condition} = $tissue;
  }

  return (\%life_stage, \%condition_species, \%sex, \%strain, \%condition_reference, \%tissue);
}
############################################################################
# this will write out an acedb tablemaker defn to a temp file
############################################################################
#  my ($life_stage, $condition_species, $strain, $condition_reference) = get_condition_details($database);
sub write_life_stage_def {
  my $def = "/tmp/Life_stages_$$.def";
  open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $species = $wormbase->full_name;
  my $txt = <<END2;


Sortcolumn 1

Colonne 1 
Subtitle Condition    
Width 50 
Optional 
Visible 
Class 
Class Condition 
From 1 
 
Colonne 2 
Subtitle Life_stage   
Width 40 
Mandatory 
Visible 
Class 
Class Life_stage 
From 1 
Tag Life_stage   
 
Colonne 3 
Subtitle Species   
Width 12 
Optional 
Visible 
Class 
Class Species 
From 1 
Tag Species   
 
Colonne 4 
Subtitle Strain   
Width 12 
Optional 
Visible 
Class 
Class Strain 
From 1 
Tag Strain   
 
Colonne 5 
Subtitle Reference   
Width 12 
Optional 
Visible 
Class 
Class Paper 
From 1 
Tag Reference   
 
Colonne 6 
Subtitle Sex  
Width 12 
Optional 
Visible 
Next_Tag 
From 1 
Tag Sex  
 
Colonne 7 
Subtitle Tissue 
Width 12 
Optional 
Visible 
Next_Tag 
From 1 
Tag Tissue 
 
END2

  print TMP $txt;
  close TMP;
  return $def;
}
############################################################################
# this will write out an acedb tablemaker defn to a temp file
############################################################################
# my ($sra_study, $sra_experiment, $condition, $analysis_reference) = get_analysis_details($database);
sub write_paper_def {
  my $def = "/tmp/Life_stages_$$.def";
  open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $species = $wormbase->full_name;
  my $txt = <<END2;

Sortcolumn 1

Colonne 1 
Subtitle Condition 
Width 80 
Optional 
Visible 
Class 
Class Analysis 
From 1 
 
Colonne 2 
Subtitle Reference 
Width 20 
Optional 
Visible 
Class 
Class Paper 
From 1 
Tag Reference  
 
Colonne 3 
Subtitle Database 
Width 12 
Optional 
Hidden 
Show_Tag 
From 1 
Tag Database 
 
Colonne 4 
Subtitle SRA 
Width 12 
Optional 
Hidden 
Class 
Class Database 
Right_of 3 
Tag Database 
 
Colonne 5 
Subtitle SRA_or_Study 
Width 12 
Optional 
Visible 
Class 
Class Database_field 
Right_of 4 
Tag  HERE  
 
Colonne 6 
Subtitle SRX_or_Study_id 
Width 12 
Optional 
Visible 
Class 
Class Text 
Right_of 5 
Tag  HERE  

Colonne 7 
Subtitle Condition 
Width 80 
Optional 
Visible 
Class 
Class Condition 
From 1 
Tag Sample 


END2

  print TMP $txt;
  close TMP;
  return $def;
}
############################################################################

sub get_CDSs {

  my %CDS;
  my $table_def = &write_CDS_def;
  my $table_query = $wormbase->table_maker_query($database, $table_def);
  while(<$table_query>) {
    chomp;
    s/\"//g;  #remove "
    next if (/acedb/ or /\/\//);
    my @data = split("\t",$_);
    my ($CDS) = @data;
    if (!defined $CDS) {next}
    $CDS{$CDS} = -1;
  }

  return %CDS;
}
############################################################################
# this will write out an acedb tablemaker defn to a temp file
############################################################################

sub write_CDS_def {
  my $def = "/tmp/CDS_$$.def";
  open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $species = $wormbase->full_name;
  my $txt = <<END3;

Sortcolumn 1

Colonne 1
//Subtitle CDS
Width 80
Optional
Visible
Class
Class CDS
From 1

END3

  print TMP $txt;
  close TMP;
  return $def;
}
############################################################################

sub get_Pseudogenes {

  my %pseudogenes;
  my $table_def = &write_Pseudogenes_def;
  my $table_query = $wormbase->table_maker_query($database, $table_def);
  while(<$table_query>) {
    chomp;
    s/\"//g;  #remove "
    next if (/acedb/ or /\/\//);
    my @data = split("\t",$_);
    my ($Pseudogene) = @data;
    if (!defined $Pseudogene) {next}
    $pseudogenes{$Pseudogene} = -1;
  }

  return %pseudogenes;
}
############################################################################
# this will write out an acedb tablemaker defn to a temp file
############################################################################

sub write_Pseudogenes_def {
  my $def = "/tmp/Life_stages_$$.def";
  open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $species = $wormbase->full_name;
  my $txt = <<END4;

Sortcolumn 1

Colonne 1
//Subtitle Pseudogene
Width 80
Optional
Visible
Class
Class Pseudogene
From 1

END4

  print TMP $txt;
  close TMP;
  return $def;
}


############################################################################
# do a simple check that the required Analysis and Condition objects are in the geneace database
# useful for checking before a Build that these objects are set up OK
# check Analysis, Condition, Life_stage, Reference (Wen needs these to exist for SPELL)

sub do_analysis_check {

  my $database = $wormbase->database('geneace');
  my ($sra_study, $sra_experiment, $analysis_condition, $analysis_reference) = get_analysis_details($database);
  my ($life_stage, $condition_species, $sex, $strain, $condition_reference, $tissue) = get_condition_details($database);
  
  my %sp = (
	    'elegans'     => 'Caenorhabditis elegans',
	    'briggsae'    => 'Caenorhabditis briggsae',
	    'japonica'    => 'Caenorhabditis japonica',
	    'brenneri'    => 'Caenorhabditis brenneri',
	    'remanei'     => 'Caenorhabditis remanei',
	    'brugia'      => 'Brugia malayi',
	    'ovolvulus'   => 'Onchocerca volvulus',
	   );

  $log->write_to("\n");
  $log->write_to("// Checking Analysis and Condition objects for species: $species\n");
  $log->write_to("\n");

  foreach my $experiment_accession (keys %{$data}) {
    my $analysis = $data->{$experiment_accession}{analysis};
    # get the components of the Analysis name
    my ($RNASeq, $species_name, $strain_name, $life_stage_name, $sex_name, $tissue_name, $sra_study_name, $sra_experiment_name) = split /\./, $analysis;

    if ($experiment_accession ne $sra_experiment_name) {
      $log->write_to("\n// ERROR: mismatch between SRX id '$experiment_accession' in table in script and the SRX id at the end of the analysis name: '$sra_experiment_name'\n//$experiment_accession\t$analysis\n\n");
      next; # if it is that wrong, then don't bother checking other details for the SRX
    }

    # Analysis object details
    if (!defined $sra_study->{$analysis}) {$log->write_to("\n//Missing SRA Study name in Analysis $analysis\n\nAnalysis : $analysis\nDatabase SRA Study $sra_study_name\n"); }
    if (!defined $sra_experiment->{$analysis}) {$log->write_to("\n//Missing SRA Experiment name in Analysis $analysis\n\nAnalysis : $analysis\nDatabase SRA SRA $sra_experiment_name\n"); }
    if (!defined $analysis_reference->{$analysis} || $analysis_reference->{$analysis} eq '') {$log->write_to("\n//Missing Reference in Analysis $analysis\n"); }
    
    if (defined $sra_study->{$analysis} && ($sra_study->{$analysis} ne $sra_study_name)) {$log->write_to("\n//WARNING: $sra_study_name doesn't match SRA Study name '$sra_study->{$analysis}' in Analysis $analysis\n"); }
    if (defined $sra_experiment->{$analysis} && ($sra_experiment->{$analysis} ne $sra_experiment_name)) {$log->write_to("\n//WARNING: $sra_experiment_name doesn't match SRA SRA name '$sra_experiment->{$analysis}' in Analysis $analysis\n"); }

    # Condition object details
    if (!defined $analysis_condition->{$analysis}) {
      $log->write_to("\n//Missing Sample Condition in Analysis $analysis\n\nAnalysis : $analysis\nSample $analysis\n"); 
      $analysis_condition->{$analysis} = $analysis;  # make a new Condition object
    }
    my $condition = $analysis_condition->{$analysis};
    if (!defined $life_stage->{$condition}) {$log->write_to("\n//Missing life_stage in Condition $condition\n\nCondition : $condition\nLife_stage \"$life_stage_name\"\n"); }
    if (defined $life_stage->{$condition} && $life_stage->{$condition} !~ /WBls:\d{7}/) {$log->write_to("\n//ERROR: Invalid life_stage in Condition $condition\n// $life_stage->{$condition}\n\n"); }
    if (defined $life_stage->{$condition} && $life_stage->{$condition} ne $life_stage_name) {$log->write_to("\n//ERROR: mismatch between Life_stage id '$life_stage_name' in table in script and in the Condition $condition\n// $life_stage->{$condition}\n\n"); }
    if (!defined $condition_species->{$condition}) {$log->write_to("\n//Missing species in Condition $condition\n\nCondition : $condition\nSpecies \"$sp{$species_name}\"\n"); }
    if (!defined $sex->{$condition}) {$log->write_to("\n//Missing sex in Condition $condition\n\nCondition : $condition\nSex \"$sex_name\"\n"); }
    if (!defined $tissue->{$condition}) {$log->write_to("\n//Missing tissue in Condition $condition\n\nCondition : $condition\nTissue \"$tissue_name\"\n"); }
    if (!defined $strain->{$condition}) {$log->write_to("\n//Missing strain in Condition $condition\n\nCondition : $condition\nStrain \"$strain_name\"\n"); }
    if (!defined $condition_reference->{$condition} && defined $analysis_reference->{$analysis}) {$log->write_to("\n//Missing Reference in Condition $condition\n\nCondition : $condition\nReference \"$analysis_reference->{$analysis}\"\n"); }

    if (defined $life_stage->{$condition} && ($life_stage->{$condition} ne $life_stage_name)) {$log->write_to("\n//WARNING: $life_stage_name doesn't match Life_stage name '$life_stage->{$condition}' in Condition $condition\n");}
    if (defined $condition_species->{$condition} && ($condition_species->{$condition} ne $sp{$species_name})) {$log->write_to("\n//WARNING: $sp{$species_name} doesn't match Species name '$condition_species->{$condition}' in Condition $condition\n");}
    if (defined $sex->{$condition} && ($sex->{$condition} ne $sex_name)) {$log->write_to("\n//WARNING: $sex_name doesn't match Sex name '$sex->{$condition}' in Condition $condition\n");}
    if (defined $tissue->{$condition} && ($tissue->{$condition} ne $tissue_name)) {$log->write_to("\n//WARNING: $tissue_name doesn't match Tissue name '$tissue->{$condition}' in Condition $condition\n");}
    if (defined $strain->{$condition} && ($strain->{$condition} ne $strain_name)) {$log->write_to("\n//WARNING: $strain_name doesn't match Strain name '$strain->{$condition}' in Condition $condition\n");}
    if (defined $condition_reference->{$condition} && ($condition_reference->{$condition} ne $analysis_reference->{$analysis})) {$log->write_to("\n//WARNING: Analysis $analysis Reference $analysis_reference->{$analysis} doesn't match Reference '$condition_reference->{$condition}' in Condition $condition\n");}
	
 }

  $log->write_to("\n");
  $log->write_to("// Finished Checking Analysis and Condition objects for species: $species\n");
  $log->write_to("\n");

}

############################################################################
