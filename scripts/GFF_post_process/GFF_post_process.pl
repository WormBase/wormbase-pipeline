#!/usr/bin/env perl
#
# GFF_post_process.pl
# 
# Coordinates the running of the suite of scripts responisble for decorating 
# and supplementing the "raw" GFF dumped from Ace with additional attributes
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2014-12-15 10:46:36 $
#
# Usage GFF_post_process [-options]


###############################################################################
# variables                                                                     #
###############################################################################

use strict;                                      
use Carp;
use Storable;
use Ace;
use Getopt::Long;

use lib $ENV{CVS_DIR};
use Wormbase;
use Log_files;

##################################################
# Script variables and command-line options      #
##################################################
my ($debug, $test, $verbose, $store, $wormbase, $gff3, $species, $database);
my ($working_dir, $force);

my ($all, 
    $prepare,
    $cleanse_gff,
    $add_landmark,
    $add_motifs,
    $add_gmap2pmap,
    $add_utr,
    $add_supplementary,
    $overload_genomic,
    $overload_transposon,
    $overload_species,
    $overload_mass_spec,
    $overload_cds,
    $overload_gene,
    $overload_snp,
    $overload_tf,
    $overload_rnai,
    $overload_operon,
    $overload_marker,
    $overload_sage,
    $overload_pcr,
    $restructure_mirna,
    $overload_mirna,
    $final,
    );

GetOptions (
  "all"                => \$all,
  "prepare"            => \$prepare,
  "cleanse"            => \$cleanse_gff,
  "overload_blat"      => \$overload_species,
  "overload_massspec"  => \$overload_mass_spec,
  "overload_cds"       => \$overload_cds,
  "overload_gene"       => \$overload_gene,
  "overload_variation" => \$overload_snp,
  "overload_tf"        => \$overload_tf,
  "overload_rnai"      => \$overload_rnai,
  "overload_operon"    => \$overload_operon,
  "overload_marker"    => \$overload_marker,
  "overload_sage"      => \$overload_sage,
  "overload_pcr"       => \$overload_pcr,
  "overload_genomic"   => \$overload_genomic,
  "overload_transposon" => \$overload_transposon,
  "add_landmark"       => \$add_landmark,
  "add_motifs"         => \$add_motifs,
  "add_gmap2pmap"      => \$add_gmap2pmap,
  "add_utr"            => \$add_utr,
  "add_supplementary"  => \$add_supplementary,
  "restructure_mirna"  => \$restructure_mirna,
  "overload_mirna"     => \$overload_mirna,
  "final"              => \$final,

  "workdir=s"   => \$working_dir,
  "force"       => \$force,
  "gff3"        => \$gff3,
  "debug=s"     => \$debug,
  "test"        => \$test,
  "verbose"     => \$verbose,
  "store:s"     => \$store,
  "species:s"   => \$species,
  "database:s"  => \$database,
    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species,
			     );
}

my $log = Log_files->make_build_log($wormbase);

$working_dir = $wormbase->sequences . "/GFF_MUNGE" if not defined $working_dir;
mkdir $working_dir if not -d $working_dir;

my $processed_gff_file = $wormbase->processed_GFF_file( $gff3 );
my $version_string = ($gff3) ? "gff3" : "gff2";
my $pragma_file = "$working_dir/pragmas.$version_string.txt";


############################################################
# Prepare a single GFF file for munging
############################################################
if ($prepare or $all) {
  &prepare();
}


############################################################
# remove unwanted lines from the GFF files
############################################################
if ($cleanse_gff or $all) {
  &run_munging_script("GFF_post_process/cleanse_gff.pl");
}

############################################################
# add Species tags to the BLAT lines
############################################################
if ($overload_species or $all) {
  &run_munging_script("GFF_post_process/overload_gff_blat_species.pl");
}

############################################################
# Adds Interpolated map positions to Genes and Alleles
############################################################
if ($overload_marker or $all) {
  &run_munging_script("GFF_post_process/overload_gff_marker_positions.pl");
}

############################################################
# Overload the Variation lines with consequence etc
############################################################
if ($overload_snp or $all) {
  &run_munging_script("GFF_post_process/overload_gff_variation.pl");
}

############################################################
# Add accessions to the Genomic_canonical  lines
############################################################
if ($cleanse_gff or $all) {
  &run_munging_script("GFF_post_process/overload_gff_genomic.pl");
}

############################################################
# Add gene, transcript and wormpep info to CDS lines
############################################################
if ($overload_cds or $all) {
  &run_munging_script("GFF_post_process/overload_gff_cds_wormpep.pl");
}

############################################################
# Add Locus, Biotype etc to gene lines
############################################################
if ($overload_gene or $all) {
  &run_munging_script("GFF_post_process/overload_gff_gene.pl");
}

############################################################
# Overload the TF-binding sites with TF namer
############################################################
if ($overload_tf or $all) {
  &run_munging_script("GFF_post_process/overload_gff_tf.pl");
}

############################################################
# Add additional information to the RNAi line
############################################################
if ($overload_rnai or $all) {
  &run_munging_script("GFF_post_process/overload_gff_rnai.pl");
}


############################################################
# Add additional information to the Operon lines
############################################################
if ($overload_operon or $all) {
  &run_munging_script("GFF_post_process/overload_gff_operon.pl");
}


############################################################
# Add data source tags to the Mass Spec lines
############################################################
if ($overload_mass_spec or $all) {
  if ($wormbase->species eq 'elegans') {
    &run_munging_script("GFF_post_process/overload_gff_mass_spec.pl");
  }
}

############################################################
# Adds decorations (counts etc) to the SAGE_tag lines
############################################################
if ($overload_sage or $all) {
  if ($wormbase->species eq 'elegans') {
    &run_munging_script("GFF_post_process/overload_gff_sage.pl");
  }
}

############################################################
# Adds decorations (Amplified etc) to the PCR_product lines
############################################################
if ($overload_sage or $all) {
  if ($wormbase->species eq 'elegans') {
    &run_munging_script("GFF_post_process/overload_gff_pcr_product.pl");
  }
}

############################################################
# Adds decorations (Family) to the Transposon span lines
############################################################
if ($overload_transposon or $all) {
    &run_munging_script("GFF_post_process/overload_gff_transposon.pl");
}

############################################################
# Restructures miRNA features into the correct hierarchy
############################################################
if ($restructure_mirna) {
    &run_munging_script("GFF_post_process/restructure_mirnas.pl");
}

############################################################
# Adds miRNA type to mature miRNA transcript lines
############################################################
if ($overload_mirna or $all) {
    &run_munging_script("GFF_post_process/overload_gff_miRNA.pl");
}

#############################################################
# Append pre-computed landmarks, gmap2pmap, protein-motifs and UTRs to the ends of the files
##############################################################
my $splits_dir = $wormbase->gff_splits; 

if ($add_landmark or $all) {
  if ($wormbase->species eq 'elegans') {
    &append_splits_files('landmarks', 'MtDNA');
  }
}

if ($add_gmap2pmap or $all) {
  if ($wormbase->species eq 'elegans') {
    &append_splits_files('gmap2pmap');
  }
}

if ($add_motifs or $all) {
  &append_splits_files('proteinmotifs');
}

if ($add_utr or $all) {
  &append_splits_files('UTR');
}

if ($add_supplementary or $all) {
  &append_supplementary_files();
}

if ($final or $all) {
  &collate_and_sort();
}

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);


###############################
sub check_complete {
  my ($marker) = @_;

  my $m_file = "$working_dir/$marker";
  if (-e $m_file) {
    return 1;
  } else {
    return 0;
  }
}


##############################
sub register_complete {
  my ($marker) = @_;

  my $m_file = "$working_dir/$marker";
  $wormbase->run_command("touch $m_file", $log);
}


#################
sub prepare {
  my (@gff_files, %pragmas, $count);

  return if &check_complete("PROGRESS.${version_string}.prepare") and not $force;

  $log->write_to("GFFmunge - prepare\n");

  if ($wormbase->assembly_type ne 'contig') {
    foreach my $chr ($wormbase->get_chromosome_names(-mito => 1, -prefix => 1)) {
      my $gff = ($gff3) ? $wormbase->GFF3_file_name($chr) : $wormbase->GFF_file_name($chr);
      push @gff_files, $gff;
    }
  } else {
    @gff_files = ($gff3) ? ($wormbase->GFF3_file_name())  : ($wormbase->GFF_file_name());
  }

  $log->write_to(" Collating files: @gff_files\n");

  open(my $out_fh, ">$processed_gff_file") 
      or $log->log_and_die("Could not open $processed_gff_file for writing\n");
  foreach my $gff (@gff_files) {
    open(my $in_fh, $gff) or $log->log_and_die("Could not open $gff for reading\n");
    while(<$in_fh>) {
      /^\#\#/ and do {
        if (not exists $pragmas{$_}) {
          $pragmas{$_} = $count++;
        }
        next;
      };
      print $out_fh $_;
    }
  }
  close($out_fh) or $log->log_and_die("Could not close $processed_gff_file after writing\n");

  open(my $p_fh, ">$pragma_file") or $log->log_and_die("Could not open $pragma_file for writing\n");
  foreach my $pragma (sort { $pragmas{$a} <=> $pragmas{$b} } keys %pragmas) {
    print $p_fh $pragma;
  }
  close($p_fh) or $log->log_and_die("Could not close $pragma_file after writing\n");

  &register_complete("PROGRESS.${version_string}.prepare");
}


#############
sub run_munging_script {
  my ($script) = @_;

  my $descriptor = $script;
  $descriptor =~ s/\.pl//; 
  $descriptor =~ s/^\S+\///;

  return if &check_complete("PROGRESS.${version_string}.$descriptor") and not $force;

  my $outfile = "$working_dir/$descriptor.${version_string}";

  my $cmd = "$script -infile $processed_gff_file -outfile $outfile";
  $cmd .= " -gff3" if $gff3;
  $cmd .= " -database $database" if defined $database;

  my $fail = $wormbase->run_script($cmd, $log);

  $log->log_and_die("Running of $script resulted in failure, so bailing\n")
      if $fail;

  my @errors = $wormbase->check_file(
    $outfile, 
    $log,
    lines => ['^##',
              "^\\S+\\s+\\S+\\s+\\S+\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+"]);

  # VERY paranoid
  if (@errors) {
    $log->log_and_die("GFF did not look clean, so bailing\n");
  }

  if ($debug) {
    $wormbase->run_command("cp -f $outfile $processed_gff_file", $log) and
        $log->log_and_die("Failed to copy $outfile into place\n");
  } else {
    $wormbase->run_command("mv -f $outfile $processed_gff_file", $log) and 
        $log->log_and_die("Failed to move $outfile into place\n");
  }

  &register_complete("PROGRESS.${version_string}.$descriptor");
}


#################################
sub append_splits_files {
  my ($suffix, $skip_pattern) = @_;

  return if &check_complete("PROGRESS.${version_string}.append_${suffix}") and not $force;

  my $outfile = "$working_dir/append_${suffix}.${version_string}";

  my %files;

  foreach my $chr ($wormbase->get_chromosome_names(-mito => 1, -prefix => 1)) {
    next if $skip_pattern and $chr =~ /$skip_pattern/;

    my $gff_file = ($gff3) 
        ? $wormbase->GFF3_file_name($chr, $suffix) 
        : $wormbase->GFF_file_name($chr, $suffix);


    $log->log_and_die("Could not find $gff_file\n")
        if not -e $gff_file;
    $files{$gff_file} = 1;
  }

  my @files = sort keys %files;

  if ($wormbase->run_command("cat $processed_gff_file @files > $outfile", $log)) {
    $log->log_and_die("Failed to append to $processed_gff_file: @files\n");
  }


  if ($debug) { 
    $wormbase->run_command("cp -f $outfile $processed_gff_file", $log) 
        and $log->log_and_die("Failed to copy the GFF_SPLITS-appended GFF file into place\n");
  } else {
    $wormbase->run_command("mv -f $outfile $processed_gff_file", $log) 
        and $log->log_and_die("Failed to move the GFF_SPLITS-appended GFF file into place\n");
  }

  &register_complete("PROGRESS.${version_string}.append_${suffix}");
}

##################################
sub append_supplementary_files {

  return if &check_complete("PROGRESS.${version_string}.append_supplementary") and not $force;

  my $outfile = "$working_dir/append_supplementary.${version_string}";

  my (@gfffiles);
  
  my $supdir1 = $wormbase->misc_dynamic . "/SUPPLEMENTARY_GFF";
  my $spec = $wormbase->species;

  if ($gff3) {
    @gfffiles = glob("$supdir1/${spec}.*.gff3");
  } else {
    @gfffiles = glob("$supdir1/${spec}.*.gff2");
  }

  if ($wormbase->run_command("cat $processed_gff_file @gfffiles > $outfile", $log)) {
    $log->log_and_die("Failed to appened supplementary GFF files to $processed_gff_file : @gfffiles\n");
  }


  if ($debug) { 
    $wormbase->run_command("cp -f $outfile $processed_gff_file", $log) 
        and $log->log_and_die("Failed to copy the SUPPLEMENTARY-appended GFF file into place\n");
  } else {
    $wormbase->run_command("mv -f $outfile $processed_gff_file", $log) 
        and $log->log_and_die("Failed to move the SUPPLEMENTARY-appended GFF file into place\n");
  }

  &register_complete("PROGRESS.${version_string}.append_supplementary");
}

############################
sub collate_and_sort {
  return if &check_complete("PROGRESS.${version_string}.collate_and_sort") and not $force;

  my $outfile = "$working_dir/collate_and_sort.${version_string}";
  open(my $out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");

  open(my $gff_in, $pragma_file) or $log->log_and_die("Could not open $pragma_file for reading\n");
  while(<$gff_in>) {
    print $out_fh $_;
  }
  
  if ($gff3) {
    open($gff_in, "sort -S 3G -k 1,1 -k4,4n -k 5,5n $processed_gff_file | perl $ENV{CVS_DIR}/GFF_post_process/sort_gff3.pl - |")
        or $log->log_and_die("Could not open sort_gff3.pl command for $processed_gff_file\n");
  } else {
    open($gff_in, "sort -S 3G -k 1,1 -k4,4n -k 5,5n $processed_gff_file |") 
        or $log->log_and_die("Could not open sort cmd for $processed_gff_file\n");
  }
  while(<$gff_in>) {
    print $out_fh $_;
  }
  close($out_fh) or $log->log_and_die("Could not close $outfile after appending\n");

  if ($debug) { 
    $wormbase->run_command("cp -f $outfile $processed_gff_file", $log) 
        and $log->log_and_die("Failed to copy the sorted/collated GFF file into place\n");
  } else {
    $wormbase->run_command("mv -f $outfile $processed_gff_file", $log) 
        and $log->log_and_die("Failed to move the sorted/collated GFF file into place\n");
  }


  #
  # Run the GFF3 post-processing script; do it here rather than earlier because this is a
  # compulsory last step
  #
  if ($gff3) {
    my $outfile = "$working_dir/post_process_gff3.gff3";
    $wormbase->run_script("GFF_post_process/post_process_gff3.pl -infile $processed_gff_file -outfile $outfile", $log)
        and $log->log_and_die("Unsuccessful post-processing of GFF3 prior to collation\n");

    if ($debug) {
      $wormbase->run_command("cp -f $outfile $processed_gff_file", $log)
          and $log->log_and_die("Failed to copy the GFF3-post-processed file into place");
    } else {
      $wormbase->run_command("mv -f $outfile $processed_gff_file", $log)
          and $log->log_and_die("Failed to move the GFF3-post-processed file into place");
    }
  }

  $wormbase->run_command("gzip -n -9 $processed_gff_file", $log)
      and $log->log_and_die("Failed to gzip the final processed GFF file\n");

  &register_complete("PROGRESS.${version_string}.collate_and_sort");
}

