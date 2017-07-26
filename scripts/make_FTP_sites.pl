#!/usr/bin/env perl
#
# make_FTP_sites.pl
#
# A PERL wrapper to automate the process of building the FTP sites 
# builds wormbase & wormpep FTP sites
# 
# Last updated by: $Author: klh $
# Last updated on: $Date: 2015-06-02 08:48:04 $
#
# see pod documentation (i.e. 'perldoc make_FTP_sites.pl') for more information.
#
##########################################################################################

=pod

=head1 NAME - make_FTP_sites.pl

=back


=head1 USAGE

=over 4

=item make_FTP_sites.pl

=back


This script does :

 [01] - make a new directory for the WS release
 [02] - copy the WS release files to the target directory
 [03] - make a new directory for the chromosome DNA/GFF/AGP files
 [04] - copy the chromosome DNA/GFF/AGP files to the target directory
 [05] - copy the models.wrm file across (also present in the database.*.tar.gz files)
 [06] - copy the relevant dbcomp file across
 [07] - copy across latest wormpep release
 [08] - make wormpep FTP site
 [09] - copy WormRNA release across
 [10] - extract confirmed genes from autoace and make a file on FTP site
 [11] - delete the old symbolic link and make the new one
 [12] - delete the old WS release version directory
 [13] - makes a file of cDNA2orf connections
 [14] - makes a file of all gene IDs with CGC names and Sequence names (where present)
 [15] - exit gracefully


=over 4

=item MANDATORY arguments:

none

=back

=over 4

=item OPTIONAL arguments:

-help (this help page)


=cut

use strict;
use lib $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'}."/Modules";
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
use IO::Handle;
use File::Path;
use Bio::SeqIO;
use JSON;

#################################################################################
# Command-line options and variables                                            #
#################################################################################

my ($help, $debug, $test, $testout, $store, $wormbase);

my $acedb;   # only copy across acedb files
my $ont;     # only copy across ontology files
my $annots;  # only copy annotations folder files
my $wormpep; # only copy wormpep files
my $homols;  # only copy best blast hits 
my $manifest;# check everything has been copied.
my $all;     # copy everything across
my $dna;
my $rna;
my $xrefs;
my $gff;
my $reports;
my $ests;
my $blastx;
my $letter;
my $md5;
my $multi_species;
my $assembly_manifest;
my (%skip_species, @skip_species, @only_species, %only_species, $WS_version, $WS_version_name);

GetOptions ("help"          => \$help,
	    "debug=s"       => \$debug,
	    "test"          => \$test,   # Copies data from the test env to a test ftp under ~/tmp/pub
            "testout=s"     => \$testout, # Copies real data to test location
	    "store:s"       => \$store,
	    "acedb"         => \$acedb,
            "blastx"        => \$blastx,
	    "multi"         => \$multi_species,
	    "dna"           => \$dna,
	    "rna"           => \$rna,
	    "wormpep"       => \$wormpep,
	    "gff"           => \$gff,
            "reports"       => \$reports,
	    "annots"        => \$annots,
            "ests"          => \$ests,
	    "homols"        => \$homols,
	    "ont"           => \$ont,
            "assmanifest"   => \$assembly_manifest,
            "letter"        => \$letter,
	    "manifest"      => \$manifest,
            "md5"            => \$md5,
            "skipspecies=s@" => \@skip_species,
            "onlyspecies=s@" => \@only_species,
	    "all"            => \$all,
            "wbversion=s"    => \$WS_version,
    )||die(&usage);


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new(-debug   => $debug,
                            -test    => $test,
      );
}

# Display help if required
&usage if ($help);

# establish log file.
my $log = Log_files->make_build_log($wormbase);
map { $skip_species{$_} = 1 } @skip_species;
map { $only_species{$_} = 1 } @only_species;

if (not $WS_version) {
  $WS_version = $wormbase->get_wormbase_version();
}
$WS_version_name = "WS${WS_version}";

my $maintainers = "All";

my $targetdir = ($testout) 
    ? "$testout/releases/$WS_version_name"
    : $wormbase->ftp_staging . "/releases/$WS_version_name";

$log->write_to("WRITING TO $targetdir\n");

#################################################################################
# Main                                                                          #
#################################################################################

&copy_acedb_files        if $acedb or $all;
&copy_blastx             if $blastx or $all;
&copy_multi_species      if $multi_species or $all;
&copy_dna_files          if $dna or $all;
&copy_rna_files          if $rna or $all;
&copy_wormpep_files      if $wormpep or $all;
&copy_gff_files          if $gff or $all;
&copy_report_files       if $reports or $all;
&copy_est_files          if $ests or $all;
&copy_annotations_files  if $annots or $all;
&copy_ontology_files     if $ont or $all;
&copy_homol_data         if $homols or $all;
&copy_assembly_manifest  if $assembly_manifest or $all;
&copy_release_letter     if $letter or $all;
&check_manifest          if $manifest or $all;
&make_md5sums            if $md5 or $all;

$log->mail;
exit(0);

##########################################################
# copy the WS acedb files across and check on the size
# The FTP disk tends to be unstable
##########################################################
sub copy_acedb_files{
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying release files\n");

  my $base_dir = $wormbase->basedir;
  my $ace_dir = $wormbase->autoace;
  my $ftp_acedb_dir = "$targetdir/acedb";

  mkpath("$ftp_acedb_dir",1,0775);

  my $filename;

  opendir (RELEASE,"$ace_dir/release") or $log->log_and_die("Could not open directory $ace_dir/release");
  while (defined($filename = readdir(RELEASE))) {
    if (($filename eq ".")||($filename eq "..")) { next;}
    if (($filename =~ /letter/)||($filename =~ /dbcomp/)) { next;}
    $wormbase->run_command("cp $ace_dir/release/$filename $ftp_acedb_dir/$filename", $log);

    my $O_SIZE = (stat("$ace_dir/release/$filename"))[7];
    my $N_SIZE = (stat("$ftp_acedb_dir/$filename"))[7];
    if ($O_SIZE != $N_SIZE) {
      $log->write_to("\tError: $filename SRC: $O_SIZE TGT: $N_SIZE - different file sizes, please check\n");
      $log->log_and_die("Couldn't copy $filename\n");
    } 
  }
  closedir RELEASE;
  
  # Copy across the models.wrm file
  $wormbase->run_command("cp $ace_dir/wspec/models.wrm $ftp_acedb_dir/models.wrm.$WS_version_name", $log);
  $wormbase->run_command("cp $ace_dir/wspec/models.wrm.annot $ftp_acedb_dir/models.wrm.$WS_version_name.annot", $log);

  # copy some miscellaneous files across
  my $old_release = $WS_version -1;
  $wormbase->run_command("cp ".	$wormbase->compare."/WS$old_release-$WS_version_name.dbcomp $ftp_acedb_dir", $log);
  $wormbase->run_command("cp $base_dir/autoace_config/INSTALL $ftp_acedb_dir", $log);  
  
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying release files\n\n");
  
}

##################################################
# copy the Non-C_elegans blastx data
##################################################
sub copy_blastx {
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying/zipping Non-elegans blastx\n");

  my $ftp_acedb_dir = "$targetdir/acedb";

  my $blastx_dir = "$ftp_acedb_dir/Non_C_elegans_BLASTX";
  mkpath("$blastx_dir",1,0775);

  my %accessors = ($wormbase->species_accessors);
  foreach my $wb (values %accessors) {
    next if exists $skip_species{$wb->species};
    next if @only_species and not exists($only_species{$wb->species});

    my $gspecies = $wb->full_name('-g_species'=>1);
    my $in_file = $wb->acefiles . "/" . $wb->species . "_blastx.ace";

    if (-e $in_file) {
      my $out_file = "$blastx_dir/$gspecies.$WS_version_name.blastx.ace.gz";

      $wormbase->run_command("cat $in_file | gzip -n -9 -c > $out_file", $log);
    }
  }

  $log->write_to("$runtime: Finished copying/zipping non-elegans blastx\n\n");
}


##################################################
# copy the DNA and agp files across
##################################################
sub copy_dna_files{
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying dna and agp files\n");

  my %accessors = ($wormbase->all_species_accessors);
  $accessors{elegans} = $wormbase;
  
  ACC: foreach my $wb (values %accessors) {
    next if exists $skip_species{$wb->species};
    next if @only_species and not exists $only_species{$wb->species};
    
    my $gspecies = $wb->full_name('-g_species'=>1);
    my $bioproj = $wb->ncbi_bioproject;
    my $chromdir = $wb->chromosomes;
    my $seqdir = $wb->sequences;
    
    my %copied_files;
    
    if (-e "$chromdir" and -e "$seqdir") {
      my $dna_dir = "$targetdir/species/$gspecies/$bioproj";
      mkpath($dna_dir,1,0775);
      #todd wants all species to have whole genome in one file

      my $species = $wb->species;
	
      my $dna_file = $wb->genome_seq;
      my $masked_file = $wb->masked_genome_seq;
      my $soft_file = $wb->softmasked_genome_seq;
      
      foreach my $f ($dna_file, $masked_file, $soft_file) {
        if (not -e $f or not -s $f) {
          $log->error("ERROR: Could not find DNA file for $gspecies ($f)\n");
          next ACC;
        }
      }
	
      my $target_dna_file =  "$dna_dir/${gspecies}.${bioproj}.${WS_version_name}.genomic.fa";
      my $target_masked = "$dna_dir/${gspecies}.${bioproj}.${WS_version_name}.genomic_masked.fa";
      my $target_soft = "$dna_dir/${gspecies}.${bioproj}.${WS_version_name}.genomic_softmasked.fa";
      
      foreach my $pair ([$dna_file, $target_dna_file],
                        [$masked_file, $target_masked],
                        [$soft_file, $target_soft]) {
        my ($src, $tgt) = @$pair;
        
        eval {
          my ($read_fh, $write_fh);
          
          open($write_fh, ">$tgt") or die "Could not open $tgt for writing\n";
          
          if ($src =~ /\.gz$/) {
            open($read_fh, "gunzip -c $src |") or die "Could not open gunzip stream to $src\n";
          } else {
            open($read_fh, $src) or die "Could not open $src for reading\n";
          }
          
          while(<$read_fh>) {
            if ($species eq 'elegans') {
              s/^\>CHROMOSOME_(\S+)/>$1/; 
            } elsif ($species eq 'briggsae') {
              s/^\>chr(\S+)/>$1/;
            }
            print $write_fh $_;
          }
          
          close($write_fh) or die "Could not close $tgt after writing\n";

          unlink "${tgt}.gz" if -e "${tgt}.gz";

          $wormbase->run_command("gzip -n -9 $tgt", $log)
              and die "Could not gzip $tgt after copying\n";
        };          
        $@ and do {
          $log->error("ERROR: Could not copy DNA file $src for $species; skipping\n");
          next ACC;
        };
      }
      
      map { $copied_files{$_} = 1 } ($dna_file, $masked_file, $soft_file);
      
      # copy over outstanding dna files
      foreach my $dna_file (glob("$seqdir/*.dna.gz"), glob("$seqdir/*.dna")) {
        if (not exists $copied_files{$dna_file}) {
          my ($prefix) = $dna_file =~ /$seqdir\/(\S+)\.dna/;
          my $target = "$dna_dir/${gspecies}.${bioproj}.${WS_version_name}.$prefix.fa.gz";
          if ($dna_file =~ /\.gz$/) {
            $wormbase->run_command("cp -f $dna_file $target", $log);
          } else {
            $wormbase->run_command("cat $dna_file | gzip -n > $target", $log);
          }
        }
      }
      
      my @agp_files = glob("$chromdir/*.agp");
      
      if (@agp_files) {
        my $target_agp_file =  "$dna_dir/${gspecies}.${bioproj}.${WS_version_name}.assembly.agp"; 
        
        if (scalar(@agp_files) == 1) {
          # just the one - assume its for whole genome and copy it across
          my ($single_file) = @agp_files;
          $wormbase->run_command("cat $single_file | gzip -n -9 -c > ${target_agp_file}.gz", $log); 
        } else {
          # assume per-chromosome
          unlink $target_agp_file;
          $wormbase->run_command("touch $target_agp_file", $log);
          foreach my $chrom ($wb->get_chromosome_names(-mito => 0, -prefix => 1)) {
            my $agp = "$chromdir/$chrom.agp";
            if (-e "$chromdir/$chrom.agp") {
              $wormbase->run_command("cat $agp >> $target_agp_file", $log);
            } else {
              $log->error("ERROR: $gspecies : missing file: $chromdir/$chrom.agp\n");
            }
          }
          $wormbase->run_command("gzip -n -9 -f $target_agp_file", $log);
        }
      }
    }
  }
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying dna\n\n");
}


##################################################
# copy the GFF
##################################################
sub copy_gff_files{
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying gff files\n");

  my %accessors = ($wormbase->species_accessors);
  $accessors{elegans} = $wormbase;
  foreach my $wb (values %accessors) {
    next if exists $skip_species{$wb->species};
    next if @only_species and not exists $only_species{$wb->species};

    my $species  = $wb->species;
    my $gspecies = $wb->full_name('-g_species' => 1);
    my $bioproj  = $wb->ncbi_bioproject;
    my $sequencesdir = $wb->sequences;

    my $source_gff2_file = $wb->processed_GFF_file();
    my $source_gff3_file = $wb->processed_GFF_file(1);
    my $source_gtf_file = "${sequencesdir}/${species}.gtf";
    
    my $gff_dir = "$targetdir/species/$gspecies/$bioproj";
    mkpath($gff_dir,1,0775);
    my $fname_prefix = "$gspecies.$bioproj.$WS_version_name";
    
    my $target_gff2_file = "$gff_dir/${fname_prefix}.annotations.gff2.gz";
    my $target_gff3_file = "$gff_dir/${fname_prefix}.annotations.gff3.gz";
    my $target_gtf_file = "$gff_dir/${fname_prefix}.canonical_geneset.gtf.gz";

    foreach my $fname_pair ([$source_gff2_file,  $target_gff2_file],
                            [$source_gff3_file,  $target_gff3_file],
                            [$source_gtf_file,   $target_gtf_file]) {
      my ($source, $target) = @$fname_pair;
      unlink $target if -e $target;

      if (-e $source) {
        $wormbase->run_command("cat $source | gzip -n > $target", $log);
      } elsif (-e "${source}.gz") {
        $wormbase->run_command("cp -f ${source}.gz $target", $log);
      } else {
        $log->error("ERROR: Missing GFF file for $species ($source)\n");    
      }
    }
  }

  # copy tierIIIs from current build dir
  my %tierIII = $wormbase->tier3_species_accessors;
  foreach my $t3 (keys %tierIII){
    next if exists $skip_species{$t3};
    next if @only_species and not exists $only_species{$t3};

    my $wb = $tierIII{$t3};
    my $species  = $wb->species;
    my $gspecies = $wb->full_name('-g_species' => 1);
    my $bioproj  = $wb->ncbi_bioproject;
    my $gff_dir = "$targetdir/species/$gspecies/$bioproj";
    mkpath($gff_dir,1,0775);

    my $source_gff3 = sprintf("%s/%s.gff3", $wb->gff, $species);
    my $target_gff3 = sprintf("%s/%s.%s.%s.annotations.gff3.gz", $gff_dir, $gspecies, $bioproj, $WS_version_name);
    my $source_gtf =  sprintf("%s/%s.gtf", $wb->gff, $species);
    my $target_gtf = sprintf("%s/%s.%s.%s.canonical_geneset.gtf.gz", $gff_dir, $gspecies, $bioproj, $WS_version_name);

    foreach my $fname_pair ([$source_gff3, $target_gff3],
                            [$source_gtf, $target_gtf]) {
      my ($source, $target) = @$fname_pair;
      unlink $target if -e $target;

      if (-e $source) {
        $wb->run_command("cat $source | gzip -n -9 > $target", $log);
      } elsif (-e "${source}.gz") {
        $wb->run_command("cp -f ${source}.gz $target", $log);
      }
      # do not throw an error for missing GFF/GTF files here. This will be picked up 
      # later by the manifest check
    }
  }

  # Finally, check the GFF3 only to generate a comparison report
  my $report_file_name = "GFF3_comparison";
  my $collated_report_file = sprintf("%s/%s", $wormbase->reports, $report_file_name);
  my @all_report_files;

  my %all_accessors = $wormbase->all_species_accessors();
  foreach my $wb ($wormbase, sort { $a->species cmp $b->species } values %all_accessors) {
    my $g_species = $wb->full_name('-g_species' => 1);
    my $species = $wb->species;
    my $bioproj = $wb->ncbi_bioproject;

    next if exists $skip_species{$species};
    next if @only_species and not exists $only_species{$species};

    my $cur_ver_gff3_file = sprintf("%s/species/%s/%s/%s.%s.WS%s.annotations.gff3.gz", 
                                    $targetdir, 
                                    $g_species,
                                    $bioproj,
                                    $g_species,
                                    $bioproj,
                                    $WS_version);
    my $cur_ver_label = sprintf("%s.%s.WS%d", $g_species, $bioproj, $WS_version);

    my $prev_ver_gff3_file = sprintf("%s/releases/WS%d/species/%s/%s/%s.%s.WS%s.annotations.gff3.gz", 
                                     $wormbase->ftp_site,
                                     $WS_version - 1,
                                     $g_species,
                                     $bioproj,
                                     $g_species,
                                     $bioproj,
                                     $WS_version - 1);
    my $prev_ver_label = sprintf("%s.%s.WS%d", $g_species, $bioproj, $WS_version - 1);

    my $report_file = sprintf("%s/%s.%s", $wb->reports, $g_species, $report_file_name);

    if (-e $prev_ver_gff3_file and -e $cur_ver_gff3_file) {
      $wormbase->run_script("generate_gff_report.pl -final "
                            . "-currentgff $cur_ver_gff3_file -currentlabel $cur_ver_label "
                            . "-previousgff $prev_ver_gff3_file -previouslabel $prev_ver_label "
                            . "-report $report_file", $log);
      push @all_report_files, $report_file;
    } else {
      $log->write_to("Skipping GFF report for $species because could not find all necessary GFF files\n");
    }                  
  }
  if (@all_report_files) {
    $wormbase->run_command("cat @all_report_files > $collated_report_file", $log);
  }

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying GFF files\n\n");
}


###############################################
# copy across miscellaneous multi-species files
###############################################
sub copy_report_files {
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying report files\n");

  my $source_report_dir = $wormbase->reports;
  my $target_report_dir = "$targetdir/REPORTS";
  mkpath($target_report_dir, 1, 0775);

  foreach my $file ("GFF3_comparison", 
                    "all_classes") {
    if (-e "$source_report_dir/$file") {
      $wormbase->run_command("cat $source_report_dir/$file | gzip -n > $target_report_dir/${file}_report.$WS_version_name.txt.gz", $log); 
    }   
  }

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying report files\n\n");
}


###############################################
# copy across miscellaneous multi-species files
###############################################
sub copy_multi_species {

  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying clustal sql dump\n");

  my $tgt_dir = "$targetdir/MULTI_SPECIES";
  mkpath($tgt_dir,1,0775);
  
  my $chromdir = $wormbase->misc_output;

  my $clust_src = "$chromdir/wormpep_clw.sql.gz";

  if (-e $clust_src) {
    my $target = "$tgt_dir/wormpep_clw.${WS_version_name}.sql.gz";
    $wormbase->run_command("cp -f -R $clust_src $target", $log);

  } else {
    $log->write_to("Warning: no clustal results found ($clust_src)\n");
  }

  #
  # RRID strain data
  #
  my $RRID_source = $wormbase->acefiles . "/${WS_version_name}_RRIDs.dat";
  my $RRID_target = "$tgt_dir/strain_RRIDs.${WS_version_name}.dat.gz"; 
  if (not -e $RRID_source) {
    $log->error("ERROR: could not find $RRID_source file\n");
  } else {
    $wormbase->run_command("cat $RRID_source | gzip -c > $RRID_target", $log);
  }
  
  #
  # Disease data for ranjana
  #
  my $disease_source = $wormbase->acefiles . "/omim_db_data.ace";
  my $disease_target = "$tgt_dir/human_disease.${WS_version_name}.ace.gz";
  
  if (not -e $disease_source) {
    $log->error("ERROR: could not find $disease_source file\n");
  } else {
    $wormbase->run_command("cat $disease_source | gzip -c > $disease_target", $log);
  }

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying multi species files\n\n");
}


############################################
# copy across ncRNA files
############################################
sub copy_rna_files{
  # $wormbase is the main build object (likely elegans) ; $wb is the species specific one
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying rna files\n");

  # run through all possible organisms
  my %accessors = ($wormbase->species_accessors);
  $accessors{elegans} = $wormbase;

  foreach my $wb (values %accessors) {
    next if exists $skip_species{$wb->species};
    next if @only_species and not exists $only_species{$wb->species};

    my $gspecies = $wb->full_name('-g_species' => 1);
    my $bioproj = $wb->ncbi_bioproject;
    my $prefix = $wb->pepdir_prefix;
    my $rnadir = $wb->wormrna;

    my $ftprna_dir = "$targetdir/species/$gspecies/$bioproj";
    mkpath($ftprna_dir,1,0775);
	
    # Note that the following points to the most recent version of the file; 
    my $src_file = "$rnadir/${prefix}rna".$wb->get_wormbase_version.".rna"; 
    my $tgt_file = "$ftprna_dir/$gspecies.$bioproj.$WS_version_name.ncRNA_transcripts.fa"; #target FTP file
    
    if (-e $src_file and -s $src_file) {
      $wormbase->run_command("gzip -n -9 -c $src_file > ${tgt_file}.gz",$log);
    } else {
      # this is a core file, needs to be present for all species; however,
      # not all species have ncRNA data. Therefore, create an empty placeholder
      $wormbase->run_command("touch $tgt_file && gzip -n $tgt_file", $log);
    }   
  }

  # copy tierIIIs if present
  my %tierIII = $wormbase->tier3_species_accessors;
  foreach my $t3 (keys %tierIII){ 
    next if exists $skip_species{$t3};
    next if @only_species and not exists $only_species{$t3};

    my $wb = $tierIII{$t3};
    my $species  = $wb->species;
    my $gspecies = $wb->full_name('-g_species' => 1);
    my $bioproj = $wb->ncbi_bioproject;
    my $gff_dir = "$targetdir/species/$gspecies/$bioproj";
    mkpath($gff_dir,1,0775);

    my $target = sprintf("%s/%s.%s.%s.ncrna_transcripts.fa.gz", $gff_dir, $gspecies, $bioproj, $WS_version_name);
    my $unzipped_source = sprintf("%s/%s.ncrna.fa", $wb->sequences, $species);

    if (-e $unzipped_source) {
      $wb->run_command("gzip -n -9 -c $unzipped_source > $target", $log);
    } else {
      $log->write_to("Warning: no ncrna data found for species $species\n");
    }
  }

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying rna\n\n");
}


############################################
# copy across ontology files
############################################
sub copy_ontology_files {
  
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying ontology files\n");
  
  my $ace_dir = $wormbase->autoace;
  
  my $obo_dir = $wormbase->primaries . "/citace/temp_unpack_dir/home/citace/Data_for_${WS_version_name}/Data_for_Ontology/";
  my $ace_ontology_dir = "$ace_dir/ONTOLOGY";
  my $ftp_ontology_dir = "$targetdir/ONTOLOGY";
  
  mkpath($ace_ontology_dir,1,0775);
  mkpath($ftp_ontology_dir,1,0775);
  
  $wormbase->run_command("cp -f $obo_dir/*.obo $ace_ontology_dir/", $log);
  foreach my $file (glob("$ace_ontology_dir/*.*")) {
    $wormbase->run_command("cp -f $file $ftp_ontology_dir/", $log);
  }
  
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying ontology files\n\n");
}


############################################
# copy across annotation files
#############################################
sub copy_annotations_files{
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying files into per-species annotations folder...\n");

  my %accessors = ($wormbase->species_accessors);
  $accessors{$wormbase->species} = $wormbase; 

  my $gspecies = $wormbase->full_name('-g_species' => 1);
  my $bioproj = $wormbase->ncbi_bioproject;
  my $srcdir = $wormbase->misc_output;
  my $annotation_dir = "$targetdir/species/$gspecies/$bioproj/annotation";

  mkpath($annotation_dir,1,0775);

  #
  # Transcription Start Site wiggle files from the Julie Ahringer, Barbara Meyer and Tom Blumenthal papers 
  #
  my $source = $wormbase->misc_dynamic . "/c_elegans.PRJNA13758.WS240.TSS.wig.tar.gz";
  my $target = "$annotation_dir/${gspecies}.${bioproj}.${WS_version_name}.TSS.wig.tar.gz";
  
  if (-e $source) {
    $wormbase->run_command("cp -f $source $target", $log);
  } else {
    $log->write_to("Warning: TSS site wiggle plot tar file not found\n");
  }


  #
  # Reuters citation index - elegans only
  #
  $source = "$srcdir/ReutersCitationIndex.xml.gz";
  $target = "$annotation_dir/${gspecies}.${bioproj}.${WS_version_name}.reuters_citation_index.xml.gz";
  
  if (-e $source) {
    $wormbase->run_command("cp -f $source $target", $log);
  } else {
    $log->write_to("Warning: Reuters citation index not found\n");
  }

  #
  # OMIM xref data
  #
  $source = "$srcdir/${WS_version_name}_OMIMXREF.dat";
  $target = "$annotation_dir/${gspecies}.${bioproj}.${WS_version_name}.omim_xrefs.txt";
  
  if (-e $source) {
    $wormbase->run_command("cp -f $source $target", $log);
  } else {
    $log->write_to("Warning: OMIM xref file not found\n");
  }

  #
  # CGC changes data
  #
  $source = "$srcdir/changed_CGC_names.dat";
  $target = "$annotation_dir/${gspecies}.${bioproj}.${WS_version_name}.changed_CGC_names.txt";
  
  if (-e $source) {
    $wormbase->run_command("cp -f $source $target", $log);
  } else {
    $log->write_to("Warning: CGC changes file not found\n");
  }
  my $readable_source = "$srcdir/readable_changed_CGC_names.dat";
  my $readable_target = "$annotation_dir/${gspecies}.${bioproj}.${WS_version_name}.readable_changed_CGC_names.txt";
  
  if (-e $readable_source) {
    $wormbase->run_command("cp -f $readable_source $readable_target", $log);
  } else {
    $log->write_to("Warning: readable CGC changes file not found\n");
  }

  #
  # Oligo mapping
  #
  foreach my $file (glob("$srcdir/*oligo_mapping")) {
    my ($pre) = $file =~ /$srcdir\/(\S+_oligo_mapping)/;
    my $target = "$annotation_dir/$gspecies.$bioproj.$WS_version_name.${pre}.txt.gz";
    
    $wormbase->run_command("cat $file | gzip -n -9 -c > $target", $log);
  }

  #
  # TAR expression data for elegans
  #
  my $TARget = "$annotation_dir/$gspecies.$bioproj.$WS_version_name.TAR_gene_expression.tar.gz";
  my $TARexpr = $wormbase->spell . "/expr.tiling_arrays.tar.gz";
  if (-e $TARexpr) {
    $wormbase->run_command("cp -f $TARexpr $TARget", $log);
  } else {
    $log->write_to("Warning: gene expression file for $gspecies not found ($TARexpr)\n");
  }

  #
  # RNASeq gene expression data for each TierII and elegans
  #  
  foreach my $species (keys %accessors){
    next if exists $skip_species{$species};
    next if @only_species and not exists $only_species{$species};

    $log->write_to("copying $species expression data to FTP site\n");
    my $wb = $accessors{$species};

    my $g_species = $wb->full_name('-g_species' => 1);
    my $bioproj = $wb->ncbi_bioproject;

    my $src = $wb->spell . "/expr.rnaseq.tar.gz";
    my $tgt_dir = "$targetdir/species/$g_species/$bioproj/annotation/";

    if (-e $src) {
      mkpath($tgt_dir,1,0775);
      my $target = "$tgt_dir/$g_species.$bioproj.$WS_version_name.SRA_gene_expression.tar.gz";
      $wormbase->run_command("cp -f $src $target", $log);
    } else {
      $log->write_to("Warning: gene expression file for $g_species not found ($src)\n");
    }

    # controls FPKM file
    my $con_src = $wormbase->misc_dynamic . "/RNASeq_controls_FPKM_${species}.dat";

    if (-e $con_src) {
      mkpath($tgt_dir,1,0775);
      my $con_target = "$tgt_dir/$g_species.$bioproj.$WS_version_name.RNASeq_controls_FPKM.dat";
      $wormbase->run_command("cp -f $con_src $con_target", $log);
    } else {
      $log->write_to("Warning: controls gene expression file for $g_species not found ($con_src)\n");
    }
  }

  #
  # Xref files
  #
  foreach my $wb (values %accessors) {
    next if exists $skip_species{$wb->species};
    next if @only_species and not exists($only_species{$wb->species});

    my $gspecies = $wb->full_name('-g_species'=>1);
    my $bioproj = $wb->ncbi_bioproject;

    my $in_prefix = $wb->reports . "/" . $wb->species . ".";
    my $out_prefix = "$targetdir/species/$gspecies/$bioproj/annotation/$gspecies.${bioproj}.${WS_version_name}.";

    foreach my $io_pair ([ "dbxrefs.txt", "xrefs.txt.gz" ],
                         [ "gene_product_info.gpi", "gene_product_info.gpi.gz" ]) {
      my ($in_suffix, $out_suffix) = @$io_pair;

      my $in_file = $in_prefix . $in_suffix;
      my $out_file = $out_prefix . $out_suffix;

      if (-e $in_file) {
        $wormbase->run_command("cat $in_file | gzip -n -9 -c > $out_file", $log);
      }
    }
  }

  #
  # repeats
  #
  foreach my $wb (values %accessors) {
    next if exists $skip_species{$wb->species};
    next if @only_species and not exists($only_species{$wb->species});

    my $gspecies = $wb->full_name('-g_species'=>1);
    my $bioproj = $wb->ncbi_bioproject;

    my $file = $wb->repeatmasker_library;

    my $out_file = "$targetdir/species/$gspecies/$bioproj/annotation/$gspecies.${bioproj}.${WS_version_name}.repeats.fa.gz";
    if (-e $file) {
      $wormbase->run_command("cat $file | gzip -n -9 -c > $out_file", $log);
    } else {
      $log->error("ERROR: Could not find file $file\n");
    }
  }

  #
  # GO annotation files
  #
  foreach my $wb (values %accessors) {
    next if exists $skip_species{$wb->species};
    next if @only_species and not exists($only_species{$wb->species});

    my $gspecies = $wb->full_name('-g_species'=>1);
    my $bioproj = $wb->ncbi_bioproject;

    my $in_file = $wormbase->autoace . 'ONTOLOGY/gene_association.'. $WS_version_name.'.wb.' . $wb->species;
    my $out_file = "$targetdir/species/$gspecies/$bioproj/annotation/$gspecies.${bioproj}.${WS_version_name}.go_annotations.gaf.gz";

    if (-e $in_file) {
        $wormbase->run_command("cat $in_file | gzip -n -9 -c > $out_file", $log);
    }
  }


  #
  # Finally, other misc annotation files
  # 
  foreach my $wb (values %accessors) {
    next if exists $skip_species{$wb->species};
    next if @only_species and not exists($only_species{$wb->species});
    
    my $gspecies = $wb->full_name('-g_species'=>1);
    my $bioproj = $wb->ncbi_bioproject;
    
    my $in_prefix = $wb->reports;
    my $out_prefix = "$targetdir/species/$gspecies/$bioproj/annotation/$gspecies.${bioproj}.${WS_version_name}.";
    
    my @files = ('functional_descriptions.txt',
                 'orthologs.txt',
                 'protein_domains.csv');
    if ($wb->species eq 'elegans') {
      push @files, ('interactions.txt',
                    'potential_promotors.fa',
                    'swissprot.txt',
                    'knockout_consortium_alleles.xml',
                    'confirmed_genes.fa',
                    'cdna2orf.txt',
                    'geneIDs.txt',
                    'geneOtherIDs.txt',
                    'pcr_product2gene.txt',
                    'interpolated_clones.txt',
                    'molecules.ace'
                    );
    }
    
    
    foreach my $file (@files) {
      
      my $in_file = $in_prefix . "/" . $wb->species .".$file";
      my $out_file = $out_prefix . $file . '.gz';
      
      if (-e $in_file) {
        $log->write_to("WARNING: $in_file is empty\n") unless -s $in_file;
        $wormbase->run_command("cat $in_file | gzip -n -9 -c > $out_file", $log);
      } else {
        $log->write_to("WARNING: can't find $in_file\n");
      }
      
      $log->write_to("WARNING: $out_file is empty\n") unless -s $out_file;
    }
  }
  
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying misc files\n\n");
}


############################################
# copy across wormpep files
#############################################
sub copy_wormpep_files {
  
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying peptide files\n");
  
  my %accessors = ($wormbase->species_accessors);
  $accessors{$wormbase->species} = $wormbase; 
  
  foreach my $species (keys %accessors){
    next if exists $skip_species{$species};
    next if @only_species and not exists $only_species{$species};

    $log->write_to("copying $species protein data to FTP site\n");
    my $wb = $accessors{$species};

    my $gspecies = $wb->full_name('-g_species' => 1);
    my $bioproj = $wb->ncbi_bioproject;
    my $peppre = $wb->pepdir_prefix;
    my $rel_last_built = $wb->get_wormbase_version;

    my $src = $wb->wormpep;
    my $tgt = "$targetdir/species/$gspecies/$bioproj";
    mkpath($tgt,1,0775);

    if ($rel_last_built == $WS_version) {
      # tar up the latest wormpep release and copy across (files added in next loop)
      my $tgz_file = "$tgt/$gspecies.$bioproj.$WS_version_name.wormpep_package.tar.gz";
      my $command = "tar -c -z -h -f $tgz_file -C $src";
      
      # for tierIIs, the *pep package often does not change between 
      # releases. Given that it includes files for differences since last
      # release etc, forward propagation of the complete package and
      # incrementing the suffices is not appropriate. 
      
      foreach my $file ( $wb->wormpep_files ){
        if (-e "$src/$file${rel_last_built}") {
          $command .= " $file${rel_last_built}";
        }
      }
      $wb->run_command("$command", $log);
    }

    #single gzipped fasta files of peptides and transcripts

    my $source_pepfile = "$src/${peppre}pep${rel_last_built}";
    my $source_cdnafile = "$src/${peppre}pep.dna${rel_last_built}";

    my (%gene_ids);
    if (-e $source_pepfile) {
      my $target_pepfile = "$tgt/$gspecies.$bioproj.$WS_version_name.protein.fa.gz";
      open(my $source_pep_fh, $source_pepfile);
      while(<$source_pep_fh>) {
        /\>(\S+)\s+.*gene=(\S+)/ and do {
          $gene_ids{$1} = $2;
        };
      }
      close($source_pep_fh);

      $wb->run_command("gzip -n -9 -c $source_pepfile > $target_pepfile", $log);
    } else {
      $log->error("ERROR: Could not find peptide file for $gspecies ($source_pepfile)\n");
    }

    if (-e $source_cdnafile) {
      my $target_cdnafile = "$tgt/$gspecies.$bioproj.${WS_version_name}.CDS_transcripts.fa.gz";
      open(my $source_cdna_fh, $source_cdnafile);
      open(my $target_cdna_fh, "| gzip -n -9 -c > $target_cdnafile");
      while(<$source_cdna_fh>) {
        if (/^\>(\S+)(.*)/) {
          if (exists $gene_ids{$1}) {
            printf $target_cdna_fh ">%s gene=%s%s\n", $1, $gene_ids{$1}, $2;
          } else {
            printf $target_cdna_fh ">%s%s\n", $1, $2;
          }
        } elsif (/(\S+)/) {
          print $target_cdna_fh uc($1), "\n";
        }
      }
      close($target_cdna_fh) or $log->log_and_die("Could not successfully close $target_cdnafile\n");
    } else {
      $log->error("ERROR: Could not find transcript file for $gspecies ($source_cdnafile)\n");
    }
  }

  #continue with Tier III processing.

  #tierIII's have no "pep" package, so lift protein and transcript files from build dir
  # Note: not all tierIIIs will have protein/cDNA sets
  my %tierIII = $wormbase->tier3_species_accessors;
  foreach my $t3 (keys %tierIII){ 
    next if exists $skip_species{$t3};
    next if @only_species and not exists $only_species{$t3};

    my $wb = $tierIII{$t3};
    my $gspecies = $wb->full_name('-g_species' => 1);
    my $bioproj = $wb->ncbi_bioproject;

    my $tgt = "$targetdir/species/$gspecies/$bioproj";
    mkpath($tgt,1,0775);

    my $src_pep_file = $wb->sequences . "/$t3.prot.fa";
    my $src_cds_file = $wb->sequences . "/$t3.cds.fa";
    my $src_mrna_file = $wb->sequences . "/$t3.mrna.fa";

    my $tgt_pep_file = "$tgt/$gspecies.$bioproj.$WS_version_name.protein.fa.gz";
    my $tgt_cds_file = "$tgt/$gspecies.$bioproj.$WS_version_name.CDS_transcripts.fa.gz";
    my $tgt_mrna_file = "$tgt/$gspecies.$bioproj.$WS_version_name.mRNA_transcripts.fa.gz";

    if (-e $src_pep_file) {
      $wb->run_command("gzip -n -9 -c $src_pep_file > $tgt_pep_file", $log);
    } else {
      $log->write_to("Did not find $src_pep_file for $t3 - proceeding\n");
    }

    if (-e $src_cds_file) {
      $wb->run_command("gzip -n -9 -c $src_cds_file > $tgt_cds_file", $log);
    } else {
      $log->write_to("Did not find $src_cds_file for $t3 - proceeding\n");
    }

    if (-e $src_mrna_file) {
      $wb->run_command("gzip -n -9 -c $src_mrna_file > $tgt_mrna_file", $log);
    } else {
      $log->write_to("Did not find $src_mrna_file for $t3 - proceeding\n");
    }
  }

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying Pep files\n\n");
}


###################################
# copy_est_files
###################################
sub copy_est_files {
  
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying est files\n");
  
  my %accessors = ($wormbase->species_accessors);
  $accessors{$wormbase->species} = $wormbase; 

  foreach my $species (keys %accessors) {
    next if exists $skip_species{$species};
    next if @only_species and not exists $only_species{$species};

    my $wb = $accessors{$species};    
    my $gspecies = $wb->full_name('-g_species' => 1);
    my $bioproj = $wb->ncbi_bioproject;

    my $tdir = "$targetdir/species/$gspecies/$bioproj";
    mkpath("$tdir",1,0775);

    my $target = "$targetdir/species/$gspecies/$bioproj/$gspecies.$bioproj.$WS_version_name.ests.fa.gz";
    my $est_dir = $wb->build_data . "/cDNA/$species";

    open(my $outfh, "| gzip -n -c -9 > $target") 
        or $log->log_and_die("Could not open $target for writing\n");

    foreach my $type (qw(EST mRNA ncRNA OST RST)) {
      my $source = "$est_dir/$type";
      if (-e $source) {
        open(my $estfh, $source);
        while(<$estfh>) {
          if (/^\>(\S+)/) {
            printf $outfh ">%s  type=%s\n", $1, $type;
          } elsif (/^(\S+)$/) {
            print $outfh "$1\n";
          }
        }
      } else {
        $log->write_to("No $type data for $species\n");
      }
    }

    close($outfh);
    
  }

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copied est files\n");

}


###################################
# copy_assembly_manifest
###################################
sub copy_assembly_manifest {
  
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying assembly manifest file\n");

  my $ftp_acedb_dir = "$targetdir/species";
  mkpath("$ftp_acedb_dir",1,0775);
  my $tgt_manifile = "$ftp_acedb_dir/ASSEMBLIES.$WS_version_name.json";
  my $src_manifile = $wormbase->reports . "ASSEMBLIES.json";

  $log->log_and_die("Assembly manifest does not exist in REPORTS directory")
      if not -e $src_manifile;

  $wormbase->run_command("cp $src_manifile $tgt_manifile", $log);
  $log->write_to("$runtime: Copied assembly manifest file\n");
}



###################################
# copy_release_letter
###################################
sub copy_release_letter {
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying release letter\n");

  my $tgt_file = "$targetdir/letter.$WS_version_name";
  my $src_file = $wormbase->reports . "/letter.$WS_version_name";

  $log->log_and_die("Release letter does not exist in the REPORTS directory")
      if not -e $src_file;

  $wormbase->run_command("cp $src_file $tgt_file", $log);

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Copied release letter\n");
}



############################################################
# copy best blast hits file to ftp site
############################################################
sub copy_homol_data {

  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying homol files\n");

  my %accessors = ($wormbase->species_accessors);
  $accessors{elegans} = $wormbase;

  foreach my $wb (values %accessors) {
    next if exists $skip_species{$wb->species};
    next if @only_species and not exists $only_species{$wb->species};

    my $blast_dir = $wb->acefiles;
    my $species = $wb->species;
    my $gspecies = $wb->full_name('-g_species' => 1);
    my $bioproj = $wb->ncbi_bioproject;

    my $source_file = "$blast_dir/${species}_best_blastp_hits";
    # this script might be run more than once if there are problems
    if(-e $source_file || -e "$source_file.gz") { 
      my $protein_dir = "$targetdir/species/$gspecies/$bioproj";
      mkpath($protein_dir,1,0775);

      my $target_file = "$protein_dir/$gspecies.$bioproj.$WS_version_name.best_blastp_hits.txt.gz";
      # this script might be run more than once if there are problems
      $wormbase->run_command("gzip -n -9 -f $source_file",$log) if (! -e "$source_file.gz"); 
      $wormbase->run_command("scp $source_file.gz $target_file", $log);
    }
    $runtime = $wormbase->runtime;
    $log->write_to("$runtime: Finished copying\n\n");
  } 
}


##################################################################################
sub make_md5sums {

  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Calculating md5sums\n");

  my $checksum_file = "$targetdir/CHECKSUMS";

  my @files;
  open(FIND, "find $targetdir -name '*.*' |");
  while(<FIND>) {
    chomp;
    s/^$targetdir\///;
    push @files, $_;
  }

  $wormbase->run_command("cd $targetdir && md5sum @files > $checksum_file", $log);

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished making md5sums\n\n");

}


##################################################################################
sub CheckSize {
  my ($first,$second)=@_;
  my $F_SIZE = (stat("$first"))[7];
  my $S_SIZE = (stat("$second"))[7];
  if ($F_SIZE != $S_SIZE) {
    $log->error("ERROR: $first SRC: $F_SIZE TGT: $S_SIZE - not same size, please check\n");
  } 
}

##################################################################################
sub usage {
    system ('perldoc',$0);
}

##################################################################################
sub check_manifest {

  my %t2_species = $wormbase->species_accessors;
  $t2_species{$wormbase->species} = $wormbase;
  
  my %t3_species = $wormbase->tier3_species_accessors;
  
  my %species;
  foreach my $sp (keys %t3_species) {
    $species{$sp} = $t3_species{$sp};
  }  
  foreach my $sp (keys %t2_species) {
    $species{$sp} = $t2_species{$sp};
  }
  
  
  my (@stanzas, $count);
  
  while(<DATA>) {
    next unless /\w/;
    next if /^\#/;
    
    /^\[(.*)\](\S+)/ and do {
      my ($sp, $par) = ($1, $2);
      push @stanzas, {
        parent => $par,
        species => {},
      };
      my (%include, %exclude);
      if ($sp =~ /^(.*):(.*)$/) {
        my ($inc, $exc) = ($1, $2);
        map { $include{$_} = 1 } split(/,/, $inc);
        map { $exclude{$_} = 1 } split(/,/, $exc);
      } else {
        map { $include{$_} = 1 } split(/,/, $sp);
      }
      
      foreach my $tsp (keys %include) {
        if ($tsp eq 'TIER3') {
          foreach my $ntsp (keys %t3_species) {
            $stanzas[-1]->{species}->{$ntsp} = 1;
          }
        } elsif ($tsp eq 'CORE') {
          foreach my $ntsp (keys %t2_species) {
            $stanzas[-1]->{species}->{$ntsp} = 1;
          }
        } elsif ($tsp eq 'ALL') {
          foreach my $ntsp (keys %species) {
            $stanzas[-1]->{species}->{$ntsp} = 1;
          }
        } else {
          $stanzas[-1]->{species}->{$tsp} = 1;
        }
      }
      foreach my $tsp (keys %exclude) {
        if (exists $stanzas[-1]->{species}->{$tsp}) {
          delete $stanzas[-1]->{species}->{$tsp};
        }
      }
      next;
    };
    
    /(\S+)/ and do {
      push @{$stanzas[-1]->{files}}, $1;
    };
  }
  
  foreach my $block (@stanzas) {
    if (not keys %{$block->{species}}) {
      my $parent = $block->{parent};
      $parent = "$targetdir/$parent";
      foreach my $ofile (@{$block->{files}}) {
        my $file = $ofile;
        $file =~ s/WSREL/$WS_version_name/g;
        $file = "$parent/$file";
        $count += &checkfile($file);
      }
    } else {
      foreach my $worm (keys %species) {
        next if exists $skip_species{$worm};
        next if @only_species and not exists $only_species{$worm};
        
        next if (keys %{$block->{species}} and not exists $block->{species}->{$worm});
        
        my $gspecies = $species{$worm}->full_name('-g_species'=>1);
        my $bioproj = $species{$worm}->ncbi_bioproject;
        my $gWORM = $species{$worm}->pepdir_prefix;
        my $species = $species{$worm}->species;
        
        my $parent = $block->{parent};
        $parent =~ s/GSPECIES/$gspecies/g;
        $parent =~ s/BIOPROJ/$bioproj/g;
        $parent =~ s/WSREL/$WS_version_name/g;
        $parent = "$targetdir/$parent";
        
        foreach my $ofile (@{$block->{files}}) {
          my $file = $ofile;
          
          $file =~ s/GSPECIES/$gspecies/g;
          $file =~ s/BIOPROJ/$bioproj/g;
          $file =~ s/WSREL/$WS_version_name/g;
          $file =~ s/WORM/$gWORM/g;
          $file = "$parent/$file";
          $count += &checkfile($file);
        }
      }
    }
  }
  
  $log->write_to("$count files in place on FTP site\n");
}


sub checkfile {
  my $file = shift;
  if( -s "$file") {
    $log->write_to("OK: $file\n");
    return 1;
  } else {
    $log->error("ERROR: $file has a problem\n");
    return 0;
  }
}


__DATA__
[elegans]species/c_elegans/BIOPROJ/annotation
GSPECIES.BIOPROJ.WSREL.affy_oligo_mapping.txt.gz
GSPECIES.BIOPROJ.WSREL.agil_oligo_mapping.txt.gz
GSPECIES.BIOPROJ.WSREL.pcr_product2gene.txt.gz
GSPECIES.BIOPROJ.WSREL.gsc_oligo_mapping.txt.gz
GSPECIES.BIOPROJ.WSREL.cdna2orf.txt.gz
GSPECIES.BIOPROJ.WSREL.confirmed_genes.fa.gz
GSPECIES.BIOPROJ.WSREL.TAR_gene_expression.tar.gz
GSPECIES.BIOPROJ.WSREL.reuters_citation_index.xml.gz
GSPECIES.BIOPROJ.WSREL.omim_xrefs.txt
GSPECIES.BIOPROJ.WSREL.changed_CGC_names.txt
GSPECIES.BIOPROJ.WSREL.readable_changed_CGC_names.txt
GSPECIES.BIOPROJ.WSREL.TSS.wig.tar.gz
GSPECIES.BIOPROJ.WSREL.geneIDs.txt.gz
GSPECIES.BIOPROJ.WSREL.interactions.txt.gz
GSPECIES.BIOPROJ.WSREL.potential_promotors.fa.gz
GSPECIES.BIOPROJ.WSREL.swissprot.txt.gz
GSPECIES.BIOPROJ.WSREL.molecules.ace.gz

[elegans]species/GSPECIES/BIOPROJ
GSPECIES.BIOPROJ.WSREL.assembly.agp.gz
GSPECIES.BIOPROJ.WSREL.wormpep_package.tar.gz
GSPECIES.BIOPROJ.WSREL.transposon_transcripts.fa.gz
GSPECIES.BIOPROJ.WSREL.transposons.fa.gz

[CORE]species/GSPECIES/BIOPROJ/annotation
GSPECIES.BIOPROJ.WSREL.functional_descriptions.txt.gz
GSPECIES.BIOPROJ.WSREL.orthologs.txt.gz
GSPECIES.BIOPROJ.WSREL.protein_domains.csv.gz
GSPECIES.BIOPROJ.WSREL.go_annotations.gaf.gz
GSPECIES.BIOPROJ.WSREL.repeats.fa.gz
GSPECIES.BIOPROJ.WSREL.gene_product_info.gpi.gz


[CORE]species/GSPECIES/BIOPROJ
GSPECIES.BIOPROJ.WSREL.best_blastp_hits.txt.gz
GSPECIES.BIOPROJ.WSREL.annotations.gff2.gz
GSPECIES.BIOPROJ.WSREL.ncRNA_transcripts.fa.gz
GSPECIES.BIOPROJ.WSREL.pseudogenic_transcripts.fa.gz
GSPECIES.BIOPROJ.WSREL.intergenic_sequences.fa.gz
GSPECIES.BIOPROJ.WSREL.xrefs.txt.gz

[CORE:pristionchus]species/GSPECIES/BIOPROJ/annotation
GSPECIES.BIOPROJ.WSREL.SRA_gene_expression.tar.gz
GSPECIES.BIOPROJ.WSREL.RNASeq_controls_FPKM.dat

[ALL]species/GSPECIES/BIOPROJ
GSPECIES.BIOPROJ.WSREL.genomic.fa.gz
GSPECIES.BIOPROJ.WSREL.genomic_masked.fa.gz
GSPECIES.BIOPROJ.WSREL.genomic_softmasked.fa.gz
GSPECIES.BIOPROJ.WSREL.annotations.gff3.gz

[ALL:elegans_hawaii]species/GSPECIES/BIOPROJ
GSPECIES.BIOPROJ.WSREL.canonical_geneset.gtf.gz
GSPECIES.BIOPROJ.WSREL.protein.fa.gz
GSPECIES.BIOPROJ.WSREL.CDS_transcripts.fa.gz
GSPECIES.BIOPROJ.WSREL.mRNA_transcripts.fa.gz

[]acedb
files_in_tar
md5sum.WSREL
models.wrm.WSREL

[]species
ASSEMBLIES.WSREL.json

[]ONTOLOGY
anatomy_association.WSREL.wb
anatomy_ontology.WSREL.obo
gene_association.WSREL.wb.c_briggsae
gene_association.WSREL.wb.c_elegans
gene_association.WSREL.wb.c_remanei
gene_association.WSREL.wb.c_japonica
gene_association.WSREL.wb.c_brenneri
gene_association.WSREL.wb.b_malayi
gene_association.WSREL.wb.p_pacificus
gene_association.WSREL.wb.o_volvulus
gene_association.WSREL.wb.s_ratti
gene_association.WSREL.wb
gene_ontology.WSREL.obo
phenotype_association.WSREL.wb
phenotype_ontology.WSREL.obo
disease_association.WSREL.wb
disease_ontology.WSREL.obo
development_association.WSREL.wb
development_ontology.WSREL.obo

[]MULTI_SPECIES
wormpep_clw.WSREL.sql.gz
strain_RRIDs.WSREL.dat.gz
human_disease.WSREL.ace.gz
