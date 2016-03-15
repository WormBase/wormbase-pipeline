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

my ($help, $debug, $test, $testout, $verbose, $store, $wormbase);

my $acedb;   # only copy across acedb files
my $chroms;  # only copy across chromosome files
my $ont;     # only copy across ontology files
my $misc;    # only copy misc files
my $wormpep; # only copy wormpep files
my $genes;   # only copy confirmed genes
my $cDNA;    # only copy cDNA2orf file
my $geneIDs; # only copy file of gene IDs
my $pcr;     # only copy file of PCR products
my $homols;  # only copy best blast hits 
my $manifest;# check everything has been copied.
my $all;     # copy everything across
my $all_nopublic;
my $dna;
my $rna;
my $xrefs;
my $gff;
my $ests;
my $gbrowse_gff;
my $blastx;
my $dump_ko;
my $md5;
my $go_public;
my $compara;
my $orthology_lists; # flatfiles of the orthologies fro CalTech
my $assembly_manifest;
my (%skip_species, @skip_species, @only_species, %only_species, $WS_version, $WS_version_name);

GetOptions ("help"          => \$help,
	    "debug=s"       => \$debug,
	    "test"          => \$test,   # Copies data from the test env to a test ftp under ~/tmp/pub
            "testout=s"     => \$testout, # Copies real data to test location
	    "verbose"       => \$verbose,
	    "store:s"       => \$store,
	    "acedb"         => \$acedb,
            "blastx"        => \$blastx,
	    "compara"       => \$compara,
	    "dna"           => \$dna,
	    "rna"           => \$rna,
	    "wormpep"       => \$wormpep,
	    "gff"           => \$gff,
	    "misc"          => \$misc,
            "ests"          => \$ests,
	    "homols"        => \$homols,
            "xrefs"         => \$xrefs,
	    "ont"           => \$ont,
	    "genes"         => \$genes,
	    "cDNAlist"      => \$cDNA,
	    "geneIDs"       => \$geneIDs,
	    "pcr"           => \$pcr,
            "gbrowsegff"    => \$gbrowse_gff,
            "knockout"      => \$dump_ko,
            "assmanifest"   => \$assembly_manifest,
	    "manifest"      => \$manifest,
            "md5"            => \$md5,
            "public"         => \$go_public,
            "skipspecies=s@" => \@skip_species,
            "onlyspecies=s@" => \@only_species,
	    "all"            => \$all,
            "allnopublic"    => \$all_nopublic,
            "wbversion=s"    => \$WS_version,
            'orthologylists'=> \$orthology_lists,
    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# establish log file.
my $log = Log_files->make_build_log($wormbase);
map { $skip_species{$_} = 1 } @skip_species;
map { $only_species{$_} = 1 } @only_species;

if ($all) {
  $compara=$acedb=$dna=$gff=$rna=$misc=$wormpep=$genes=$cDNA=$ests=$geneIDs=$pcr=$homols=$manifest=$ont=$xrefs=$blastx=$dump_ko=$md5=$assembly_manifest=$go_public=$orthology_lists=1;
} elsif ($all_nopublic) {
  $compara=$acedb=$dna=$gff=$rna=$misc=$wormpep=$genes=$cDNA=$ests=$geneIDs=$pcr=$homols=$manifest=$ont=$xrefs=$blastx=$dump_ko=$md5=$assembly_manifest=$orthology_lists=1;
}

if (not $WS_version) {
  $WS_version = $wormbase->get_wormbase_version();
}
$WS_version_name = "WS${WS_version}";

my $maintainers     = "All";

my $targetdir = ($testout) 
    ? "$testout/releases/$WS_version_name"
    : $wormbase->ftp_site . "/staging/releases/$WS_version_name";

$log->write_to("WRITING TO $targetdir\n");

#################################################################################
# Main                                                                          #
#################################################################################

#if this fails its really important that the next step doesnt run - create lockat start and remove at end
# will be checked by finish_build.pl

my $lockfile = $wormbase->autoace . "/FTP_LOCK";
$log->write_to("writing lock file\n");

open (FTP_LOCK,">$lockfile") or $log->log_and_die("cant write lockfile $!\n");

print FTP_LOCK "If this file exists something has failed in make_ftp_site.pl\n DO NOT continue until you know what happend and have fixed it\n\n";
close FTP_LOCK;

&copy_acedb_files if ($acedb);    # make a new directory for the WS release and copy across release files

&copy_blastx if ($blastx);

&copy_compara if ($compara);          # copies the comparative data to the FTP site

&copy_dna_files if ($dna);

&copy_rna_files if ($rna);

&copy_wormpep_files if ($wormpep);    # copied from ~wormpub/WORMPEP

&copy_gff_files if ($gff);
  		  
&copy_est_files if ($ests);

&copy_misc_files if ($misc);

&copy_ontology_files if ($ont);       # make a new /ontology directory and copy files across

&copy_homol_data if ($homols);        # copies best blast hits files across

&copy_xrefs if ($xrefs);              # copies xref file for elegans and briggsae

&make_assembly_manifest if $assembly_manifest;

&extract_confirmed_genes if ($genes); # make file of confirmed genes from autoace and copy across

&make_cDNA2ORF_list if ($cDNA);       # make file of cDNA -> ORF connections and add to FTP site

&make_geneID_list if ($geneIDs);      # make file of WBGene IDs -> CGC name & Sequence name and add to FTP site

&make_pcr_list if ($pcr);             # make file of PCR products -> WBGene IDs, CDS, CGC name

&extract_ko if ($dump_ko);

&make_orthologylists if ($orthology_lists);

&make_gbrowse_gff if ($gbrowse_gff);

&check_manifest if ($manifest);       # compares whats on the FTP site with what should be

&make_md5sums if ($md5);              # creates a file of md5sums for containing entries for all files in release

$wormbase->run_command("chmod -R ug+w $targetdir", $log);  

if ($go_public) {
  if ($testout) {
    $log->write_to("You cannot go public having written to a test location; you must do it manually\n");
  } else {
    if ($log->report_errors == 0) {
      # moves release folder from staging to final resting place
      &go_public();
    }             
  }
}


# warn about errors in subject line if there were any
$wormbase->run_command("rm -f $lockfile",$log) if ($log->report_errors == 0);
$log->mail;
print "Finished.\n" if ($verbose);
exit (0);




#################################################################################
# Subroutines                                                                   #
#################################################################################


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

  opendir (RELEASE,"$ace_dir/release") or croak ("Could not open directory $ace_dir/release");
  while (defined($filename = readdir(RELEASE))) {
    if (($filename eq ".")||($filename eq "..")) { next;}
    if (($filename =~ /letter/)||($filename =~ /dbcomp/)) { next;}
    $wormbase->run_command("cp $ace_dir/release/$filename $ftp_acedb_dir/$filename", $log);

    my $O_SIZE = (stat("$ace_dir/release/$filename"))[7];
    my $N_SIZE = (stat("$ftp_acedb_dir/$filename"))[7];
    if ($O_SIZE != $N_SIZE) {
      $log->write_to("\tError: $filename SRC: $O_SIZE TGT: $N_SIZE - different file sizes, please check\n");
      croak "Couldn't copy $filename\n";
    } 
  }
  closedir RELEASE;
  
  # Copy across the models.wrm file
  $wormbase->run_command("cp $ace_dir/wspec/models.wrm $ftp_acedb_dir/models.wrm.$WS_version_name", $log);

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
# copy the xref files 
##################################################
sub copy_xrefs {
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying/zipping xref files\n");

  my %accessors = ($wormbase->species_accessors);
  $accessors{elegans} = $wormbase;

  foreach my $wb (values %accessors) {
    next if exists $skip_species{$wb->species};
    next if @only_species and not exists($only_species{$wb->species});

    my $gspecies = $wb->full_name('-g_species'=>1);
    my $bioproj = $wb->ncbi_bioproject;

    my $xref_file = $wb->acefiles . "/DBXREFs.txt";
    my $out_file = "$targetdir/species/$gspecies/$bioproj/$gspecies.${bioproj}.${WS_version_name}.xrefs.txt.gz";

    if (-e $xref_file) {
      $wormbase->run_command("cat $xref_file | gzip -n -9 -c > $out_file", $log);
    } else {
      $log->error("ERROR: Could not find xref file $xref_file\n");
    }
  }

  $log->write_to("$runtime: Finished copying/zipping xrefs\n\n");
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

    my $gff_dir = "$targetdir/species/$gspecies/$bioproj";
    mkpath($gff_dir,1,0775);
    my $fname_prefix = "$gspecies.$bioproj.$WS_version_name";

    my $target_gtf_file  = "$gff_dir/${fname_prefix}.canonical_geneset.gtf.gz";
    
    if (-e $source_gff2_file) {
      #concatenated whole genome file require for all species
      my $target = "$gff_dir/${fname_prefix}.annotations.gff2";
      $wormbase->run_command("rm -f $target", $log);
      $wormbase->run_command("cp -f -R $source_gff2_file $target", $log);
      $wormbase->run_command("gzip -n -9 -f $target",$log);
    } elsif (-e "${source_gff2_file}.gz") {
      my $target = "$gff_dir/${fname_prefix}.annotations.gff2.gz";
      $wormbase->run_command("cp -f -R ${source_gff2_file}.gz $target", $log);
    }else {
      $log->write_to("WARNING: No GFF2 file for $species ($source_gff2_file)\n");
    }
    
    if (-e $source_gff3_file) {
      #concatenated whole genome file require for all species
      my $target = "$gff_dir/${fname_prefix}.annotations.gff3";
      $wormbase->run_command("rm -f $target", $log);
      $wormbase->run_command("cp -f -R $source_gff3_file $target", $log);
      $wormbase->run_command("gzip -n -9 -f $target",$log);      

      $wb->run_command("perl $ENV{CVS_DIR}/GFF_post_process/extract_canonical_geneset.pl -infile $source_gff3_file -outfile $target_gtf_file", $log);

    } elsif (-e "${source_gff3_file}.gz") {
      my $target = "$gff_dir/${fname_prefix}.gff3.gz";
      $wormbase->run_command("cp -f -R ${source_gff3_file}.gz $target", $log);
      $wb->run_command("perl $ENV{CVS_DIR}/GFF_post_process/extract_canonical_geneset.pl -infile ${source_gff3_file}.gz -outfile $target_gtf_file", $log);
    } else {
      $log->error("ERROR: No GFF3 file for $species ($source_gff3_file)\n");
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

    my $source = sprintf("%s/%s.gff3", $wb->gff, $species);
    my $zipped_source = $source . ".gz";
    my $target = sprintf("%s/%s.%s.%s.annotations.gff3.gz", $gff_dir, $gspecies, $bioproj, $WS_version_name);

    if (-e $source) {
      $wb->run_command("cat $source | gzip -n -9 > $target", $log);
    } elsif (-e $zipped_source) {
      $wb->run_command("cp -f $zipped_source $target", $log);
    } else {
      $log->error("ERROR: No GFF3 data found for species $species\n");
    }
  }

  # Finally, check the GFF3 only to generate a comparison report
  my $report_file_name = sprintf("GFF3_comparison.WS%d-WS%d.txt", $WS_version - 1, $WS_version);
  my $collated_report_file = sprintf("%s/%s", $wormbase->reports, $report_file_name);
  my @all_report_files;

  my %all_accessors = $wormbase->all_species_accessors();
  foreach my $wb ($wormbase, sort { $a->species cmp $b->species } values %all_accessors) {
    my $g_species = $wb->full_name('-g_species' => 1);
    my $species = $wb->species;
    my $bioproj = $wb->ncbi_bioproject;

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

  my $target_report_dir = "$targetdir/REPORTS";
  mkpath($target_report_dir, 1, 0775);
  $wormbase->run_command("cp $collated_report_file $target_report_dir/", $log);


  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying GFF files\n\n");
}

############################################
# copy across comparative analysis
############################################

sub copy_compara {

  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying clustal sql dump\n");

  mkpath("$targetdir/COMPARATIVE_ANALYSIS",1,0775);
  
  my $chromdir = $wormbase->misc_output;

  my $clust_src = "$chromdir/wormpep_clw.sql.gz";

  if (-e $clust_src) {
    my $target = "$targetdir/COMPARATIVE_ANALYSIS/wormpep_clw.${WS_version_name}.sql.gz";
    $wormbase->run_command("cp -f -R $clust_src $target", $log);

  } else {
    $log->write_to("Warning: no clustal results found ($clust_src)\n");
  }


  my @compara_tar_args;

  my $acedir = $wormbase->acefiles;
  my $compara_source = "$acedir/compara.ace";
  if (not -e $compara_source) {
    $log->error("ERROR: could not find compara.ace file\n");
  } else {
    push @compara_tar_args, " -C $acedir compara.ace";
  }
  
  my $genomic_align_dir = $wormbase->misc_dynamic . "/SUPPLEMENTARY_GFF";
  my @files = glob("$genomic_align_dir/*.genomic_alignment.gff3");
  @files = map { $_ =~ /$genomic_align_dir\/(\S+)$/ and $1 } @files;
  if (@files) {
    push @compara_tar_args, " -C $genomic_align_dir @files";
  }

  if (@compara_tar_args) {
    my $target = "$targetdir/COMPARATIVE_ANALYSIS/compara.$WS_version_name.tar.gz";
    $wormbase->run_command("tar zcvf $target @compara_tar_args", $log);
  } else {
    $log->write_to("Warning: no compara results moved to FTP site\n");
  }

  # Copy hal file
  my $cactus=$wormbase->misc_static . '/core_worms.hal';
  if (not -e $cactus) {
    $log->error("ERROR: could not find $cactus file\n");
  } else {
    $wormbase->run_command("cp $cactus $targetdir/COMPARATIVE_ANALYSIS/");
  }

  # Copy RRID strain data
  my $RRID_source = "$acedir/${WS_version_name}_RRIDs.dat";
  if (not -e $RRID_source) {
    $log->error("ERROR: could not find $RRID_source file\n");
  } else {
    $wormbase->run_command("cp $RRID_source $targetdir/COMPARATIVE_ANALYSIS/");
  }
  
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying comparative results\n\n");
}

############################################
# create orthology lists
############################################

sub make_orthologylists{
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

       my $ftp_dir = "$targetdir/species/$gspecies/$bioproj/annotation";
       mkpath($ftp_dir,1,0775);

       my $ofile = "$ftp_dir/$gspecies.$bioproj.WS".$wormbase->get_wormbase_version.'.orthologs.txt';
       $wormbase->run_script("dump_ortholog_flatfile.pl -out $ofile -fullspecies \"${\$wb->full_name}\"",$log);
       if (-e $ofile and -s $ofile){
         $wormbase->run_command("touch $ofile && gzip -9 -n $ofile", $log);
       }
  }
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

sub copy_misc_files{
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying misc files\n");
  # zip and copy the microarray oligo mapping files.

  my $gspecies = $wormbase->full_name('-g_species' => 1);
  my $bioproj = $wormbase->ncbi_bioproject;
  my $srcdir = $wormbase->misc_output;
  my $annotation_dir = "$targetdir/species/$gspecies/$bioproj/annotation";

  mkpath($annotation_dir,1,0775);

  #
  # Transcription Start Site wiggle files from the Julie Ahringer, Barbara Meyer and Tom Blumenthal papers 
  #
  if ($wormbase->species eq 'elegans') {
    my $source = $wormbase->misc_dynamic . "/c_elegans.PRJNA13758.WS240.TSS.wig.tar.gz";
    my $target = "$annotation_dir/${gspecies}.${bioproj}.${WS_version_name}.TSS.wig.tar.gz";

    if (-e $source) {
      $wormbase->run_command("cp -f $source $target", $log);
    } else {
      $log->write_to("Warning: TSS site wiggle plot tar file not found\n");
    }
  }


  #
  # Reuters citation index - elegans only
  #
  if ($wormbase->species eq 'elegans') {
    my $source = "$srcdir/ReutersCitationIndex.xml.gz";
    my $target = "$annotation_dir/${gspecies}.${bioproj}.${WS_version_name}.reuters_citation_index.xml.gz";

    if (-e $source) {
      $wormbase->run_command("cp -f $source $target", $log);
    } else {
      $log->write_to("Warning: Reuters citation index not found\n");
    }
  }

  #
  # OMIM xref data
  #
  if ($wormbase->species eq 'elegans') {
    my $source = "$srcdir/${WS_version_name}_OMIMXREF.dat";
    my $target = "$annotation_dir/${gspecies}.${bioproj}.${WS_version_name}.omim_xrefs.txt";
    
    if (-e $source) {
      $wormbase->run_command("cp -f $source $target", $log);
    } else {
      $log->write_to("Warning: OMIM xref file not found\n");
    }
  }

  #
  # CGC changes data
  #
  if ($wormbase->species eq 'elegans') {
    my $source = "$srcdir/changed_CGC_names.dat";
    my $target = "$annotation_dir/${gspecies}.${bioproj}.${WS_version_name}.changed_CGC_names.txt";
    
    if (-e $source) {
      $wormbase->run_command("cp -f $source $target", $log);
    } else {
      $log->write_to("Warning: CGC changes file not found\n");
    }
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
  # Disease data for ranjana
  #
  # make the specific places for the files to be copied, probably overkill do create the local DISEASE dir?
  my $ace_dir = $wormbase->autoace;
  my $ace_disease_dir = "$ace_dir/DISEASE";
  my $ftp_disease_dir = "$targetdir/DISEASE";

  mkpath($ace_disease_dir,1,0775);
  mkpath($ftp_disease_dir,1,0775);
  my $Disease_source = $wormbase->acefiles."/omim_db_data.ace";
  my $Disease_targ = "Human_disease_data";
  
  if (-e $Disease_source) {
    $wormbase->run_command("cp -f $Disease_source $ace_disease_dir/$Disease_targ", $log);
    $wormbase->run_command("cp -f $Disease_source $ftp_disease_dir/$Disease_targ", $log);
  }
  else {
    $log->write_to("Warning: Disease data file for $gspecies not found ($Disease_source)\n");
  }

  #
  # RNASeq gene expression data for each TierII and elegans
  #
  my %accessors = ($wormbase->species_accessors);
  $accessors{$wormbase->species} = $wormbase; 
  
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
        /^\>(\S+)\s+\S+\s+(\S+)/ and do {
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
          print $target_cdna_fh "$1\n";
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

    my $src_pep_file = $wb->sequences . "/$t3.prot.fas";
    my $src_cdna_file = $wb->sequences . "/$t3.cdna.fas";

    my $tgt_pep_file = "$tgt/$gspecies.$bioproj.$WS_version_name.protein.fa.gz";
    my $tgt_cdna_file = "$tgt/$gspecies.$bioproj.$WS_version_name.CDS_transcripts.fa.gz";

    if (-e $src_pep_file) {
      $wb->run_command("gzip -n -9 -c $src_pep_file > $tgt_pep_file", $log);
    } else {
      $log->write_to("Did not find $src_pep_file for $t3 - proceeding\n");
    }

    if (-e $src_cdna_file) {
      $wb->run_command("gzip -n -9 -c $src_cdna_file > $tgt_cdna_file", $log);
    } else {
      $log->write_to("Did not find $src_cdna_file for $t3 - proceeding\n");
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
        or croak("Could not open $target for writing\n");

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
# make_assembly_manifest
###################################
sub make_assembly_manifest {
  
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: making assembly manifest file\n");

  my $ftp_acedb_dir = "$targetdir/species";
  mkpath("$ftp_acedb_dir",1,0775);
  my $manifile = "$ftp_acedb_dir/ASSEMBLIES.$WS_version_name.json";

  open(my $fh, ">$manifile\n") or $log->log_and_die("Could not open $manifile for writing\n");
  my $db = Ace->connect(-path => $wormbase->autoace );

  my (%accessors, %accessors_by_species, %json);

  %accessors = ($wormbase->all_species_accessors);
  $accessors{$wormbase->species} = $wormbase;
  
  foreach my $acc (values %accessors) {
    my $species_name = $acc->full_name;
    push @{$accessors_by_species{$species_name}}, $acc;
  }
  
  foreach my $species (sort keys %accessors_by_species) {
    my (@accs) = @{$accessors_by_species{$species}};
    
    my $species_obj = $db->fetch(-class => 'Species', -name => "$species");
    my @seq_col = $species_obj->at('Assembly');
    
    my $g_species = $accs[0]->full_name(-g_species => 1); 
    
    my $obj = {
      full_name => $species,
      assemblies => [],
    };
    
    $json{$g_species} = $obj;
    
    foreach my $seq_col_name (@seq_col) {
      my $seq_col = $seq_col_name->fetch;
      
      my $strain = $seq_col->Strain;
      my $assembly_name = $seq_col->Name;
      unless (defined $assembly_name) {
	$log->write_to("ERROR: The $g_species assembly appears to be unnamed so the JSON file will reflect this....this DOES have implications for the website, so best to check and patch the build if there is an official assembly name and you have time (a patch can be supplied to the web team so you don't have to re-package from scratch).\n");
      }
      my $first_ws_rel = $seq_col->First_WS_release;
      my @labs;
      
      if ($seq_col->at('Origin.Laboratory')) {      
        my @laboratory = $seq_col->at('Origin.Laboratory');
        foreach my $lab (@laboratory) {
          push @labs, $lab->fetch->Mail->name;
        }
      }
      
      my ($bioproj, $gc_acc);
      
      my @db = $seq_col->at('Origin.DB_info.Database');
      foreach my $db (@db) {
        if ($db->name eq 'NCBI_BioProject') {
          $bioproj = $db->right->right->name;
        } elsif ($db->name eq 'NCBI_Genome_Assembly') {
          $gc_acc  = $db->right->right->name;
        }
      }
      
      # need to find the corresponding accessor, because only that
      # will tell us whether the bioproject is the canonical one
      my ($rel_acc) = grep { $_->ncbi_bioproject eq $bioproj } @accs;
      next unless $rel_acc;
      my $is_canonical = $rel_acc->is_canonical;
      my $bioproj_desc = $rel_acc->bioproject_description;

      push @{$obj->{assemblies}}, {
        bioproject => $bioproj,
        bioproject_description => $bioproj_desc,
        assembly_accession => $gc_acc,
        assembly_name => (defined $assembly_name) ? $assembly_name->name : "$g_species Assembly unnamed",
        appeared_in => 'WS'.$first_ws_rel->name,
        is_canonical => ($is_canonical) ? JSON::true : JSON::false,
        strain => (defined $strain) ? $strain->name : "Unknown strain",
        laboratory => \@labs,
      };
    }
  }
  
  $db->close;

  my $json_obj = JSON->new;
  my $string = $json_obj->allow_nonref->canonical->pretty->encode(\%json);
  
  print $fh $string;
  close($fh);
  

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: written assembly manifest file\n");

}


###################################
# extract confirmed genes
###################################

sub extract_confirmed_genes{

  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Extracting elegans confirmed genes\n");

  my $gspecies = $wormbase->full_name('-g_species' => 1);
  my $bioproj = $wormbase->ncbi_bioproject;

  my $ace_dir = $wormbase->autoace;
  my $annotation_dir = "$targetdir/species/$gspecies/$bioproj/annotation";
  my $wormdna = $wormbase->wormpep."/".$wormbase->pepdir_prefix."pep.dna${WS_version}";
  my $seqio  = Bio::SeqIO->new('-file' => "$wormdna", '-format' => 'fasta');

  my $db = Ace->connect(-path => "$ace_dir/") || die print "ACE connection failure: ", Ace->error;
  my $query = "Find elegans_CDS; Confirmed";

  my @confirmed_genes   = $db->fetch(-query=>$query);
  my %conf_genes;
  map { $conf_genes{$_->name} = 1 } @confirmed_genes;

  mkpath("$annotation_dir",1,0775);

  my $out_file = "$annotation_dir/$gspecies.$bioproj.$WS_version_name.confirmed_genes.fa";

  open(OUT,">$out_file") or 
      $log->write_to("Couldn't write to $out_file\n");

  while(my $seq = $seqio->next_seq){
    if($conf_genes{$seq->id}) {
      print OUT "\n>".$seq->id."\n";
      print OUT $wormbase->format_sequence($seq->seq);
    }
  }
  close(OUT);

  $wormbase->run_command("gzip -n -9 -f $out_file", $log);

  $db->close;

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished extracting confirmed genes\n\n");

  return(0);
}


################################################################################
# dump KO data
################################################################################

sub extract_ko {
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Extracting Knockout Consortium Data\n");

  my $gspecies = $wormbase->full_name('-g_species' => 1);
  my $bioproj = $wormbase->ncbi_bioproject;

  my $annotation_dir = "$targetdir/species/$gspecies/$bioproj/annotation";
  mkpath("$annotation_dir",1,0775);

  my $outfile = "$annotation_dir/${gspecies}.${bioproj}.${WS_version_name}.knockout_consortium_alleles.xml.gz";

  $wormbase->run_script("dump_ko.pl -file $outfile",$log);
  
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished dumping KO data\n\n");

}


################################################################################
# make list of cDNA -> orf connections
################################################################################

sub make_cDNA2ORF_list {

  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: making cDNA2ORF files\n");
  # simple routine to just get cDNA est names and their correct ORFs and make an FTP site file
  # two columns, second column supports multiple ORF names

  my $gspecies = $wormbase->full_name('-g_species' => 1);
  my $bioproj = $wormbase->ncbi_bioproject;
  my $annotation_dir = "$targetdir/species/$gspecies/$bioproj/annotation";
  my $ace_dir = $wormbase->autoace;

  mkpath("$annotation_dir",1,0775);

  my $tace = $wormbase->tace;
  my $command=<<EOF;
Table-maker -p "$ace_dir/wquery/cDNA2CDS.def" quit
EOF

  my $dir = "$ace_dir";

  my %cDNA2orf;
  open (TACE, "echo '$command' | $tace $dir | ") || croak "Couldn't access $dir\n";  
  while (<TACE>){
    chomp;
    if (m/^\"/){
      s/\"//g;
      m/(.+)\s+(.+)/;     
      $cDNA2orf{$1} .= "$2 ";
    }
  }
  close(TACE);
  
  # output to ftp site

  
  my $out = "$annotation_dir/$gspecies.$bioproj.$WS_version_name.cdna2orf.txt";
  open(OUT, ">$out") or $log->write_to("Couldn't open $out\n");

  foreach my $key (keys %cDNA2orf){
    print OUT "$key,$cDNA2orf{$key}\n";
  }
  close(OUT);

  $wormbase->run_command("gzip -n -9 -f $out", $log);

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished making cdna2orf list\n\n");  
  
}

################################################################################
# make list of WBGene IDs to CGC name and Sequence name
################################################################################

sub make_geneID_list {
  # For each 'live' Gene object, extract 'CGC_name' and 'Sequence_name' fields (if present)

  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: making Gene ID list\n");
 
  my %accessors = ($wormbase->species_accessors);
  foreach my $wb ($wormbase, values %accessors) {
    next if exists $skip_species{$wb->species};
    next if @only_species and not exists $only_species{$wb->species};

    my $gspecies = $wb->full_name('-g_species' => 1);
    my $bioproj = $wb->ncbi_bioproject;
    my $full_name = $wb->full_name;
    my $tax_id = $wb->ncbi_tax_id;

    my $annotation_dir = "$targetdir/species/$gspecies/$bioproj/annotation";
    my $ace_dir = $wormbase->autoace;

    my $out     = "$annotation_dir/$gspecies.$bioproj.$WS_version_name.geneIDs.txt";
    my $otherNames = "$annotation_dir/$gspecies.$bioproj.$WS_version_name.geneOtherIDs.txt";

    mkpath("$annotation_dir",1,0775);
    open GENEID,">$out" ||die($!);
    open OTHERID, ">$otherNames" ||die($!);

    my $db = Ace->connect(-path => "$ace_dir/") || die (Ace->error);
    my $gene_it = $db->fetch_many(-query => "Find Gene; Species=\"${full_name}\"");
    while(my $gene=$gene_it->next){
      print GENEID join(",", 
                         $tax_id,
                         $gene,
                         ($gene->CGC_name||''),
                         ($gene->Sequence_name||''),
                         $gene->Status), "\n";
      print OTHERID join("\t",$gene,$gene->Status,$gene->Sequence_name,$gene->CGC_name,$gene->Other_name),"\n";
    }

    close GENEID;
    close OTHERID;
    $db->close();
    $wormbase->run_command("gzip -n -9 -f $out", $log);
    $wormbase->run_command("gzip -n -9 -f $otherNames", $log);
  } 

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished making gene_id list\n\n");
}


################################################################################
# make list of PCR_product connections to CDS and Gene ID plus CGC name
################################################################################

sub make_pcr_list {

  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: making PCR product 2 gene list list\n");

  my $gspecies = $wormbase->full_name('-g_species' => 1);
  my $bioproj = $wormbase->ncbi_bioproject;
  my $annotation_dir = "$targetdir/species/$gspecies/$bioproj/annotation";
  my $ace_dir = $wormbase->autoace;

  my $tace    = $wormbase->tace;
  my $command = "Table-maker -p $ace_dir/wquery/pcr_product2gene.def\nquit\n";
  my $dir     = "$ace_dir";
  my $out     = "$annotation_dir/$gspecies.$bioproj.$WS_version_name.pcr_product2gene.txt";
  mkpath("$annotation_dir",1,0775);

  # hashes needed because one pcr product may hit two or more genes
  my %pcr2gene;
  my %gene2sequence;
  my %gene2cgc;

  open (TACE, "echo '$command' | $tace $dir | ") || croak "Couldn't access $dir\n";  
  while (<TACE>){
    if (m/^\"/){
      s/\"//g;
      chomp;
      # split the line into various fields
      my ($pcr,$gene,$cgc,$sequence) = split(/\t/, $_) ;

      # fill hashes and arrays, allow for multiple genes per PCR product
      $pcr2gene{$pcr}      .= "$gene,";
      ($gene2cgc{$gene} = $cgc) if($cgc);
      $gene2sequence{$gene} = "$sequence";
    }
  }
  close(TACE);

  # Now write output, cycling through list of pcr products in %pcr2gene

  open (OUT, ">$out") or $log->log_and_die("Couldn't open $out\n");

  foreach my $pcr (keys %pcr2gene){
    my @genes = split(/,/,$pcr2gene{$pcr});
    my $counter =0;
    print OUT "$pcr";
    foreach my $gene (@genes){

      # remove next element if it is the same gene ID and start loop again
      if (defined($genes[$counter+1]) && $genes[$counter+1] eq $genes[$counter]){
	splice(@genes, $counter+1,1);
	redo;
      }
      $counter++;
      print OUT "\t$gene";
      print OUT "($gene2cgc{$gene})" if (exists($gene2cgc{$gene}));
      
      # now print sequence name
      print OUT ",$gene2sequence{$gene}";      
    }
    print OUT "\n";
  }

  close(OUT);

  $wormbase->run_command("gzip -n -9 -f $out", $log);

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished making list\n\n");
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

##############################################################
sub make_gbrowse_gff {
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Making Gbrowse GFF\n");

  my %accessors = ($wormbase->species_accessors);
  $accessors{elegans} = $wormbase;
  
  foreach my $wb (values %accessors) {
    my $species = $wb->species;

    next if exists $skip_species{$species};
    next if @only_species and not exists $only_species{$species};

    my $gspecies = $wb->full_name('-g_species' => 1);
    my $bioproj = $wb->ncbi_bioproject;
    my $tgt_dir = "$targetdir/species/$gspecies/$bioproj";

    my $in_filename = "$tgt_dir/$gspecies.$bioproj.$WS_version_name.annotations.gff2.gz"; 

    if (not -e $in_filename) {
      $log->error("ERROR: Could not make GBrowse-ready GFF for $species; in file $in_filename not found\n");
      next;
    }

    my $out_filename = "$tgt_dir/$gspecies.$bioproj.$WS_version_name.GBrowse.gff2";
    my $command = "perl $ENV{CVS_DIR}/web_data/process_elegans_gff-standalone.pl -species $species -outfile $out_filename";
    if ($debug) {
      $command .= " -debug $debug ";
    }

    $wb->run_command("gunzip -c $in_filename | $command", $log);
    $wb->run_command("gzip -n -9 $out_filename", $log);
  }

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished generating Gbrowse GFF\n\n");
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
sub go_public {
  my $base_path = $wormbase->ftp_site;
  my $staging_path = "$base_path/staging/releases/$WS_version_name";
  my $final_path   = "$base_path/releases/$WS_version_name";
  my $frozen_path = "$base_path/FROZEN_RELEASES";

  $wormbase->run_command("mv $staging_path $final_path", $log);

  if ($WS_version =~ /\d+5$/ || $WS_version =~ /\d+0$/) {
    $wormbase->run_command("cd $frozen_path && ln -sf ../releases/$WS_version_name $WS_version_name", $log);
  }
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
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
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

#for manifest checks
# gspecies = $wb->full_name('-gspecies' => 1)
# WORM = $wb->pepdir_prefix

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
GSPECIES.BIOPROJ.WSREL.TSS.wig.tar.gz

[elegans]species/GSPECIES/BIOPROJ
GSPECIES.BIOPROJ.WSREL.assembly.agp.gz
GSPECIES.BIOPROJ.WSREL.wormpep_package.tar.gz
GSPECIES.BIOPROJ.WSREL.transposon_transcripts.fa.gz

[CORE]species/GSPECIES/BIOPROJ/annotation
GSPECIES.BIOPROJ.WSREL.geneIDs.txt.gz
GSPECIES.BIOPROJ.WSREL.orthologs.txt.gz

[CORE]species/GSPECIES/BIOPROJ
GSPECIES.BIOPROJ.WSREL.best_blastp_hits.txt.gz
GSPECIES.BIOPROJ.WSREL.annotations.gff2.gz
GSPECIES.BIOPROJ.WSREL.ncRNA_transcripts.fa.gz
GSPECIES.BIOPROJ.WSREL.mRNA_transcripts.fa.gz
GSPECIES.BIOPROJ.WSREL.pseudogenic_transcripts.fa.gz
GSPECIES.BIOPROJ.WSREL.intergenic_sequences.fa.gz
GSPECIES.BIOPROJ.WSREL.canonical_geneset.gtf.gz
GSPECIES.BIOPROJ.WSREL.xrefs.txt.gz

[CORE:pristionchus]species/GSPECIES/BIOPROJ/annotation
GSPECIES.BIOPROJ.WSREL.SRA_gene_expression.tar.gz
GSPECIES.BIOPROJ.WSREL.RNASeq_controls_FPKM.dat

# for all species
[ALL]species/GSPECIES/BIOPROJ
GSPECIES.BIOPROJ.WSREL.genomic.fa.gz
GSPECIES.BIOPROJ.WSREL.genomic_masked.fa.gz
GSPECIES.BIOPROJ.WSREL.genomic_softmasked.fa.gz
GSPECIES.BIOPROJ.WSREL.annotations.gff3.gz
GSPECIES.BIOPROJ.WSREL.protein.fa.gz
GSPECIES.BIOPROJ.WSREL.CDS_transcripts.fa.gz


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

[]COMPARATIVE_ANALYSIS
compara.WSREL.tar.gz
wormpep_clw.WSREL.sql.gz
core_worms.hal
WSREL_RRIDs.dat
