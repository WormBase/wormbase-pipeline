#!/usr/local/bin/perl5.8.0 -w
#
# make_FTP_sites.pl
#
# A PERL wrapper to automate the process of building the FTP sites 
# builds wormbase & wormpep FTP sites
# 
# Last updated by: $Author: klh $
# Last updated on: $Date: 2011-06-01 15:37:20 $
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


#################################################################################
# Command-line options and variables                                            #
#################################################################################

my ($help, $debug, $test, $testout, $verbose, $store, $wormbase);

my $release; # only copy across release files
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
my $dna;
my $rna;
my $gff;
my $ests;
my $gbrowse_gff;
my $blastx;
my $dump_ko;
my $md5;
my $supplementary;
my $compara;
my (%skip_species, @skip_species, @only_species, %only_species);

GetOptions ("help"          => \$help,
	    "debug=s"       => \$debug,
	    "test"          => \$test,   # Copies data from the test env to a test ftp under ~/tmp/pub
            "testout=s"     => \$testout, # Copies real data to test location
	    "verbose"       => \$verbose,
	    "store:s"       => \$store,
	    "release"       => \$release,
            "blastx"         => \$blastx,
	    "compara"       => \$compara,
	    "dna"           => \$dna,
	    "rna"           => \$rna,
	    "wormpep"       => \$wormpep,
	    "gff"           => \$gff,
	    "supplementary" => \$supplementary,
	    "misc"          => \$misc,
            "ests"          => \$ests,
	    "homols"        => \$homols,
	    "ont"           => \$ont,
	    "genes"         => \$genes,
	    "cDNAlist"      => \$cDNA,
	    "geneIDs"       => \$geneIDs,
	    "pcr"           => \$pcr,
            "gbrowsegff"    => \$gbrowse_gff,
            "knockout"      => \$dump_ko,
	    "manifest"      => \$manifest,
            "md5"            => \$md5,
            "skipspecies=s@" => \@skip_species,
            "onlyspecies=s@" => \@only_species,
	    "all"            => \$all);


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

# using -all option?
#($clustal=$release=$dump_ko=$dna=$gff=$supplementary=$rna=$misc=$wormpep=$genes=$cDNA=$geneIDs=$pcr=$homols=$manifest=$ont = 1 ) if ($all);
# remove supplementary from above list, now off by default
($compara=$release=$dna=$gff=$rna=$misc=$wormpep=$genes=$cDNA=$ests=$geneIDs=$pcr=$homols=$manifest=$ont=$blastx=$dump_ko=$gbrowse_gff=$md5=1 ) if ($all);

my $wormpep_ftp_root = ($testout) 
    ? "$testout/databases/wormpep"
    : "/nfs/disk69/ftp/pub/databases/wormpep";

my $WS              = $wormbase->get_wormbase_version();      # e.g.   132
my $WS_name         = $wormbase->get_wormbase_version_name(); # e.g. WS132
my $maintainers     = "All";

my $targetdir = ($testout) 
    ? "$testout/releases/$WS_name"
    : $wormbase->ftp_site . "/releases/$WS_name";

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

&copy_release_files if ($release);    # make a new directory for the WS release and copy across release files

&copy_blastx if ($blastx);

&copy_compara if ($compara);          # copies the comparative data to the FTP site

&copy_dna_files if ($dna);

&copy_rna_files if ($rna);

&copy_wormpep_files if ($wormpep);    # copied from ~wormpub/WORMPEP

&copy_gff_files if ($gff);
  		  
&copy_supplementary_gff_files if ($supplementary);  		  

&copy_est_files if ($ests);

&copy_misc_files if ($misc);

&copy_ontology_files if ($ont);       # make a new /ontology directory and copy files across

&copy_homol_data if ($homols);        # copies best blast hits files across

&extract_confirmed_genes if ($genes); # make file of confirmed genes from autoace and copy across

&make_cDNA2ORF_list if ($cDNA);       # make file of cDNA -> ORF connections and add to FTP site

&make_geneID_list if ($geneIDs);      # make file of WBGene IDs -> CGC name & Sequence name and add to FTP site

&make_pcr_list if ($pcr);             # make file of PCR products -> WBGene IDs, CDS, CGC name

&extract_ko if ($dump_ko);

&make_gbrowse_gff if ($gbrowse_gff);

&check_manifest if ($manifest);       # compares whats on the FTP site with what should be

&make_md5sums if ($md5);              # creates a file of md5sums for containing entries for all files in release

&make_frozen_links;                   # creates a frozen symbolic link if it is a 5 or 0 release



################################
#
# Tidy up and exit
#
################################
$wormbase->run_command("chgrp -R  worm $targetdir", $log);  
$wormbase->run_command("chmod -R ug+w $targetdir", $log);  

# warn about errors in subject line if there were any
$wormbase->run_command("rm -f $lockfile",$log) if ($log->report_errors == 0);
$log->mail;
print "Finished.\n" if ($verbose);
exit (0);




#################################################################################
# Subroutines                                                                   #
#################################################################################


##########################################################
# copy the WS release files across and check on the size
# The FTP disk tends to be unstable
##########################################################

sub copy_release_files{
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
  $wormbase->run_command("cp $ace_dir/wspec/models.wrm $ftp_acedb_dir/models.wrm.$WS_name", $log);

  # copy some miscellaneous files across
  my $old_release = $WS -1;
  $wormbase->run_command("cp ".	$wormbase->compare."/WS$old_release-$WS_name.dbcomp $ftp_acedb_dir", $log);
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
      my $out_file = "$blastx_dir/$gspecies.$WS_name.blastx.ace.gz";

      $wormbase->run_command("cat $in_file | gzip -9 -c > $out_file", $log);
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

  foreach my $wb (values %accessors) {
    next if exists $skip_species{$wb->species};
    next if @only_species and not exists $only_species{$wb->species};

    my $gspecies = $wb->full_name('-g_species'=>1);
    my $chromdir = $wb->chromosomes;
    
    my %copied_files;

    if (-e "$chromdir") {
      my $dna_dir = "$targetdir/species/$gspecies";
      mkpath($dna_dir,1,0775);
      #todd wants all species to have whole genome in one file
      if ($wb->assembly_type eq 'contig') {
	my $species = $wb->species;

	my $dna_file = "$chromdir/supercontigs.fa";
	my $masked_file = "$chromdir/".$species."_masked.dna.gz";
	my $soft_file = "$chromdir/".$species."_softmasked.dna.gz";

        foreach my $f ($dna_file, $masked_file, $soft_file) {
          if (not -e $f or not -s $f) {
            $log->error("Could not find DNA file for $gspecies ($f)\n");
          }
        }
		
        my $target_dna_file =  "$dna_dir/".$gspecies.".$WS_name.genomic.fa.gz";
        my $target_masked = "$dna_dir/".$gspecies.".$WS_name.genomic_masked.fa.gz";
        my $target_soft = "$dna_dir/".$gspecies.".$WS_name.genomic_softmasked.fa.gz";

	$wormbase->run_command("gzip -9 -c $dna_file > $target_dna_file",$log);
	$wormbase->run_command("cp -f $soft_file $target_soft", $log);
	$wormbase->run_command("cp -f $masked_file $target_masked", $log);

        map { $copied_files{$_} = 1 } ($dna_file, $masked_file, $soft_file);
	
      } elsif ($wb->assembly_type eq 'chromosome') {
	# copy the chromosome files across
	#$wormbase->run_command("cp -f -R $chromdir/*.dna* $dna_dir/", $log);
	#$wormbase->run_command("gzip -f -9 $dna_dir/*.dna", $log);
        
	# do for each type of unmasked and masked sequence
	foreach my $type_pair (["", "genomic"],  
                               ["_masked", "genomic_masked"],
                               ["_softmasked", "genomic_softmasked"]) {
        
	  # delete any existing whole-genome files produced by running this script more than once
          my ($src_type, $tgt_type) = @$type_pair;

	  unlink "$dna_dir/${gspecies}.${WS_name}.${tgt_type}.fa";         
          $wormbase->run_command("touch $dna_dir/${gspecies}.${WS_name}.${tgt_type}.fa",$log);

	  foreach my $chrom ($wb->get_chromosome_names(-mito => 1, -prefix => 1)) {
	    my $chrom_file = "$chromdir/$chrom"; # basic form of the dna file
            
	    # is the data gzipped?
	    if (-e "$chrom_file${src_type}.dna") {
              my $source = "${chrom_file}${src_type}.dna";
	      $wormbase->run_command("cat $source >> $dna_dir/${gspecies}.${WS_name}.${tgt_type}.fa", $log);
              $copied_files{$source} = 1;
	    } elsif (-e "${chrom_file}${src_type}.dna.gz") {
              my $source = "${chrom_file}${src_type}.dna.gz";
	      $wormbase->run_command("gunzip -c $source >> $dna_dir/${gspecies}.${WS_name}.${tgt_type}.fa", $log);
              $copied_files{$source} = 1;
	    } else {
              $log->error("$gspecies : missing file: $chrom_file${src_type}.dna\n");
            }
	  }
          
	  # gzip the resulting file
	  $wormbase->run_command("gzip -9 -f $dna_dir/${gspecies}.${WS_name}.${tgt_type}.fa", $log);	
	}
      }

      # copy over outstanding dna files
      foreach my $dna_file (glob("$chromdir/*.dna.gz")) {
        if (not exists $copied_files{$dna_file}) {
          my ($prefix) = $dna_file =~ /$chromdir\/(\S+)\.dna.gz/;
          my $target = "$dna_dir/${gspecies}.${WS_name}.$prefix.fa.gz";
          $wormbase->run_command("cp -f $dna_file $target", $log);
        }
      }
        
      my @agp_files = glob("$chromdir/*.agp");
      
      if (@agp_files) {
        my $target_agp_file =  "$dna_dir/${gspecies}.${WS_name}.assembly.agp"; 
        
        if (scalar(@agp_files) == 1) {
          # just the one - assume its for whole genome and copy it across
          my ($single_file) = @agp_files;
          $wormbase->run_command("cat $single_file | gzip -9 -c > ${target_agp_file}.gz", $log); 
        } else {
          # assume per-chromosome
          unlink $target_agp_file;
          $wormbase->run_command("touch $target_agp_file", $log);
          foreach my $chrom ($wb->get_chromosome_names(-mito => 1, -prefix => 1)) {
            my $agp = "$chromdir/$chrom.agp";
            if (-e "$chromdir/$chrom.agp") {
              $wormbase->run_command("cat $agp >> $target_agp_file", $log);
            } else {
              $log->error("$gspecies : missing file: $chromdir/$chrom.agp\n");
            }
          }
          $wormbase->run_command("gzip -9 -f $target_agp_file", $log);
        }
     
      } else {
        $log->error("$gspecies : unknown assembly_type\n");
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
    my $chromdir = $wb->chromosomes;

    if (-e "$chromdir") {
      my $gff_dir = "$targetdir/species/$gspecies";
      mkpath($gff_dir,1,0775);

      #concatenated whole genome file require for all species
      my $whole_filename = "$gspecies.$WS_name.annotations.gff2"; 
      $wormbase->run_command("rm -f $gff_dir/$whole_filename", $log);

      if($wb->assembly_type eq 'contig') {
        if (-e "$chromdir/$species.gff") { # tierII does it this way
          $wormbase->run_command("cp -f -R $chromdir/$species.gff $gff_dir/$whole_filename", $log);
        } else { 
          $log->error("$chromdir/$species.gff missing\n");
        }
      } else {
        $wormbase->run_command("cat $chromdir/*.gff* > $gff_dir/$whole_filename", $log);
        #$wormbase->run_command("cp -f $chromdir/*.gff* $gff_dir/", $log); #individual files too
      }

      # add supplementary and nGASP GFF
      my $ngaspdir;
      if($species eq 'elegans') {
	my $supdir = "$chromdir/SUPPLEMENTARY_GFF";
	my @gfffiles = glob("$supdir/*.gff");
	foreach my $sup (@gfffiles){
          $wb->run_command("cat $sup >> $gff_dir/$whole_filename", $log);
	}
	$ngaspdir = $supdir;
      }
      $ngaspdir = $wb->database("$species")."/nGASP" unless $ngaspdir;
	
      # nGASP - zcat files stored under primaries (or BUILD_DATA for C_ele) on to FTP full GFF file.
      if(-e $ngaspdir){
	my @ngasp_methods = qw(augustus fgenesh jigsaw mgene);
	foreach my $method(@ngasp_methods){
          my $file = "$ngaspdir/$species.$method.gff2.gz";
          if(-e $file){
            $wb->run_command("zcat $file >> $gff_dir/$whole_filename", $log);
          } else { 
            $log->error("$file missing\n");
          }
	}
      } else { 
        $log->write_to("Warning: no ngasp for $gspecies\n");
      }
      
      $wormbase->run_command("gzip -9 -f $gff_dir/$whole_filename",$log);

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
    my $gff_dir = "$targetdir/species/$gspecies";
    mkpath($gff_dir,1,0775);

    my $source = sprintf("%s/%s.gff3", $wb->chromosomes, $species);
    my $target = sprintf("%s/%s.%s.annotations.gff3", $gff_dir, $gspecies, $WS_name);

    if (-e $source) {
      $wb->run_command("cp -f $source $target", $log);
      $wb->run_command("gzip -9 -f $target",$log);
    } else {
      $log->error("No GFF3 data found for species $species\n");
    }

  }


  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying GFF files\n\n");
}


############################################
# copy across supplementary GFF files
############################################

sub copy_supplementary_gff_files{
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying supplementary gff files\n");

  if($wormbase->species eq'elegans') {
    my $chromdir = $wormbase->chromosomes;
    if (-e "$chromdir/SUPPLEMENTARY_GFF") {
      my $sgff_dir = "$targetdir/species/c_elegans/SUPPLEMENTARY_GFF";
      mkpath($sgff_dir,1,0775);
      foreach my $gff_file (glob("$chromdir/SUPPLEMENTARY_GFF/*.gff")) {
        my ($fname) = $gff_file =~ /SUPPLEMENTARY_GFF\/(\S+)\.gff/;
        $wormbase->run_command("cat $gff_file | gzip -9 -c > $sgff_dir/c_elegans.$WS_name.$fname.gff2.gz", $log);
      }
    }
  }
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying supplementary GFF\n\n");
}


############################################
# copy across comparative analysis
############################################

sub copy_compara {

  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: copying clustal sql dump\n");

  mkpath("$targetdir/COMPARATIVE_ANALYSIS",1,0775);
  
  my $chromdir = $wormbase->autoace;
  my $clust_src = "$chromdir/wormpep_clw.sql.bz2";
  my $compara_src = "$chromdir/compara.tar.bz2";

  if (-e $clust_src) {
    my $target = "$targetdir/COMPARATIVE_ANALYSIS/wormpep_clw.${WS_name}.sql.bz2";
    $wormbase->run_command("cp -f -R $clust_src $target", $log);

  } else {
    $log->write_to("Warning: no clustal results found ($clust_src)\n");
  }
  
  if (-e $compara_src) {
    my $target = "$targetdir/COMPARATIVE_ANALYSIS/compara.$WS_name.tar.gz";
    $wormbase->run_command("cp -f -R $compara_src $target", $log);
  } else {
    $log->write_to("Warning: no clustal results found ($compara_src)\n");
  }
    
  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished copying clustalw\n\n");
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
    my $prefix = $wb->pepdir_prefix;
    my $rnadir = $wb->wormrna;

    my $ftprna_dir = "$targetdir/species/$gspecies";
    mkpath($ftprna_dir,1,0775);
	
    # Note that the following points to the most recent version of the file; 
    my $src_file = "$rnadir/${prefix}rna".$wb->get_wormbase_version.".rna"; 
    my $tgt_file = "$ftprna_dir/$gspecies.$WS_name.ncrna_transcripts.fa"; #target FTP file
    
    if (-e $src_file and -s $src_file) {
      $wormbase->run_command("gzip -9 -c $src_file > ${tgt_file}.gz",$log);
    } else {
      # this is a core file, needs to be present for all species; however,
      # not all species have ncRNA data. Therefore, create an empty placeholder
      $wormbase->run_command("touch $tgt_file && gzip $tgt_file", $log);
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
    my $gff_dir = "$targetdir/species/$gspecies";
    mkpath($gff_dir,1,0775);

    my $target = sprintf("%s/%s.%s.ncrna_transcripts.fa.gz", $gff_dir, $gspecies, $WS_name);
    my $unzipped_source = sprintf("%s/SEQUENCES/%s.ncrna.fa", $wb->autoace, $species);

    if (-e $unzipped_source) {
      $wb->run_command("gzip -9 -c $unzipped_source > $target", $log);
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

  my $obo_dir = $wormbase->primaries . "/citace/temp_unpack_dir/home/citace/Data_for_${WS_name}/Data_for_Ontology/";
  my $ace_ontology_dir = "$ace_dir/ONTOLOGY";
  my $ftp_ontology_dir = "$targetdir/ONTOLOGY";

  mkpath($ace_ontology_dir,1,0775);
  mkpath($ftp_ontology_dir,1,0775);

  $wormbase->run_command("cp -f $obo_dir/*.obo $ace_ontology_dir/", $log);
  foreach my $file (glob("$ace_ontology_dir/*.*")) {
    $wormbase->run_command("cp -f $file $ftp_ontology_dir/", $log);
  }

  # change group ownership
  $wormbase->run_command("chgrp -R  worm $ace_ontology_dir", $log);  

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
  my $srcdir = $wormbase->autoace;
  my $annotation_dir = "$targetdir/species/$gspecies/annotation";

  mkpath($annotation_dir,1,0775);

  #
  # Oligo mapping
  #
  foreach my $file (glob("$srcdir/*oligo_mapping")) {
    my ($pre) = $file =~ /$srcdir\/(\S+_oligo_mapping)/;
    my $target = "$annotation_dir/$gspecies.$WS_name.${pre}.txt.gz";
    
    $wormbase->run_command("cat $file | gzip -9 -c > $target", $log);
  }
  
  #
  # RNASeq gene expression data
  #
  my $expr = "$srcdir/expr.tar.gz";
  if (-e $expr) {
    my $target = "$annotation_dir/$gspecies.$WS_name.SRA_gene_expression.tar.gz";
    $wormbase->run_command("cp -f $expr $target", $log);
  } else {
    $log->write_to("Warning: gene expression file for $gspecies not found ($expr)\n");
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
    my $peppre = $wb->pepdir_prefix;
    my $rel_last_built = $wb->get_wormbase_version;

    my $src = $wb->wormpep;
    my $tgt = "$targetdir/species/$gspecies";
    mkpath($tgt,1,0775);

    if ($rel_last_built == $WS) {
      # tar up the latest wormpep release and copy across (files added in next loop)
      my $tgz_file = "$tgt/$gspecies.$WS_name.${peppre}pep_package.tar.gz";
      my $command = "tar -c -z -h -f $tgz_file -C $src";
      
      # for tierIIs, the *pep package often does not change between 
      # releases. Given that it includes files for differences since last
      # release etc, forward propagation of the complete package and
      # incrementing the suffices is not appropriate. 
      
      my @wormpep_files = $wb->wormpep_files;
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
      my $target_pepfile = "$tgt/$gspecies.${WS_name}.protein.fa.gz";
      open(my $source_pep_fh, $source_pepfile);
      while(<$source_pep_fh>) {
        /^\>(\S+)\s+\S+\s+(\S+)/ and do {
          $gene_ids{$1} = $2;
        };
      }
      close($source_pep_fh);

      $wb->run_command("gzip -9 -c $source_pepfile > $target_pepfile", $log);
    } else {
      $log->error("Could not find peptide file for $gspecies ($source_pepfile)\n");
    }

    if (-e $source_cdnafile) {
      my $target_cdnafile = "$tgt/$gspecies.${WS_name}.cds_transcripts.fa.gz";
      open(my $source_cdna_fh, $source_cdnafile);
      open(my $target_cdna_fh, "| gzip -9 -c > $target_cdnafile");
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
      $log->error("Could not find transcript file for $gspecies ($source_cdnafile)\n");
    }


    if ($WS == $rel_last_built) {
      my $source =$wb->wormpep . "/${peppre}pep$WS";
      my $target = "$wormpep_ftp_root/${peppre}pep$WS";

      $wormbase->run_command("cp -f $source $target", $log);
      $wormbase->run_command("chgrp worm $target", $log);

      # Uniprot require symlinks to "wormpep" and "brigpep"
      if ($species eq 'elegans' or $species eq 'briggsae') {
        my $link_target = "$wormpep_ftp_root/${peppre}pep";
        $wormbase->run_command("ln -sf $target $link_target", $log);
        $wormbase->run_command("chgrp worm $link_target", $log);
      }
    } else {
      $log->write_to("Skipping $gspecies databases cp as NOT rebuilt for WS$WS.\n");
    }


  }

  #tierIII's have no "pep" package, so lift protein and transcript files from build dir
  my %tierIII = $wormbase->tier3_species_accessors;
  foreach my $t3 (keys %tierIII){ 
    next if exists $skip_species{$t3};
    next if @only_species and not exists $only_species{$t3};

    my $wb = $tierIII{$t3};
    my $gspecies = $wb->full_name('-g_species' => 1);

    my $src =$wb->basedir . "/WORMPEP/".$wb->pepdir_prefix."pep$WS/".$wb->pepdir_prefix."pep.$WS_name.fa.gz";
    my $tgt = "$targetdir/species/$gspecies";
    mkpath($tgt,1,0775);

    my $src_pep_file = $wb->autoace . "/SEQUENCES/$t3.prot.fas";
    my $src_cdna_file = $wb->autoace . "/SEQUENCES/$t3.cdna.fas";

    my $tgt_pep_file = "$tgt/$gspecies.$WS_name.protein.fa.gz";
    my $tgt_cdna_file = "$tgt/$gspecies.$WS_name.cds_transcripts.fa.gz";

    if (-e $src_pep_file) {
      $wb->run_command("gzip -9 -c $src_pep_file > $tgt_pep_file", $log);
    } else {
      $log->error("Could not find $src_pep_file for $t3; worrysome\n");
    }

    if (-e $src_cdna_file) {
      $wb->run_command("gzip -9 -c $src_cdna_file > $tgt_cdna_file", $log);
    } else {
      $log->error("Could not find $src_cdna_file for $t3; worrysome\n");
    }
  }

  #
  # finally, need to update ftp/databases/wormpep 
  #

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

    my $tdir = "$targetdir/species/$gspecies";
    mkpath("$tdir",1,0775);

    my $target = "$targetdir/species/$gspecies/$gspecies.$WS_name.ests.fa.gz";
    my $est_dir = $wb->build_data . "/cDNA/$species";

    open(my $outfh, "| gzip -c -9 > $target") 
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
# extract confirmed genes
###################################

sub extract_confirmed_genes{

  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Extracting elegans confirmed genes\n");

  my $gspecies = $wormbase->full_name('-g_species' => 1);

  my $ace_dir = $wormbase->autoace;
  my $annotation_dir = "$targetdir/species/$gspecies/annotation";
  my $wormdna = $wormbase->wormpep."/".$wormbase->pepdir_prefix."pep.dna".$wormbase->get_wormbase_version;
  my $seqio  = Bio::SeqIO->new('-file' => "$wormdna", '-format' => 'fasta');

  my $db = Ace->connect(-path => "$ace_dir/") || die print "ACE connection failure: ", Ace->error;
  my $query = "Find elegans_CDS; Confirmed";

  my @confirmed_genes   = $db->fetch(-query=>$query);
  my %conf_genes;
  map { $conf_genes{$_->name} = 1 } @confirmed_genes;

  mkpath("$annotation_dir",1,0775);

  my $out_file = "$annotation_dir/$gspecies.$WS_name.confirmed_genes.fa";

  open(OUT,">$out_file") or 
      $log->write_to("Couldn't write to $out_file\n");

  while(my $seq = $seqio->next_seq){
    if($conf_genes{$seq->id}) {
      print OUT "\n>".$seq->id."\n";
      print OUT $wormbase->format_sequence($seq->seq);
    }
  }
  close(OUT);

  $wormbase->run_command("gzip -9 -f $out_file", $log);

  $db->close;

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished extracting confirmed genes\n\n");

  return(0);
}

sub extract_ko {
  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Extracting Knockout Consortium Data\n");

  my $gspecies = $wormbase->full_name('-g_species' => 1);
  my $annotation_dir = "$targetdir/species/$gspecies/annotation";
  mkpath("$annotation_dir",1,0775);

  my $outfile = "$annotation_dir/${gspecies}.${WS_name}.knockout_consortium_alleles.xml.bz2";

  $wormbase->run_script("dump_ko.pl -file $outfile",$log);
  $wormbase->run_command("chgrp -R  worm $annotation_dir", $log);  
  
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
  my $annotation_dir = "$targetdir/species/$gspecies/annotation";
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

  
  my $out = "$annotation_dir/$gspecies.$WS_name.cdna2orf.txt";
  open(OUT, ">$out") or $log->write_to("Couldn't open $out\n");

  foreach my $key (keys %cDNA2orf){
    print OUT "$key,$cDNA2orf{$key}\n";
  }
  close(OUT);

  $wormbase->run_command("gzip -9 -f $out", $log);

  $runtime = $wormbase->runtime;
  $log->write_to("$runtime: Finished making cdna2orf list\n\n");  
  
}

################################################################################
# make list of WBGene IDs to CGC name and Sequence name
################################################################################

sub make_geneID_list {

  my $runtime = $wormbase->runtime;
  $log->write_to("$runtime: making Gene ID list\n");
  # For each 'live' Gene object, extract 'CGC_name' and 'Sequence_name' fields (if present)

  my $gspecies = $wormbase->full_name('-g_species' => 1);
  my $annotation_dir = "$targetdir/species/$gspecies/annotation";
  my $ace_dir = $wormbase->autoace;

  my $tace    = $wormbase->tace;
  my $command = "Table-maker -p $ace_dir/wquery/gene2cgc_name_and_sequence_name.def\nquit\n";
  my $dir     = "$ace_dir";
  my $out     = "$annotation_dir/$gspecies.$WS_name.geneIDs.txt";
  mkpath("$annotation_dir",1,0775);

  open (OUT, ">$out") || croak "Couldn't open $out\n";
  open (TACE, "echo '$command' | $tace $dir | ") or $log->write_to("cant access $dir\n");  
  while (<TACE>){
      if (m/^\"/){
	  s/\"//g;
	  tr/\t/,/;
	  print OUT "$_";
      }
  }
  close(TACE);
  close(OUT);

  $wormbase->run_command("gzip -9 -f $out", $log);

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
  my $annotation_dir = "$targetdir/species/$gspecies/annotation";
  my $ace_dir = $wormbase->autoace;

  my $tace    = $wormbase->tace;
  my $command = "Table-maker -p $ace_dir/wquery/pcr_product2gene.def\nquit\n";
  my $dir     = "$ace_dir";
  my $out     = "$annotation_dir/$gspecies.$WS_name.pcr_product2gene.txt";
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

  $wormbase->run_command("gzip -9 -f $out", $log);

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

    my $source_file = "$blast_dir/${species}_best_blastp_hits";
    # this script might be run more than once if there are problems
    if(-e $source_file || -e "$source_file.gz") { 
      my $protein_dir = "$targetdir/species/$gspecies";
      mkpath($protein_dir,1,0775);

      my $target_file = "$protein_dir/$gspecies.$WS_name.best_blastp_hits.txt.gz";
      # this script might be run more than once if there are problems
      $wormbase->run_command("gzip -9 -f $source_file",$log) if (! -e "$source_file.gz"); 
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
    my $tgt_dir = "$targetdir/species/$gspecies";

    my $in_filename = "$tgt_dir/$gspecies.$WS_name.annotations.gff2.gz"; 

    if (not -e $in_filename) {
      $log->error("Could not make GBrowse-ready GFF for $species; in file $in_filename not found\n");
      next;
    }

    my $out_filename = "$tgt_dir/$gspecies.$WS_name.GBrowse.gff2.gz";
    my $command = "perl $ENV{CVS_DIR}/web_data/process_elegans_gff-standalone.pl -species $species -database ". $wormbase->autoace;
    if ($debug) {
      $command .= " -debug $debug ";
    }

    my @pipe = ("gunzip -c $in_filename",
                
                "gzip -c -9 > $out_filename");

    $wb->run_command( "gunzip -c $in_filename | $command | gzip -c -9 > $out_filename", $log);
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
sub make_frozen_links {
  my $lnpath = $wormbase->ftp_site;
  if ($WS =~ /\d+5$/ || $WS =~ /\d+0$/) {
    $wormbase->run_command("ln -sf $lnpath/releases/$WS_name $lnpath/releases/FROZEN_RELEASES/$WS_name", $log);
  }
}

##################################################################################

sub CheckSize {
  my ($first,$second)=@_;
  my $F_SIZE = (stat("$first"))[7];
  my $S_SIZE = (stat("$second"))[7];
  if ($F_SIZE != $S_SIZE) {
    $log->error("\tERROR: $first SRC: $F_SIZE TGT: $S_SIZE - not same size, please check\n");
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
      foreach my $tsp (split(/,/, $sp)) {
        if ($tsp eq 'TIER3') {
          foreach my $ntsp (keys %t3_species) {
            $stanzas[-1]->{species}->{$ntsp} = 1;
          }
        } elsif ($tsp eq 'TIER2') {
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
        $file =~ s/WSREL/$WS_name/g;
        $file = "$parent/$file";
        $count += &checkfile($file);
      }
    } else {
      foreach my $worm (keys %species) {
        next if exists $skip_species{$worm};
        next if @only_species and not exists $only_species{$worm};
        
        next if (keys %{$block->{species}} and not exists $block->{species}->{$worm});
        
        my $gspecies = $species{$worm}->full_name('-g_species'=>1);
        my $gWORM = $species{$worm}->pepdir_prefix;
        my $species = $species{$worm}->species;
        
        my $parent = $block->{parent};
        $parent =~ s/GSPECIES/$gspecies/g;
        $parent =~ s/WSREL/$WS_name/g;
        $parent = "$targetdir/$parent";
        
        #print "PARENT = $parent\n";
        
        foreach my $ofile (@{$block->{files}}) {
          my $file = $ofile;
          
          $file =~ s/GSPECIES/$gspecies/g;
          $file =~ s/WSREL/$WS_name/g;
          $file =~ s/WORM/$gWORM/g;
          $file = "$parent/$file";
          #print "FILE = $file\n";
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
[elegans]species/c_elegans/annotation
c_elegans.WSREL.affy_oligo_mapping.txt.gz
c_elegans.WSREL.agil_oligo_mapping.txt.gz
c_elegans.WSREL.pcr_product2gene.txt.gz
c_elegans.WSREL.geneIDs.txt.gz
c_elegans.WSREL.gsc_oligo_mapping.txt.gz
c_elegans.WSREL.cdna2orf.txt.gz
c_elegans.WSREL.confirmed_genes.fa.gz

[elegans]species/c_elegans
GSPECIES.WSREL.WORMpep_package.tar.gz

[elegans,briggsae]species/GSPECIES
GSPECIES.WSREL.assembly.agp.gz

# tierII specific stuff
[TIER2]species/GSPECIES
GSPECIES.WSREL.best_blastp_hits.txt.gz
GSPECIES.WSREL.intergenic_sequences.fa.gz
GSPECIES.WSREL.GBrowse.gff2.gz
GSPECIES.WSREL.ncrna_transcripts.fa.gz
GSPECIES.WSREL.annotations.gff2.gz

[TIER3]species/GSPECIES
GSPECIES.WSREL.annotations.gff3.gz

# for all species
[ALL]species/GSPECIES
GSPECIES.WSREL.genomic.fa.gz
GSPECIES.WSREL.genomic_masked.fa.gz
GSPECIES.WSREL.genomic_softmasked.fa.gz
GSPECIES.WSREL.protein.fa.gz
GSPECIES.WSREL.cds_transcripts.fa.gz


[]acedb
files_in_tar
md5sum.WSREL
models.wrm.WSREL

[]ONTOLOGY
anatomy_association.WSREL.wb
anatomy_ontology.WSREL.obo
gene_association.WSREL.wb.cb
gene_association.WSREL.wb.ce
gene_association.WSREL.wb
gene_ontology.WSREL.obo
phenotype_association.WSREL.wb
phenotype_ontology.WSREL.obo


[]COMPARATIVE_ANALYSIS
compara.WSREL.tar.gz
wormpep_clw.WSREL.sql.bz2
