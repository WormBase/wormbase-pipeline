
#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $WORMBASE_GENOMES = {
  'brugia_malayi_prjna10729'             => 1,
  'onchocerca_volvulus_prjeb513'         => 1,
  'strongyloides_ratti_prjeb125'         => 1,
  'pristionchus_pacificus_prjna12644'    => 1,
  'caenorhabditis_angaria_prjna51225'    => 1,
  'caenorhabditis_elegans_prjna13758'    => 1,
  'caenorhabditis_briggsae_prjna10731'   => 1,
  'caenorhabditis_brenneri_prjna20035'   => 1,
  'caenorhabditis_remanei_prjna53967'    => 1,
  'caenorhabditis_japonica_prjna12591'   => 1,
  'caenorhabditis_sinica_prjna194557'    => 1,  
  'caenorhabditis_tropicalis_prjna53597' => 1,
  'panagrellus_redivivus_prjna186477'    => 1,
}; 

my (
  $wb_rel_num,
  $rel_num,
  $checksum, 
  $verbose,
  $staging_tl_dir,
  $production_tl_dir,
  $staging,
  $production,
  $copy_previous,
    );

$rel_num = "666";
$wb_rel_num = "666";

&GetOptions(
  'relnum=s'         => \$rel_num,
  'wbrelnum=s'       => \$wb_rel_num,
  'checksum'         => \$checksum,
  'verbose'          => \$verbose,
  'stagingdir=s'     => \$staging_tl_dir,
  'releasedir=s'     => \$production_tl_dir,
  'staging'          => \$staging,
  'release'          => \$production,
    );

my $release = "WBPS${rel_num}";
my $prev_release = "WBPS" . ($rel_num - 1);
my $wb_release = "WS${wb_rel_num}";

my $ftp_root = "/ebi/ftp/pub/databases/wormbase";

$staging_tl_dir            = "$ftp_root/staging/parasite/releases" if not defined $staging_tl_dir;
$production_tl_dir         = "$ftp_root/parasite/releases" if not defined $production_tl_dir;
my $wb_ftp_root_staging    = "../../../../../../../releases";
my $wb_ftp_root_production = "../../../../../../releases";

if ($staging) {
  &make_core_symlinks($wb_ftp_root_staging, $staging_tl_dir );
  &make_md5sums($staging_tl_dir);
} 

if ($production) {
  if (-d "$staging_tl_dir/$release" and not -d "$production_tl_dir/$release") {
    system("mv $staging_tl_dir/$release $production_tl_dir") and die "Could not mv $staging_tl_dir/$release to $production_tl_dir\n";
  }
  if (-d "$production_tl_dir/$release") {
    &make_core_symlinks($wb_ftp_root_production, $production_tl_dir );
    &make_md5sums($production_tl_dir);
  } 
}


#####################
sub make_md5sums {
  my ($tl_dir) = @_;

  my $targetdir = "$tl_dir/$release";
  my $checksum_file = "CHECKSUMS";

  my @files;
  open(FIND, "find $targetdir/species -name '*.*' | sort |");
  while(<FIND>) {
    chomp;
    s/^$targetdir\///;
    $verbose and print STDERR "Will checksum: $_\n";
    push @files, $_;
  }

  system("cd $targetdir && md5sum @files > $checksum_file") and die "Could not calc checksums\n";

}


#############################################
sub make_core_symlinks {
  my ($core_ftp_root, $this_tl_dir) = @_;

  foreach my $wb_genome (keys %$WORMBASE_GENOMES) {

    my ($genus_pre, $genus_suf, $spe, $bioproject) = $wb_genome =~ /(\S)(\S+)_(\S+)_(\S+)/; 
    $bioproject = uc($bioproject);
    
    my $ps_species_name = "${genus_pre}${genus_suf}_${spe}";
    my $wb_species_name = "${genus_pre}_${spe}";

    my $link_dir_dest = join("/", $this_tl_dir, $release, "species", $ps_species_name, $bioproject);
    my $link_dir_source = "$core_ftp_root/$wb_release/species/$wb_species_name/$bioproject";
    mkpath $link_dir_dest if not -d $link_dir_dest;

    foreach my $fsuffix( "genomic.fa.gz", 
                         "genomic_softmasked.fa.gz", 
                         "genomic_masked.fa.gz",
                         "CDS_transcripts.fa.gz", 
                         "mRNA_transcripts.fa.gz", 
                         "protein.fa.gz", 
                         "ncRNA_transcripts.fa.gz",
                         "annotations.gff3.gz",
                         "canonical_geneset.gtf.gz") {
      my $link_fname_dest = join(".", $ps_species_name, $bioproject, $release, $fsuffix);
      my $link_fname_source = join(".", $wb_species_name, $bioproject, $wb_release, $fsuffix);

      unlink "$link_dir_dest/$link_fname_dest" 
          if -e "$link_dir_dest/$link_fname_dest" or -l "$link_dir_dest/$link_fname_dest";

      system("cd $link_dir_dest && ln -s $link_dir_source/$link_fname_source $link_fname_dest") 
          and die "Could not create symlink to $link_dir_source/$link_fname_source in $link_dir_dest\n";
    }
  }
}
