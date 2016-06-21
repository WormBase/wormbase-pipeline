
#!/usr/bin/env perl

use strict;
use Getopt::Long;
use DBI;
use FindBin;
use File::Path;

my $SCRIPT_LOC = "$FindBin::Bin/../../../scripts";

my $DUMP_GENOME_SCRIPT      = "ENSEMBL/scripts/dump_genome.pl";
my $DUMP_TRANSCRIPTS_SCRIPT = "ENSEMBL/scripts/dump_transcripts.pl";
my $DUMP_GFF3_SCRIPT        = "ENSEMBL/scripts/dump_gff3.pl";
my $DUMP_GTF_SCRIPT         = "ENSEMBL/scripts/dump_gtf_from_ensembl.pl";


foreach my $script ($DUMP_GENOME_SCRIPT, $DUMP_TRANSCRIPTS_SCRIPT, $DUMP_GFF3_SCRIPT) {
  die "Could not find $script in $SCRIPT_LOC\n" if not -e "$SCRIPT_LOC/$script";
}


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
  @genomes,
  $host,
  $user,
  $port,
  $g_nomask, 
  $g_smask,
  $g_hmask,
  $prot,
  $cds_tran,
  $mrna_tran,
  $gff3,
  $gtf,
  $all,
  $tl_out_dir,
  $wb_rel_num,
  $wb_ftp_dir,
  $rel_num,
  $core_symlinks,
  $checksum, 
  $verbose,
  $prev_rel_tl_dir,
    );

$user = 'ensro';
$rel_num = "666";
$wb_rel_num = "666";

&GetOptions(
  'user=s'           => \$user,
  'host=s'           => \$host,
  'port=s'           => \$port,
  'genome=s'         => \@genomes,
  'unmaskedgenome'   => \$g_nomask,
  'softmaskedgenome' => \$g_smask,
  'hardmaskedgenome' => \$g_hmask,
  'proteins'         => \$prot,
  'cdstranscripts'   => \$cds_tran,
  'mrnatranscripts'  => \$mrna_tran,
  'gff3'             => \$gff3,
  'gtf'              => \$gtf,
  'all'              => \$all,
  'outdir=s'         => \$tl_out_dir,
  'relnum=s'         => \$rel_num,
  'wbrelnum=s'       => \$wb_rel_num,
  'wbftpdir=s'       => \$wb_ftp_dir,
  'checksum'         => \$checksum,
  'verbose'          => \$verbose,
  'coresymlinks'     => \$core_symlinks,
  'prevreldir=s'     => \$prev_rel_tl_dir,  # copy data from previous release dir (if it exists)
    );

my $release = "WBPS${rel_num}";
my $prev_release = "WBPS" . ($rel_num - 1);
my $wb_release = "WS${wb_rel_num}";

my $ftp_root = "/nfs/ftp/pub/databases/wormbase";

$tl_out_dir      = "$ftp_root/staging/parasite/releases" if not defined $tl_out_dir;
$prev_rel_tl_dir = "$ftp_root/parasite/releases" if not defined $prev_rel_tl_dir;
$wb_ftp_dir      = "$ftp_root/releases" if not defined $wb_ftp_dir;


if ($g_nomask or $g_smask or $g_hmask or $cds_tran or $mrna_tran or $prot or $gff3 or $gtf or $all) {

  my @genome_db_pairs = &get_databases($rel_num, @genomes);

  foreach my $pair (@genome_db_pairs) {
    my ($genome, $dbname) = @$pair;

    my ($species, $bioproject) = $genome =~ /^(\S+_\S+)_(\S+)/;
    $bioproject = uc($bioproject);

    if ($bioproject !~ /^PRJ/) {
      my ($pre, $suf) = $bioproject =~ /^(\S+)(PRJ\w+\d+)/;
      $bioproject = "${pre}_${suf}";
    }

    my $outdir = join("/", $tl_out_dir, $release, "species", $species, $bioproject);
    my $prevdir;

    if ($prev_rel_tl_dir) {
      $prevdir = join("/", $prev_rel_tl_dir, $prev_release, "species", $species, $bioproject);
      if (not -d $prevdir) {
        warn("Previous release data for $prevdir  not found. Will generate.\n");
        undef $prevdir;
      }
    }


    &write_file($species, 
                $bioproject,
                $outdir, 
                "genomic.fa", 
                "$DUMP_GENOME_SCRIPT  -host $host -port $port -user $user -dbname $dbname",
                $prevdir) if $g_nomask or $all;
    
    &write_file($species, 
                $bioproject,
                $outdir, 
                "genomic_softmasked.fa", 
                "$DUMP_GENOME_SCRIPT  -host $host -port $port -user $user -dbname $dbname -softmask",
                $prevdir) if $g_smask or $all;
    
    &write_file($species, 
                $bioproject,
                $outdir, 
                "genomic_masked.fa", 
                "$DUMP_GENOME_SCRIPT  -host $host -port $port -user $user -dbname $dbname -mask",
                $prevdir) if $g_hmask or $all;
    
    &write_file($species, 
                $bioproject,
                $outdir, 
                "CDS_transcripts.fa", 
                "$DUMP_TRANSCRIPTS_SCRIPT  -host $host -port $port -user $user -dbname $dbname -cds",
                $prevdir) if $cds_tran or $all;
    
    &write_file($species, 
                $bioproject,
                $outdir, 
                "mRNA_transcripts.fa", 
                "$DUMP_TRANSCRIPTS_SCRIPT  -host $host -port $port -user $user -dbname $dbname -mrna",
                $prevdir) if $mrna_tran or $all;
    
    &write_file($species, 
                $bioproject,
                $outdir, 
                "protein.fa", 
                "$DUMP_TRANSCRIPTS_SCRIPT  -host $host -port $port -user $user -dbname $dbname -pep",
                $prevdir) if $prot or $all;
    
    &write_file($species, 
                $bioproject,
                $outdir, 
                "annotations.gff3", 
                "$DUMP_GFF3_SCRIPT  -host $host -port $port -user $user -dbname $dbname",
                $prevdir) if $gff3 or $all;
    
    &write_file($species,
                $bioproject,
                $outdir,
                "canonical_geneset.gtf",
                "$DUMP_GTF_SCRIPT -host $host -port $port -user $user -dbname $dbname", 
                $prevdir) if $gtf or $all;
  }  
}

&make_md5sums() if $checksum;
&make_core_symlinks(@genomes) if $core_symlinks;

###############################

###############################
sub write_file {
  my ($species, $bioproject, $outdir, $suffix, $script, $prevdir) = @_;

  $verbose and print STDERR "  Creating $suffix file for $species/$bioproject\n";

  mkpath $outdir if not -d $outdir;
  my $fname = join(".", $species, $bioproject, $release, $suffix);
  
  if (defined $prevdir) {
    my $prevfile = join("/", $prevdir, join(".", $species, $bioproject, $prev_release, $suffix, "gz"));
    die "Could not find $prevfile\n" if not -e $prevfile;
    $fname .= ".gz";
    $verbose and print STDERR "    Copying file from previous release ($prevfile)\n";
    unlink "$outdir/$fname" if -e "$outdir/$fname";
    system("cp $prevfile $outdir/$fname") and die "Could not copy $prevfile to $outdir/$fname\n";
  } else {
    $verbose and print STDERR "    Dumping file using $script\n";

    my $this_cmd = "perl $SCRIPT_LOC/$script -outfile $outdir/$fname";

    unlink "$outdir/$fname" if -e "$outdir/$fname";
    system($this_cmd) and die "Could not successfully run $this_cmd\n";
    unlink "$outdir/${fname}.gz" if -e "$outdir/${fname}.gz";
    system("gzip -9 -n $outdir/$fname") and die "Could not successfully gzip $outdir/$fname\n";
  }
}



#####################
sub get_databases {
  my ($ps_rel_number, @genomes) = @_;

  my @dbs;


  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%d', $host, $port );
  my $dbh = DBI->connect( $dsn, $user, "",
                          { 'PrintError' => 1, 'RaiseError' => 1 } );

  
  my $sth = $dbh->prepare('SHOW DATABASES');
  $sth->execute;
  while(my ($db_name) = $sth->fetchrow_array) {
    if ($db_name =~ /^(\S+_\S+_\S+)_core_${ps_rel_number}_/) {

      my $genome = $1;
      next if @genomes and not grep { $genome eq $_ } @genomes;
      next if exists $WORMBASE_GENOMES->{$genome};

      push @dbs, [$genome, $db_name];
    }
  }

  return @dbs;
}


#####################
sub make_md5sums {

  my $targetdir = "$tl_out_dir/$release";
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
  my (@genomes) = @_;

  foreach my $wb_genome (keys %$WORMBASE_GENOMES) {
    next if @genomes and not grep { $wb_genome eq $_ } @genomes;

    my ($genus_pre, $genus_suf, $spe, $bioproject) = $wb_genome =~ /(\S)(\S+)_(\S+)_(\S+)/; 
    $bioproject = uc($bioproject);
    
    my $ps_species_name = "${genus_pre}${genus_suf}_${spe}";
    my $wb_species_name = "${genus_pre}_${spe}";

    my $link_dir_dest = join("/", $tl_out_dir, $release, "species", $ps_species_name, $bioproject);
    my $link_dir_source = "$wb_ftp_dir/$wb_release/species/$wb_species_name/$bioproject";

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

      unlink "$link_dir_dest/$link_fname_dest" if -l "$link_dir_dest/$link_fname_dest";

      system("cd $link_dir_dest && ln -s $link_dir_source/$link_fname_source $link_fname_dest") 
          and die "Could not create symlink to $link_dir_source/$link_fname_source in $link_dir_dest\n";
    }
  }
}
