
#!/usr/bin/env perl

use strict;
use Getopt::Long;
use DBI;
use FindBin;
use File::Path;

my $DUMP_GENOME_SCRIPT = "dump_genome.pl";
my $DUMP_TRANSCRIPTS_SCRIPT = "dump_transcripts.pl";
my $DUMP_GFF3_SCRIPT = "dump_gff3.pl";

my $WORMBASE_CORE = {
  'brugia_malayi_prjna10729'           => 1,
  'onchocerca_volvulus_prjeb513'       => 1,
  'strongyloides_ratti_prjeb125'       => 1,
  'pristionchus_pacificus_prjna12644'  => 1,
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
  $all,
  $tl_out_dir,
  $wb_rel_num,
  $rel_num,
  $core_symlinks,
  $checksum, 
  $verbose,
  $prev_rel_tl_dir,
    );

$user = 'ensro';
$tl_out_dir = ".";
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
  'all'              => \$all,
  'outdir=s'         => \$tl_out_dir,
  'relnum=s'         => \$rel_num,
  'wbrelnum=s'       => \$wb_rel_num,
  'checksum'         => \$checksum,
  'verbose'          => \$verbose,
  'coresymlinks'     => \$core_symlinks,
  'prevreldir=s'     => \$prev_rel_tl_dir,  # copy data from previous release dir (if it exists)
    );

my $release = "WBPS${rel_num}";
my $prev_release = "WBPS" . ($rel_num - 1);
my $wb_release = "WS${wb_rel_num}";


if ($g_nomask or $g_smask or $g_hmask or $cds_tran or $mrna_tran or $prot or $gff3 or $all) {

  my @genome_db_pairs = &get_databases($rel_num, @genomes);

  foreach my $pair (@genome_db_pairs) {
    my ($genome, $dbname) = @$pair;

    my ($species, $bioproject) = $genome =~ /^(\S+_\S+)_(\S+)/;
    $bioproject = uc($bioproject);
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
                "$DUMP_GENOME_SCRIPT -dbname $dbname",
                $prevdir) if $g_nomask or $all;

    &write_file($species, 
                $bioproject,
                $outdir, 
                "genomic_softmasked.fa", 
                "$DUMP_GENOME_SCRIPT -dbname $dbname -softmask",
                $prevdir) if $g_smask or $all;

    &write_file($species, 
                $bioproject,
                $outdir, 
                "genomic_masked.fa", 
                "$DUMP_GENOME_SCRIPT -dbname $dbname -mask",
                $prevdir) if $g_hmask or $all;

    &write_file($species, 
                $bioproject,
                $outdir, 
                "CDS_transcripts.fa", 
                "$DUMP_TRANSCRIPTS_SCRIPT -dbname $dbname -cds",
                $prevdir) if $cds_tran or $all;

    &write_file($species, 
                $bioproject,
                $outdir, 
                "mRNA_transcripts.fa", 
                "$DUMP_TRANSCRIPTS_SCRIPT -dbname $dbname -mrna",
                $prevdir) if $mrna_tran or $all;

    &write_file($species, 
                $bioproject,
                $outdir, 
                "protein.fa", 
                "$DUMP_TRANSCRIPTS_SCRIPT -dbname $dbname -pep",
                $prevdir) if $prot or $all;

    &write_file($species, 
                $bioproject,
                $outdir, 
                "annotations.gff3", 
                "$DUMP_GFF3_SCRIPT -dbname $dbname",
                $prevdir) if $gff3 or $all;
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
    system("cp $prevfile $outdir/$fname") and die "Could not copy $prevfile to $outdir/$fname\n";
  } else {
    $verbose and print STDERR "    Dumping file using $script\n";

    my $this_cmd = "$FindBin::Bin/$script -host $host -port $port -user $user -outfile $outdir/$fname";

    system($this_cmd) and die "Could not successfully run $this_cmd\n";
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
      next if exists $WORMBASE_CORE->{$genome};

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

  foreach my $core_genome (keys %$WORMBASE_CORE) {
    next if @genomes and not grep { $core_genome eq $_ } @genomes;

    my ($genus_pre, $genus_suf, $spe, $bioproject) = $core_genome =~ /(\S)(\S+)_(\S+)_(\S+)/; 
    $bioproject = uc($bioproject);
    
    my $ps_species_name = "${genus_pre}${genus_suf}_${spe}";
    my $wb_species_name = "${genus_pre}_${spe}";

    my $link_dir_dest = join("/", $tl_out_dir, $release, "species", $ps_species_name, $bioproject);
    my $link_dir_source = "../../../../../../releases/$wb_release/species/$wb_species_name/$bioproject";

    mkpath $link_dir_dest if not -d $link_dir_dest;

    foreach my $fsuffix( "genomic.fa.gz", 
                         "genomic_softmasked.fa.gz", 
                         "genomic_masked.fa.gz",
                         "CDS_transcripts.fa.gz", 
                         "mRNA_transcripts.fa.gz", 
                         "protein.fa.gz", 
                         "ncRNA_transcripts.fa.gz",
                         "annotations.gff3.gz") {
      my $link_fname_dest = join(".", $ps_species_name, $bioproject, $release, $fsuffix);
      my $link_fname_source = join(".", $wb_species_name, $bioproject, $wb_release, $fsuffix);

      system("cd $link_dir_dest && ln -s $link_dir_source/$link_fname_source $link_fname_dest") 
          and die "Could not create symlink for $core_genome $fsuffix\n";
    }
  }
}
