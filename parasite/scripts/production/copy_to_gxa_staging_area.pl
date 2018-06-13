#!/usr/bin/env perl
#
# Note: this script should only be run after the data has been put in the FTP staging area, and the correct symlinks are in place
#

use strict;
use Getopt::Long;
use File::Path qw (make_path);
use Bio::EnsEMBL::Registry;

my ($gxa_staging_dir,
    $ftp_source_dir,
    $rel_num,
    $reg_conf,
    @entries,
    $force_copy
    );

&GetOptions("gxastaging=s"    => \$gxa_staging_dir,
            "ftpsourcedir=s"  => \$ftp_source_dir,
            "reg_conf=s"      => \$reg_conf,
            "relnum=s"        => \$rel_num,
            "force_copy=s"    => \$force_copy 
    );

die "You must supply a ParaSite release number with -relnum\n" if not defined $rel_num;
die "You must supply a Registry file with -reg_conf\n" if not defined $reg_conf;

Bio::EnsEMBL::Registry->load_all( $reg_conf );

$gxa_staging_dir = "/nfs/ftp/pub/databases/wormbase/collaboration/EBI/GxA" if not defined $gxa_staging_dir;
$ftp_source_dir = "/nfs/ftp/pub/databases/wormbase/staging/parasite/releases/WBPS${rel_num}/species" if not defined $ftp_source_dir;


my $all_dbas = Bio::EnsEMBL::Registry->get_all_DBAdaptors(-GROUP => 'core');

foreach my $dba (@$all_dbas) {
  
  my $mc = $dba->get_MetaContainer();
  
  my $dbname   = $mc->dbc->dbname();
  my ($species, $db_rel_num) = $dbname =~ /^([^_+]+_[^_]+)_[^_]+_core_(\d+)/;
  next if not defined $db_rel_num;
  next if $rel_num ne $db_rel_num;

  next if $mc->get_division() ne "EnsemblParasite";
  

  my $assembly   = $mc->single_value_by_key('assembly.name');
  my $bioproject = $mc->single_value_by_key('species.ftp_genome_id');
  my $taxon      = $mc->single_value_by_key('species.taxonomy_id');
  
  my $prefix = join(".", $species, $bioproject, "WBPS${rel_num}");

  push @entries, {
    species    => $species,
    bioproject => $bioproject,
    taxon      => $taxon,
    assembly   => $assembly,
    to_copy    => [
      $prefix . ".canonical_geneset.gtf.gz",
      $prefix . ".genomic.fa.gz",
      $prefix . ".mRNA_transcripts.fa.gz",
    ]
  };
}

#
# copy the Genome and GTF files
# 
my $outdir = "$gxa_staging_dir/WBPS${rel_num}";
make_path $outdir;

foreach my $entry (@entries) {
  foreach my $file_name (@{$entry->{to_copy}}) {
    my $f = join ("/", $ftp_source_dir, $entry->{species}, $entry->{bioproject}, $file_name);
    die "Missing: $f" unless -f $f;
    next if -f "$outdir/$file_name" and not $force_copy;
    system("cp -v $f $outdir/$file_name") and die "Could not copy $f to $outdir\n";
  }
}

#
# Write the summary file
#
open(my $out_fh, ">$outdir/assembly_names.txt") or die "Could not open the assembly_names.txt file for writing\n";
foreach my $entry (sort { $a->{species} cmp $b->{species} or $a->{bioproject} cmp $b->{bioproject} } @entries) {
  printf $out_fh "%s\t%s\t%d\t%s\n", $entry->{species}, $entry->{bioproject}, $entry->{taxon}, $entry->{assembly};
}
close($out_fh);

#
# And finally, update the "latest" symlink
#
system("cd $gxa_staging_dir && rm latest && ln -sv WBPS${rel_num} latest") and die "Could not update latest symlink\n";
