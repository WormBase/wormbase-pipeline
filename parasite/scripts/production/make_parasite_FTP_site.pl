#!/usr/bin/env perl
# Create target FTP directory
# If possible link to entries in wormbase FTP
# Otherwise copy from source
# Then, recreate the checksums
use strict;
use Getopt::Long;
use File::Path;
use File::Basename;
use File::Spec;
my (
  $wormbase_release_ftp_dir,
  $source_dir,
  $wbps_release_ftp_dir,
  $wbps_version
   );
&GetOptions(
 'wormbase_release_ftp_dir=s' => \$wormbase_release_ftp_dir,
 'source_dir=s' => \$source_dir,
 'wbps_release_ftp_dir=s' => \$wbps_release_ftp_dir,
 'wbps_version=i' => \$wbps_version,
) ;
my $usage = " Usage: $0 --wbps_version=\$PARASITE_VERSION --source_dir=<where folders with individual species are> --wormbase_release_ftp_dir=<release tied to this WBPS version> --wbps_release_ftp_dir=<target directory>";
die ("--source_dir not a directory: $source_dir . $usage") unless -d $source_dir;
die ("--wormbase_release_ftp_dir not a directory: $wormbase_release_ftp_dir . $usage") unless -d $wormbase_release_ftp_dir;
die ($usage) unless $wbps_release_ftp_dir;
die ($usage) unless $wbps_version;

for my $path_species (glob "$source_dir/*") {
  my $species = basename $path_species;
  my ($spe, $cies) = split(/_/, $species);
  for my $this_source_dir ( glob "$source_dir/$species/*" ) {
     my $bioproject = basename $this_source_dir;
     my $this_target_dir = "$wbps_release_ftp_dir/species/$species/$bioproject";
     mkpath $this_target_dir if not -d $this_target_dir;
     my $putative_wormbase_dir = join("/", $wormbase_release_ftp_dir,"species", lc((substr $spe, 0, 1 ) . "_" . $cies) , uc($bioproject));
     if ( -d $putative_wormbase_dir ) {
        print localtime ." ". $species . " making symlinks $putative_wormbase_dir -> $this_target_dir \n";
        &make_symlinks_to_wormbase_species (
          "$species.$bioproject.WBPS$wbps_version",
          $putative_wormbase_dir,
          $this_target_dir
        );
     } else {
        my $cp_cmd = "rsync -a --include='*.gz' --exclude '*' $this_source_dir/ $this_target_dir/";
        print localtime . " $species $cp_cmd\n";
        system($cp_cmd) and die("Failed: $cp_cmd");
     }
  }
}
print localtime . " Finished moving species files, remaking checksums file \n" ;
my $checksum_file = "CHECKSUMS";

my @files;
open(FIND, "find $wbps_release_ftp_dir/species -name '*.*' | sort |");
while(<FIND>) {
  chomp;
    s/^$wbps_release_ftp_dir\///;
    push @files, $_;
  }

system("cd $wbps_release_ftp_dir && md5sum @files > $checksum_file") and die "Could not calc checksums\n";
print localtime . " Completed \n";
#####################

sub make_symlinks_to_wormbase_species {
  my ($name_root, $wormbase_source_dir, $target_dir) = @_;
  for (glob("$target_dir/*")){
    unlink;
  }
  for my $file_type (
qw/protein.fa.gz
mRNA_transcripts.fa.gz
genomic_softmasked.fa.gz
genomic_masked.fa.gz
genomic.fa.gz
CDS_transcripts.fa.gz
canonical_geneset.gtf.gz
annotations.gff3.gz/) {
  my ($wormbase_link_target, @others) = glob("$wormbase_source_dir/*.$file_type");
  die ("No unambiguous file of type $file_type found in $wormbase_source_dir - WormBase, what happened?!") 
      unless $wormbase_link_target and not @others;
  my $l = "$target_dir/$name_root.$file_type";
  # Need relative links - the FTP site is mirrored over from the NFS directory we operate in
  my $wormbase_relative_link_target =File::Spec->abs2rel ($wormbase_link_target, $target_dir); 
  symlink ($wormbase_relative_link_target , $l) or die "Could not create symlink: $l -> $wormbase_relative_link_target ($wormbase_link_target )";
  }
}
