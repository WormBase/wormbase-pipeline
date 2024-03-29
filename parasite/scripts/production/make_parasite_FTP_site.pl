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
  $wbps_version,
  $sync_files_skip,
  $checksums_skip,
  $repeats_dir,
   );
&GetOptions(
 'wormbase_release_ftp_dir=s' => \$wormbase_release_ftp_dir,
 'source_dir=s' => \$source_dir,
 'wbps_release_ftp_dir=s' => \$wbps_release_ftp_dir,
 'wbps_version=i' => \$wbps_version,
 'sync_files_skip' => \$sync_files_skip,
 'checksums_skip' => \$checksums_skip,
 'repeats_dir=s' => \$repeats_dir,
) ;
my $usage = " Usage: $0 --wbps_version=\$PARASITE_VERSION --source_dir=<where folders with individual species are> --wormbase_release_ftp_dir=<release tied to this WBPS version> --wbps_release_ftp_dir=<target directory> --repeats_dir=<repeats directory>";
die ("--source_dir not a directory: $source_dir . $usage") unless -d $source_dir;
die ("--wormbase_release_ftp_dir not a directory: $wormbase_release_ftp_dir . $usage") unless -d $wormbase_release_ftp_dir;
die ("--repeats_dir not a directory: $repeats_dir . $usage") unless -d $repeats_dir;
die ($usage) unless $wbps_release_ftp_dir;
die ($usage) unless $wbps_version;

for my $path_species (glob "$source_dir/*") {
  next if $sync_files_skip;
  my $species = basename $path_species;
  my ($spe, $cies) = split(/_/, $species);
  for my $this_source_dir ( glob "$source_dir/$species/*" ) {
     my $bioproject = basename $this_source_dir;
     my $this_target_dir = "$wbps_release_ftp_dir/species/$species/$bioproject";
     mkpath $this_target_dir if not -d $this_target_dir;
     my $putative_wormbase_dir = join("/", $wormbase_release_ftp_dir,"species", lc((substr $spe, 0, 1 ) . "_" . $cies) , uc($bioproject));
     if ( -d $putative_wormbase_dir and $wbps_version > 14 ) { # wait for WormBase to update C. remanei PX356
        print localtime ." ". $species . " making symlinks $putative_wormbase_dir -> $this_target_dir \n";
        &make_symlinks_to_wormbase_species (
          "$species.$bioproject.WBPS$wbps_version",
          $putative_wormbase_dir,
          $this_target_dir
        );
        my $cp_cmd = "rsync -a --include='*.paralogs.tsv.gz' --include='*.orthologs.tsv.gz' --include='*phenotypes.gaf.gz' --exclude '*' $this_source_dir/ $this_target_dir/";
        print localtime . " $species $cp_cmd\n";
        system($cp_cmd) and die("Failed: $cp_cmd");
     } else {
        my $cp_cmd = "rsync -a --include='*.gz' --exclude '*' $this_source_dir/ $this_target_dir/";
        print localtime . " $species $cp_cmd\n";
        system($cp_cmd) and die("Failed: $cp_cmd");
     }
    # find and copy repeat modeler2 files
    (my $bioproject_no_underscore = $bioproject) =~ s/_//g;
    my $repeats_file = $species."_".lc($bioproject_no_underscore)."-families.fa.gz";
    my $target_file_name = $species.".".uc($bioproject).".WBPS".$wbps_version.".repeat-families.fa.gz";	
    my $cp_repeats = "cp $repeats_dir/$repeats_file $this_target_dir/$target_file_name";
    print localtime . " $species $cp_repeats\n";
    system($cp_repeats) and die("Failed: $cp_repeats");
  }
}
unless ($checksums_skip){
  print localtime . " Remaking checksums file \n" ;
  system("cd $wbps_release_ftp_dir && find -L species -type f -exec md5sum \"{}\" + | sort -k 2,2 > CHECKSUMS") and die "Could not calc checksums\n";
}
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
