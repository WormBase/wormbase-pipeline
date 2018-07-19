
package GenomeBrowser::Deployment;
use LWP;
use File::Basename;

# Be in EBI
# Have tunnels enabled
# Then you can do everything

our $SANGER_HOST="sangerngs"; # Made up ssh alias

our $EBI_PATH="/nfs/ftp/pub/databases/arrayexpress/data/atlas/rnaseq";
our $EBI_URL="ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq";
our $SANGER_PATH="/data/production/parasites/rnaseqer";
our $SANGER_URL="http://ngs.sanger.ac.uk/production/parasites/rnaseqer";

sub location {
 my $run_id = shift;
 (my $prefix = $run_id) =~ s/(.{6}).*/$1/;
 return join "/", $prefix, $run_id, "$run_id.bw";
}

sub file_present_at_ebi {
  my $path = join "/", $EBI_PATH, location(@_);
  return -f $path;
}

sub sync_ebi_to_sanger {
  my $source_url = join "/", $EBI_URL, location(@_);
  my $target_path = join "/", $SANGER_PATH, location(@_);
  my $target_dir = dirname $target_path;
  print qx/ssh $SANGER_HOST mkdir -p $target_dir/;
  print qx/ssh $SANGER_HOST wget --continue --no-verbose -O $target_path $source_url/; 
}

sub file_present_at_sanger {
  my $path = join "/", $SANGER_PATH, location(@_);
  system("ssh $SANGER_HOST ls $path > /dev/null 2>&1");
  return 0 == $?;
}

sub file_is_online {
  my $path = join "/", $SANGER_URL, location(@_);
  return LWP::UserAgent->new->head($path)->is_success;
}
1;
