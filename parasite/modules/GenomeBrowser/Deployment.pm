use strict;
use warnings;
package GenomeBrowser::Deployment;
use LWP;
use File::Basename;
use Log::Any qw($log);

# Be in EBI
# Have tunnels enabled
# Then you can do everything

my $SANGER_HOST="sangerngs"; # Made up ssh alias

my $EBI_PATH="/nfs/ftp/pub/databases/arrayexpress/data/atlas/rnaseq";
my $EBI_URL="ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq";
my $SANGER_PATH="/data/production/parasites/rnaseqer";
my $SANGER_URL="https://ngs.sanger.ac.uk/production/parasites/rnaseqer";

sub location {
  my ( $root, $species, $assembly, $run_id ) = @_;
  (my $prefix = $run_id) =~ s/(.{6}).*/$1/;
  return join "/", $root, $species, $assembly, $prefix, "$run_id.bw";
}
sub run_in_sanger {
  my $cmd = "ssh $SANGER_HOST ".shift;
  my $output = `$cmd`;
  die $log->fatal("Failed: $cmd, output: $output") if $?;
}
sub sync_ebi_to_sanger {
  my ($species, $assembly, $run_id, $source_url, %opts) = @_;
  my $target_path = location ($SANGER_PATH,$species, $assembly, $run_id);
  my $target_url  = location ($SANGER_URL, $species, $assembly, $run_id);
  my $target_dir = dirname $target_path;

  if ($opts{do_sync} // not LWP::UserAgent->new->head($target_url)->is_success){
    $log->info("Initiating remote download: $source_url -> $SANGER_HOST:$target_path");
    run_in_sanger("mkdir -p $target_dir");
    run_in_sanger("wget --continue --no-verbose -O $target_path $source_url");
  } 
  return $target_url;
}
1;
