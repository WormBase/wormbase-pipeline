use strict;
use warnings;
package GenomeBrowser::Deployment;
use LWP;
use File::Basename;
use Log::Any qw($log);
use Try::Tiny;
use File::Temp qw/tempdir/;

# Be in EBI

# For embassy:
# load module embassy

# For sanger:
# Enable tunnels

my $SANGER_HOST="sangerngs"; # Made up ssh alias

my $EBI_PATH="/nfs/ftp/pub/databases/arrayexpress/data/atlas/rnaseq";
my $EBI_URL="ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq";
my $SANGER_PATH="/data/production/parasites/rnaseqer";
my $SANGER_URL="https://ngs.sanger.ac.uk/production/parasites/rnaseqer";
my $EMBASSY_PATH=$ENV{EMBASSY_RNASEQER_PATH};
my $EMBASSY_URL=$ENV{EMBASSY_ACCESS_URL_RNASEQER};
my $EMBASSY_BASE_URL=$ENV{EMBASSY_URL};
my $EMBASSY_COMMAND=$ENV{EMBASSY_COMMAND};
my $tmpdir = tempdir(CLEANUP => 1);

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

sub solve_ftp_file_format {
    my ($path, $tmp) = @_;
    my $bzpath = $path.".tar.bz2";
    my $fh;
    if (-e $path) {
       return $path
    } elsif (-e $bzpath ) {
        my $tmpfile_basename = `tar -xvf ${bzpath} -C ${tmp}`;
        chomp($tmpfile_basename);
        my $tmpfile = "${tmp}/${tmpfile_basename}";
        if (-e $tmpfile) {
            return $tmpfile
        } else {
            $log->fatal("Failed: $tmpfile does not exist");
            return "ERROR";
        }
    } else {
        $log->fatal("Failed: Could not resolve $path or $bzpath");
        return "ERROR";
    }
}


sub sync_ebi_externally {
    my ($species, $assembly, $run_id, $source_url, %opts) = @_;
    if ($opts{deploy_to_sanger}) {
        return(sync_ebi_to_sanger($species, $assembly, $run_id, $source_url, %opts));
    }
    return(sync_ebi_to_embassy($species, $assembly, $run_id, $source_url, %opts));
}

sub sync_ebi_to_sanger {
  my ($species, $assembly, $run_id, $source_url, %opts) = @_;
  my $target_path = location ($SANGER_PATH, $species, $assembly, $run_id);
  my $target_url  = location ($SANGER_URL, $species, $assembly, $run_id);
  my $target_dir = dirname $target_path;

  unless ($opts{deployment_skip} or LWP::UserAgent->new->head($target_url)->is_success){
    $log->info("Initiating remote download: $source_url -> $SANGER_HOST:$target_path");
    run_in_sanger("mkdir -p $target_dir");
    try {
      run_in_sanger("wget --continue --no-verbose -O $target_path $source_url");
    } catch {
        (my $path = $source_url) =~ s/ftp:\/\/ftp\.ebi\.ac\.uk\/pub/\/nfs\/ftp\/public/;
        $log->info("Trying scp $path $SANGER_HOST:$target_path");
        my $cp_cmd = "scp $path $SANGER_HOST:$target_path";
        my $output = `$cp_cmd`;
        die $log->fatal("Failed: $cp_cmd, output: $output") if $?;
    }
  } 
  return $target_url;
}

sub sync_ebi_to_embassy {
  my ($species, $assembly, $run_id, $source_url, %opts) = @_;
  my $target_path = location ($EMBASSY_PATH, $species, $assembly, $run_id);
  my $target_url  = location ($EMBASSY_URL, $species, $assembly, $run_id);
  my $target_dir = dirname $target_path;
  unless ($opts{deployment_skip} or LWP::UserAgent->new->head($target_url)->is_success or index($source_url, $EMBASSY_BASE_URL) != -1) {
      $log->info("deployment: $opts{deployment_skip}");
      $log->info("Initiating remote download: $source_url -> $target_path");
      (my $path = $source_url) =~ s/ftp:\/\/ftp\.ebi\.ac\.uk\/pub/\/nfs\/ftp\/public/;
      my $resolved_path = solve_ftp_file_format($path, $tmpdir);
      unless ($resolved_path eq "ERROR") {
          $log->info("Trying $EMBASSY_COMMAND s3 scp $resolved_path $target_path");
          my $cp_cmd = "$EMBASSY_COMMAND s3 cp $resolved_path $target_path";
          my $output = `$cp_cmd`;
          die $log->fatal("Failed: $cp_cmd, output: $output") if $?;
      }
  }
  if (index($source_url, $EMBASSY_BASE_URL) != -1){
      return $source_url;
  } else {
      return $target_url;
  }
}


1;
