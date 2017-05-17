#!/software/bin/perl -w

use strict;
use Getopt::Long;
use Carp;
use DBI;
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;

$| = 1;

my $db = "worm_ensembl_elegans";
my $user = 'wormro';
my $host = $ENV{'WORM_DBHOST'};
my $port = $ENV{'WORM_DBPORT'};
my $segsize = 1000000;

my ($dump_loc, $dump_one_script, $out_file_prefix);
my $resource;
my $test;

&GetOptions(
	    'host:s'       => \$host,
	    'user:s'       => \$user,
	    'db:s'         => \$db,
	    'port:n'       => \$port,
	    'segsize:n'    => \$segsize,
	    'resource:s'   => \$resource,
	    'test'         => \$test,
	    'dump_script:s'=> \$dump_one_script,
            'dumploc:s'    => \$dump_loc,
            'prefix:s'     => \$out_file_prefix,
	   );

croak("You must supply a valid location for the dump") 
    if not defined $dump_loc or not -d $dump_loc;
croak("You must supply an executable script for dumping one batch with -dump_script")
    if not defined $dump_one_script or not -x $dump_one_script;

my $dbh=DBI->connect("dbi:mysql:database=$db;host=$host;port=$port", "wormro");

my $sth = $dbh->prepare('select count(*) from protein_feature where analysis_id in (select analysis_id from analysis where module="BlastPep")');
$sth->execute;
my ($nrow) = $sth->fetchrow;


my $nseg = int(($nrow/$segsize))+1;

my $job_name = "worm_${db}_dump";

my $lsf=LSF::JobManager->new(-q => $ENV{'LSB_DEFAULTQUEUE'},
#                             -P => 'wormbase', 
			     -R => "select[mem>4000] rusage[mem=4000]", 
			     -M => 4000, 
#			     -F => 2000000, 
			     -J => $job_name);

for (my $i=0; $i<$nseg; $i++) {
  my $start = $i*$segsize;

  my @bsub_options = (-o => "junk$i.log", -e => "junk$i.err");

  my $cmd = "perl $dump_one_script -host $host -user $user -port $port -db $db -start $start -count $segsize -out $dump_loc/$out_file_prefix.$i.srt";
  print "$cmd\n";

  $lsf->submit(@bsub_options, $cmd);

  last if $test;
}

$lsf->wait_all_children( history => 1 );
print "All children have completed!\n";

for my $job ($lsf->jobs){ # much quicker if history is pre-cached
  print "Job $job (" . $job->history->command . ") exited ". $job->history->exit_status ."\n" if ($job->history->exit_status != 0);
}
$lsf->clear; # clear out the job manager to reuse.

exit(0);
