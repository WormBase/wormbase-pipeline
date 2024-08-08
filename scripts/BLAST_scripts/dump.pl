#!/software/bin/perl -w

use strict;
use Getopt::Long;
use Carp;
use DBI;
use Modules::WormSlurm;
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

my $sth = $dbh->prepare('select count(*) from protein_feature where analysis_id in (select analysis_id from analysis where program="blastp")');
$sth->execute;
my ($nrow) = $sth->fetchrow;


my $nseg = int(($nrow/$segsize))+1;

my $job_name = "worm_${db}_dump";

my %slurm_jobs;
for (my $i=0; $i<$nseg; $i++) {
  my $start = $i*$segsize;

  my $cmd = "perl $dump_one_script -host $host -user $user -port $port -db $db -start $start -count $segsize -out $dump_loc/$out_file_prefix.$i.srt";
  print "$cmd\n";

  my $job_id = WormSlurm::submit_job_with_name($cmd, 'production', '4g', '4:00:00', "junk$i.log", "junk$i.err");
  $slurm_jobs{$job_id} = $cmd;

  last if $test;
}

WormSlurm::wait_for_jobs(keys %slurm_jobs);
print "All children have completed!\n";

for my $job_id (keys %slurm_jobs){
    my $exit_code = WormSlurm::get_exit_code($job_id);
    print "Slurm job $job_id (" . $slurm_jobs{$job_id} . ") exited " . $exit_code ."\n" if ($exit_code != 0);
}

exit(0);
