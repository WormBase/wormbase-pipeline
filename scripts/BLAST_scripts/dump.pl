#!/usr/local/ensembl/bin/perl -w

use strict;
use Getopt::Long;
use Carp;
use DBI;

$| = 1;

my $db = "worm_ensembl_elegans";
my $user = 'wormro';
my $host = 'ia64d';
my $port = 3306;
my $segsize = 1000000;
my $resource;
#my $resource = 'model=ES40_667';
my $machines = undef;
my $test;
my $dump_one_script;

&GetOptions(
	    'host:s'       => \$host,
	    'user:s'       => \$user,
	    'db:s'         => \$db,
	    'port:n'       => \$port,
	    'segsize:n'    => \$segsize,
	    'machines:s'   => \$machines,
	    'resource:s'   => \$resource,
	    'test'         => \$test,
	    'dump_script:s'  => \$dump_one_script,
	   );
$dump_one_script = $dump_one_script?glob($dump_one_script):'/lustre/work1/ensembl/wormpipe/script/dump_one_new.pl';

my $dbh=DBI->connect( "dbi:mysql:database=$b;host=$host;port=$port", "wormro");

my $sth = $dbh->prepare('select count(*) from protein_feature where analysis_id in (select analysis_id from analysis where module="BlastPep")');
$sth->execute;
my ($nrow) = $sth->fetchrow;


my $nseg = int(($nrow/$segsize))+1;

print "N seg = $nseg\n";

for (my $i=0; $i<$nseg; $i++) {
  my $start = $i*$segsize;
  print "start = $start\n";

 # my $comstr = "bsub -P wormbase -q normal -R $resource -o junk$i.log ";
  my $comstr = "bsub -P wormbase -q normal -o junk$i.log "; 
  $comstr .= "-m \"$machines\" " if ($machines); 
  $comstr .= "$dump_one_script -host $host -user $user -port $port -db $db -start $start -count $segsize -out junk$i.srt";

  print $comstr . "\n";

  system($comstr);
  last if $test;
}
