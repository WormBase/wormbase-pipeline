#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  submit_clustal.pl
#
#        USAGE:  ./submit_clustal.pl 
#
#  DESCRIPTION:  a script to submit clustalw jobs into LSF after cleaning
#                up the postgres database
#
#       AUTHOR:   (Michael Han), <mh6@sanger.ac.uk>
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  10/27/08 10:40:29 GMT
#     REVISION:  ---
#         NOTE:  * clustalw needs to be in the search path (most times in /software/worm/bin)
#                * the postgresql server needs to run
#===============================================================================

use lib '/software/worm/lib/site_perl';
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use DBI;
use strict;

my ($debug,$store,$species,$wb,$user,$pass,$test,$ace_database,$dontclean);
GetOptions(
 'debug=s'   => \$debug,
 'store=s'   => \$store,
 'species=s' => \$species,
 'user=s'    => \$user,
 'password=s'=> \$pass,
 'test=s'    => \$test,
 'acedb=s'   => \$ace_database,
 'dontclean' => \$dontclean,
) ||die(@!);

if ($store) {
    $wb = Storable::retrieve($store) or die('cannot load from storable');
}
else { $wb = Wormbase->new( -debug => $debug, -test => $test, -organism => $species) }

my $log = Log_files->make_build_log($wb);

# clean out the database
my $prefix = $wb->wormpep_prefix;
my $pgdb = DBI->connect('dbi:Pg:dbname=worm_clw;host=pgsrv1;port=5436',$user,$pass);
$pgdb->do("DELETE FROM clustal WHERE peptide_id LIKE \'$prefix\%\'") unless $dontclean;

# get wormpep dir and file name

my $infile = $wb->wormpep;
$infile.='/'.$wb->pepdir_prefix . 'pep'.$wb->version;
# my $acedb = $wb->autoace;
my $acedb = ($ace_database || glob('~wormpub/BUILD/autoace'));

# get wormpep number of proteins
my $count = `fgrep -c ">" $infile`;
chomp $count;
my $job = 1; 

$species ||= $wb->species;

# submit jobs, one for every 5k of sequences
while ($job <= ($count/5000)+1){
	my $script="/software/worm/perl_510/bin/perl ~/wormbase/scripts/clustal_runner.pl -window 5000 -offset $job -user $user -pass $pass -pepfile $infile -database $acedb -species $species";
	print `bsub -q long -e /dev/null -o /dev/null -J clw -R "select[mem>3000] rusage[mem=3000]" -M 3000000 $script`;
	$job++;
}

$log->mail;
