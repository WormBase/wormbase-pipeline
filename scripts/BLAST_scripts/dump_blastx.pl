#!/usr/bin/env perl
#
# DESCRIPTION:
#  script to submit blastx dumping scripts onto the farm
#  and concatenate them at the end
# 
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2008-08-08 13:41:36 $
# 


my $usage = <<USAGE;
dump_blastx.pl options:
        -debug USER_NAME    sets email address and debug mode
        -store FILE_NAME    use a Storable wormbase configuration file
	-database DB_NAME   Ensembl datbase to dump from. If not used deafults to the database of the storable
	-species SPECIES_NAME which species you want to dump (a.e. elegans, briggsae,....)
	-test               if you want to use TEST_BUILD instead of BUILD
        -dumpdir DIRECTORY_NAME id you want to dump it to a different directory
USAGE
									

use Getopt::Long;
use lib $ENV{CVS_DIR};
use lib '/software/worm/lib/site_perl';
use LSF;
use LSF::JobManager;
use Wormbase;
use strict;


my ($database,$store,$debug,$species,$test,$dumpdir);
GetOptions(
	'database=s'  => \$database,
	'store=s'     => \$store,
	'debug=s'     => \$debug,
	'species=s'   => \$species,
	'test'        => \$test,
	'dumpdir=s'   => \$dumpdir,
) || die($usage);


my $wormbase;
if ($store) {
    $wormbase = retrieve($store) or croak("Can't restore wormbase from $store\n");
} else {
    $wormbase = Wormbase->new(
        -debug    => $debug,
        -test     => $test,
	-organism => $species,
    );
}


my $log = Log_files->make_build_log($wormbase);

# that might work until we change logic_names
my %logic2type = (
	remaneiX => 'Remanei',
	brigpepX => 'Briggsae',
	wormpepX => 'Elegans',
	ppapepX  => 'Pristionchus',
	GadflyX  => '1',
	ipi_humanX => '1',
	slimtremblX => '1',
	yeastX  => '1',
	slimswissprotX => '1',
);

my $m=LSF::JobManager->new(-q => 'normal',-o => '/dev/null',-e=>'/dev/null',-R => '"select[mem>4000] rusage[mem=4000]"',-M => 4000000, -F => 400000);

my $storable =  $wormbase->autoace . '/'. ref($wormbase).'.store';
$dumpdir ||= '/lustre/work1/ensembl/wormpipe/dumps';
my $organism = lc (ref($wormbase));

$database ||= "worm_ensembl_$organism";

# here goes the main bit:
foreach my $db(keys %logic2type){
	my $options="-database $database -logicname $db -outfile $dumpdir/${organism}_$db.ace -store $storable";
	$options.=' -self' if $logic2type{$db} eq ref $wormbase; # set selfhit removal for the self-blasts
	$options.=' -toplevel';# unless ref $wormbase eq 'Elegans'; # elegans dumps on clone level
	$m->submit("/software/bin/perl $ENV{CVS_DIR}/BLAST_scripts/blastx_dump.pl $options");
}

$m->wait_all_children( history => 1 );
print "All children have completed!\n";

for my $job ($m->jobs){ # much quicker if history is pre-cached
       $log->write_to("$job exited non zero\n") if $job->history->exit_status != 0;
}
$m->clear; # clear out the job manager to reuse.

# concatenate the ace files into a big blob for later parsing with ensembl/ipi scripts
my $outfile="$dumpdir/${organism}_blastx.ace";

# in case of Elegans do something else
if (ref $wormbase eq 'Elegans'){
  my @files = glob("$dumpdir/$organism*X.ace");
  unlink $outfile if -e $outfile;
  foreach my $file (@files ){
          system ("cat $file |/software/bin/perl $ENV{CVS_DIR}/BLAST_scripts/convert_chromblast2clone.pl >> $outfile") 
                 && die("cannot concatenate $file to $outfile\n" );
	 }
} else {
   system ("cat $dumpdir/$organism*X.ace > $outfile") && die("cannot concatenate dumpdir/$organism*X.ace to $outfile\n" );
}

# $wormbase->run_command("rm -f $dumpdir/$organism*X.ace",$log);

$log->mail();
								    
