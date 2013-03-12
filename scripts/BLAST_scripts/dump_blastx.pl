#!/usr/bin/env perl
#
# DESCRIPTION:
#  script to submit blastx dumping scripts onto the farm
#  and concatenate them at the end
# 
# Last edited by: $Author: gw3 $
# Last edited on: $Date: 2013-03-12 13:31:28 $
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
if (defined $ENV{'SANGER'}) {
  use lib '/software/worm/lib/site_perl';
} else {
  use lib "$ENV{'WORM_SW_ROOT'}/lib/perl5/site_perl"; 
}
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
	jappepX  => 'Japonica',
	brepepX  => 'Brenneri',
	GadflyX  => '1',
	ipi_humanX => '1',
	slimtremblX => '1',
	yeastX  => '1',
	slimswissprotX => '1',
);

my $m;
if (defined $ENV{'SANGER'}) {
  $m=LSF::JobManager->new(-q => 'normal',-o => '/dev/null',-e=>'/dev/null',-R => '"select[mem>4000] rusage[mem=4000]"',-M => 4000000, -F => 400000);
} else {
  $m=LSF::JobManager->new(-q => $ENV{'LSB_DEFAULTQUEUE'},-o => '/dev/null',-e=>'/dev/null',-R => '"select[mem>4000] rusage[mem=4000]"',-M => 4000, -F => 400000);
}

my $storable =  $wormbase->autoace . '/'. ref($wormbase).'.store';
$dumpdir ||= "$ENV{'PIPELINE'}/dumps";
my $organism = lc (ref($wormbase));

$database ||= "worm_ensembl_$organism";

my $name = 1;
my @outfiles;
# here goes the main bit:
foreach my $db(keys %logic2type){
  my @chroms = @{$wormbase->get_binned_chroms('5')}; # mo MtDNA
  $log->write_to("bsub commands . . . . \n\n");
  foreach my $chrom ( @chroms ) {
    $name ++;
    my $err = "/tmp/dump_blastx.$db.$name.err.$$";
    my $outfile = "$dumpdir/${organism}_$db.$name.ace";
    push @outfiles,$outfile;
    my $options="-database $database -logicname $db -outfile $outfile -store $storable -sequence $chrom";
    $options.=' -self' if $logic2type{$db} eq ref $wormbase; # set selfhit removal for the self-blasts
    my $cmd = "perl $ENV{'CVS_DIR'}/BLAST_scripts/blastx_dump.pl $options";
    $m->submit($cmd);
  }
}
$m->wait_all_children( history => 1 );
$log->write_to("All children have completed!\n");
for my $job ( $m->jobs ) {
  $log->error("Job $job (" . $job->history->command . ") exited non zero\n") if $job->history->exit_status != 0;
}
$m->clear;   

# check that all files end with a blank line, 
# otherwise the the job that created them was probably terminated prematurely by LSF
$log->write_to("\nTesting output ace files to see if they completed\n");
foreach my $file (@outfiles) {
   
    $log->write_to("$file ");
    my $endline = `tail -1 $file`;
    if (-e $file && $endline =~ /^\s*\n$/) {
      $log->write_to("- looks OK\n");
    } else {
      $log->write_to("- appears to be prematurely terminated.\nlast line in $file is:\n$endline\n");
      $log->error;
    }
  }

# concatenate the ace files into a big blob for later parsing with ensembl/ipi scripts
my $outfile="$dumpdir/${organism}_blastx.ace";
$log->write_to("Concatenating the ace files to create $outfile\n");

# in case of Elegans do something else
if ($wormbase->species eq 'elegans' or $wormbase->species eq 'briggsae'){
  my @files = glob("$dumpdir/$organism*X.*.ace");
  unlink $outfile if -e $outfile;
  foreach my $file (@files ){
    $log->write_to("\tcat $file\n");
    system ("cat $file | perl $ENV{'CVS_DIR'}/BLAST_scripts/convert_chromblast2clone.pl >> $outfile") 
      && die("cannot concatenate $file to $outfile\n" );
  }
} else {
  system ("cat $dumpdir/$organism*X*.ace > $outfile") && die("cannot concatenate dumpdir/$organism*X.ace to $outfile\n" );
}


$log->mail();
								    
