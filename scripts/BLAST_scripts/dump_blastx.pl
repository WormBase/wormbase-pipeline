#!/usr/bin/env perl
#
# DESCRIPTION:
#  script to submit blastx dumping scripts onto the farm
#  and concatenate them at the end
# 
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2014-04-24 10:28:50 $
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
									
use strict;
use Getopt::Long;

use lib $ENV{CVS_DIR};
use threads;
use Thread::Queue;
use Wormbase;
use Coords_converter;



my ($database,$store,$debug,$species,$test,$dumpdir);
GetOptions(
  'database=s'  => \$database,
  'store=s'     => \$store,
  'debug=s'     => \$debug,
  'species=s'   => \$species,
  'test'        => \$test,
  'dumpdir=s'   => \$dumpdir,
) || die($usage);

## create the worker threads before the log object is created and copied around
my @workers;
my $queue=Thread::Queue->new;

sub blastx_worker{
 while (my $string = $queue->dequeue){
   system($string) && die(@!);
 }
}


for (1..12){ # 12 methods/threads
  push @workers,threads->create('blastx_worker');
}

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
$species = $wormbase->species;

# that might work until we change logic_names
my %logic2type = (
  wormpepx       => 'elegans',
  ppapepx        => 'pristionchus',
  jappepx        => 'japonica',
  brugpepx       => 'brugia',
  brigpepx       => 'briggsae',
  remapepx       => 'remanei',
  brepepx        => 'brenneri',
  ovolpepx       => 'ovolvulus',
  gadflyx        => '1',
  ipi_humanx     => '1',
  yeastx         => '1',
  slimswissprotx => '1',
  );

my $storable =  $wormbase->autoace . '/'. ref($wormbase).'.store';
$dumpdir ||= "$ENV{PIPELINE}/dumps";
$database ||= "worm_ensembl_${species}";

my @outfiles;

foreach my $db (keys %logic2type){
  my $outfile = "$dumpdir/${species}_$db.1.ace";
  unlink $outfile if -e $outfile; # cleanup leftover bits
  push @outfiles,$outfile;
  my $options="-database $database -logicname $db -outfile $outfile -store $storable -toplevel";
  $options .= ' -self' if $logic2type{$db} eq $species;
  my $cmd = "perl $ENV{CVS_DIR}/BLAST_scripts/blastx_dump.pl $options";
  $queue->enqueue($cmd);
}
$queue->enqueue(undef) for 1..12; # 12 methods

map {$_->join} @workers; # wait for the worker threads to finish

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
my $outfile="$dumpdir/${species}_blastx.ace";
unlink $outfile if -e $outfile;
$log->write_to("Concatenating the ace files to create $outfile\n");

# in case of Elegans do something else
if ($wormbase->species eq 'elegans' or 
    $wormbase->species eq 'briggsae'){
  my @files = glob("$dumpdir/$species*x.*.ace");

  # the following will trigger the generation of SL_coords and clone_coords files if they do not exists
  Coords_converter->invoke($wormbase->orgdb, undef, $wormbase);

  open(my $out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");

  foreach my $file (@files ){
    $log->write_to("\tcat $file\n");
    open(my $conv_fh, "perl $ENV{CVS_DIR}/BLAST_scripts/convert_chromblast2clone.pl -species $species $file |")
        or $log->log_and_die("Could not open convert_chromblast2clone command\n");
    while(<$conv_fh>) {
      print $out_fh $_;
    }
  }
  close($out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
} else {
  $wormbase->run_command("cat $dumpdir/$species*x*.ace > $outfile", $log);
}

$log->mail();								    
