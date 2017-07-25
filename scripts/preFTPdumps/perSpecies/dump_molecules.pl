#!/usr/bin/env perl
# dumps the whole molecule class into the reports directory

use strict;

use IO::File;
use Storable;
use Getopt::Long;

use lib $ENV{CVS_DIR};

use Wormbase;
use Log_files;

my ($store,$debug,$test,$database,$species,$outfile);
GetOptions(
       'store=s' => \$store,
       'debug=s' => \$debug,
       'test'    => \$test,
       'species=s'  => \$species,
       'database=s' => \$database,
       'outfile=s'  => \$outfile,
)||die(@!);

my $wormbase;
if ($store) { $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")} 
else {$wormbase = Wormbase->new( -debug => $debug, 
                                 -test => $test,
                                 -organism => $species);
}

my $log = Log_files->make_build_log($wormbase);

$database = $wormbase->autoace if not defined $database;

$outfile = $wormbase->reports . '/molecules.ace' if not defined $outfile;
unlink $outfile if -e $outfile;

$log->write_to("writing to $outfile\n");

open TACE, '|'.$wormbase->tace ." $database";
print TACE "Find Molecule\nShow -a -f $outfile\nquit\n";

close TACE;

$log->mail;
exit(0);
