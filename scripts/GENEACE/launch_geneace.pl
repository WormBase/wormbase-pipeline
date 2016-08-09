#!/software/bin/perl -w
#
# A little script to check for write access on geneace and if it is free 
# launch geneace with the custom ry mt3 uses.

use strict;
use lib $ENV{'CVS_DIR'};
use Getopt::Long;
use File::stat;
use Time::localtime;
use Wormbase;
use Storable;
my $c;
my $debug;

GetOptions(
	   'debug:s' => \$debug,
	   'c:s'  => \$c,
	  );

  my $wormbase = Wormbase->new( -debug => $debug,
                           );

my $file = "/nfs/wormpub/DATABASES/geneace/database/lock.wrm";
if (-e $file) {
  print "\nTesting for write access.\n\nStatus: Busy ";
  my $timestamp = ctime(stat($file)->mtime);
  print "(Write access taken $timestamp)\nBy Session User - ";
  open my $filehandle, "<", $file or die "could not open $file: $!";
    print while <$filehandle>;
  exit;
} 
else {
  print "\nTesting for write access.\n\nStatus: Free\n\n";
$wormbase->run_command("/software/worm/acedb/current/bin/xace -nosplash -tsuser mt3 ~wormpub/DATABASES/geneace",);
}
exit;
