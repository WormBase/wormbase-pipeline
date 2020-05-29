#!/software/bin/perl -w
use lib $ENV{'CVS_DIR'};
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
#use lib '/nfs/users/nfs_g/gw3/Nameserver-API';

use lib '/software/worm/lib/perl';

use NameDB_handler;
use Wormbase;
use Storable;
use Getopt::Long;

my ($help, $debug, $verbose, $store, $wormbase, $species);

GetOptions (
	    'debug=s'    => \$debug,
	    'store=s'    => \$store,
	    'species=s'  => \$species,
	   )||die();
$species || die "Species option is mandatory\n";

my $wormbase;
if ($store) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -organism => $species
			   );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

print "The new nameserver REST API will not allow you to simply make a new Gene ID without a CGC name or a sequence name.\nTry using 'batch_gene -action new' instead.\n\n";

$log->mail;
