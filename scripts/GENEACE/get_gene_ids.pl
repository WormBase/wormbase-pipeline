#!/software/bin/perl -w
use lib $ENV{'CVS_DIR'};
use lib '/software/worm/lib/perl';

use NameDB_handler;
use Wormbase;
use Storable;
use Getopt::Long;

# connect to name server and set domain to 'Gene'
my $DOMAIN  = 'Gene';
my ($PASS,$USER,$debug, $test, $store, $species);

GetOptions (
	    'debug=s'    => \$debug,
	    'test'       => \$test,
	    'store=s'    => \$store,
	    'species=s'  => \$species,
	    'user=s'	 => \$USER,
	    'pass=s'	 => \$PASS,
	   )||die();
$species || die "Species option is mandatory\n";

my $wormbase;
if ($store) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);
my $DB = $wormbase->test ? 'test_wbgene_id;utlt-db:3307' : 'nameserver_live;web-wwwdb-core-02:3449';
my $db = NameDB_handler->new($DB,$USER,$PASS);
$db->setDomain($DOMAIN);

while (<>) {
	chomp;
	my @F=split;
	my $gene_id = $db->new_gene($F[0], 'Sequence', $species);
	print "$F[0] $gene_id\n";
}
$log->mail;
