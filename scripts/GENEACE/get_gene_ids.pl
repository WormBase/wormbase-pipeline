#!/software/bin/perl -w

use lib $ENV{'CVS_DIR'};
use lib '/nfs/WWWdev/SANGER_docs/lib/Projects/C_elegans';
use lib '/software/worm/lib/perl';
use NameDB_handler;
use Wormbase;
use Storable;
use Getopt::Long;

#connect to name server and set domain to 'Gene'

my $PASS;
my $USER;
my $DOMAIN  = 'Gene';
my ($debug, $test, $store, $species, $database);

GetOptions (
	    "debug=s"    => \$debug,
	    "test"       => \$test,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "database:s" => \$database,
	    "user:s"	 => \$USER,
	    "pass:s"	 => \$PASS,
	   );


my $wormbase;

die "Species option is mandatory\n" unless $species;


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species
			     );
}
# establish log file.
my $log = Log_files->make_build_log($wormbase);
my $DB = $wormbase->test ? 'test_wbgene_id;utlt-db:3307' : 'wbgene_id;shap:3303';
my $db = NameDB_handler->new($DB,$USER,$PASS,'/nfs/WWWdev/SANGER_docs/htdocs/');
$db->setDomain($DOMAIN);

while (<>) {
	chomp;
	my @F=split;
	my $gene_id = $db->new_gene($F[0], 'Sequence', $species);
	print "$F[0] $gene_id\n";
}

$log->mail;
exit;
