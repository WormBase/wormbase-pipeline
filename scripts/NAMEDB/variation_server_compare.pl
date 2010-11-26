#!/software/bin/perl -w

use lib '/nfs/WWWdev/SANGER_docs/lib/Projects/C_elegans/';
use lib $ENV{'CVS_DIR'};

use strict;
use NameDB_handler;
use Wormbase;
use Log_files;
use Ace;
use Carp;
use Getopt::Long;
use Storable;

my ($help, $debug, $test, $store, $database);

GetOptions (    'help'       => \$help,
                'debug=s'    => \$debug,
                'test'       => \$testt,
                'store:s'    => \$store,
                'database:s' => \$database,
           );

my $wormbase;
if ($store) {
    $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
    $wormbase = Wormbase->new(-debug => $debug, -test => $test);
}


# establish log file.
my $log = Log_files->make_build_log($wormbase);

#connect to database and read in data
my $acedb = ($database or $wormbase->database('geneace'));
$log->write_to("Checking $acedb for errors\n");

my $def = "$ENV{CVS_DIR}/wquery/geneace/variation_nameserver_comm.def";

my $TABLE = $wormbase->table_maker_query($acedb, $def);
my %ace_ids;
while(<$TABLE>) {
    next unless /WBVar/;
    s/\"//g;
    my ($id,$name) = split;
    $ace_ids{$id} = $name;
}
delete $ace_ids{'acedb>'};
close $TABLE;

#connect to name server and set domain to 'Gene'
my $DB      = 'wbgene_id;shap;3303';
my $PASS    = "wormpub";
my $USER    = "wormpub";
my $DOMAIN  = 'Variation';
my $db = NameDB_handler->new($DB,$USER,$PASS,'/nfs/WWWdev/SANGER_docs/htdocs');

# get nameserver data
my $query = 'SELECT primary_identifier.object_public_id, 
                    secondary_identifier.object_name
             FROM primary_identifier,secondary_identifier 
             WHERE  primary_identifier.object_id = secondary_identifier.object_id 
             AND   (secondary_identifier.name_type_id = 5) 
             AND   (primary_identifier.object_live = 1) 
             ORDER BY object_public_id';

#RESULTS . . 
#+------------------+-------------+
#| object_public_id | object_name |
#+------------------+-------------+
#| WBVar00000001    | a83         |
#| WBVar00000002    | ac1         |
#| WBVar00000003    | act-123     |
#+------------------+-------------+

my $sth = $db->dbh->prepare($query);
$sth->execute();

my %server_ids;
while (my ( $var, $name ) = $sth->fetchrow_array){
    $server_ids{$var} = $name;
}

foreach my $var (keys %ace_ids) {
    if($server_ids{$var}){
      if($server_ids{$var} eq $ace_ids{$var}) {
        #print "$var ok\n";
      }
      else {
        $log->error("$var - acedb = $ace_ids{$var} : var-server = $server_ids{$var}\n");
      }
      delete($server_ids{$var});
    } else {
        $log->error("$var missing from var-server\n");
    }
}

foreach my $var (keys %server_ids) {
    $log->error("$var missing Public_name or absent from acedb\n");
}

$log->write_to("Work Done!\n");
$log->mail();
exit(0);