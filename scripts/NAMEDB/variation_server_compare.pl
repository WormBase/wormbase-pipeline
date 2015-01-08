#!/software/bin/perl -w

use lib $ENV{'CVS_DIR'};
use lib "$ENV{CVS_DIR}/NAMEDB/lib";

use strict;
use NameDB_handler;
use Wormbase;
use Log_files;
use Ace;
use Carp;
use Getopt::Long;
use Storable;

my ($help, $debug, $test, $store, $database, $paranoid);

GetOptions (    'help'       => \$help,
                'debug=s'    => \$debug,
                'test'       => \$test,
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

my $def = "$ENV{CVS_DIR}/../wquery/geneace/variation_nameserver_comm.def";

my $TABLE = $wormbase->table_maker_query($acedb, $def);
my %ace_ids;
while(<$TABLE>) {
  next unless /WBVar/;
  s/\"//g;
  my ($id,$status,$name) = split;
  
  $ace_ids{$id} = {
    name => $name,
    status => $status,
  };
}
close($TABLE) or $log->log_and_die("Could not successfully complete TM query\n");

#connect to name server and set domain to 'Gene'
my $DB      = 'wbgene_id;shap;3303';
my $PASS    = 'wormpub';
my $USER    = 'wormpub';
my $DOMAIN  = 'Variation';
my $db = NameDB_handler->new($DB,$USER,$PASS,$wormbase->wormpub . '/DATABASES/NameDB');

# get nameserver data
my $query = 'SELECT primary_identifier.object_public_id, 
                    secondary_identifier.object_name,
                    primary_identifier.object_live
             FROM primary_identifier,secondary_identifier 
             WHERE  primary_identifier.object_id = secondary_identifier.object_id 
             AND   (secondary_identifier.name_type_id = 5) 
             ORDER BY object_public_id';

my $sth = $db->dbh->prepare($query);
$sth->execute();

my %server_ids;
while (my ( $var, $name, $live ) = $sth->fetchrow_array){
  $server_ids{$var} = {
    name => $name,
    status => ($live) ? 'Live' : 'Dead',
  };
}

my %ids = map { $_ => 1 } (keys %server_ids, keys %ace_ids);

foreach my $var (sort keys %ids) {
  if (exists $server_ids{$var} and not exists $ace_ids{$var}) {
    if ($server_ids{$var}->{status} eq 'Live') {
      $log->error(sprintf("%s (%s in server) missing from Acedb\n", 
                          $var, 
                          $server_ids{$var}->{status}));
    }
  } elsif (exists $ace_ids{$var} and not exists $server_ids{$var}) {
    if ($ace_ids{$var}->{status} eq 'Live') {
      $log->error(sprintf("%s (%s in Acedb) missing from Server or has no secondary namer\n", 
                          $var, 
                          $ace_ids{$var}->{status}));
    }
  } else {
    # must be present in both

    if(exists $server_ids{$var}){
      next if $ace_ids{$var}->{status} eq 'Dead' and $server_ids{$var}->{status} eq 'Dead';
      
      if ($server_ids{$var}->{name} ne $ace_ids{$var}->{name}) {
        $log->error(sprintf("%s (%s in Acedb) public names disagree: acedb = %s : var-server = %s\n", 
                            $var,
                            $ace_ids{$var}->{status},
                            $ace_ids{$var}->{name},
                            $server_ids{$var}->{name}));
        
      } elsif ($ace_ids{$var}->{status} ne $server_ids{$var}->{status} and
               ($server_ids{$var}->{status} eq 'Dead' and $ace_ids{$var}->{status} eq 'Live' or
                $server_ids{$var}->{status} eq 'Live' and $ace_ids{$var}->{status} eq 'Dead')) {
        
        $log->error(sprintf("%s acedb status = %-10s : var-server status = %-10s\n", 
                            $var,
                            #$ace_ids{$var}->{name},
                            $ace_ids{$var}->{status}, 
                            $server_ids{$var}->{status}, 
                            ));
        
      }
    }  
  }
}
  
  
$log->write_to("Work Done!\n");
$log->mail();
exit(0);
