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

my ($help, $debug, $test, $store, $database, $def, $dump,$old,);

my $ndb = 'nameserver_live';
my $nhost = 'web-wwwdb-core-02';
my $nport = 3449;
my $nuser = 'wormpub';
my $npass = 'wormpub';

GetOptions (	"help"       => \$help,
            	"debug=s"    => \$debug,
		"test"       => \$test,
		"store:s"    => \$store,
		"database:s" => \$database,
		"def:s"      => \$def,
                "ndb:s"      => \$ndb,
                "nhost:s"    => \$nhost,
                "nport:s"    => \$nport,
                "nuser:s"    => \$nuser,
                "npass:s"    => \$npass,
                "old" => \$old,
                'dumpfile=s' => \$dump,     # file dumped from the datomic Name service check the confluence guide on how to produce this
		);

my $wormbase;
if ( $store ) {
    $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
    $wormbase = Wormbase->new( -debug   => $debug,
			       -test    => $test,
			       );
}

my (%ace_genes, %server_genes, %all_gene_ids);

my %exceptions = ('C05G5.6'=>1, 'T22C1.13'=>1, 'F22D3.9'=>1, 'C54G10.6'=>1, 'C11H1.12'=>1, 'T01B6.12'=>1, 'C40A11.10'=>1);

log->log_and_die("you need to specify -dumpfile and point it to a datomic nameservice dump file\n") unless ($dump);

# establish log file.
my $log = Log_files->make_build_log($wormbase);

#connect to database and read in data
my $acedb = ($database or $wormbase->database('camace'));

$log->write_to("Checking $acedb for errors against $ndb;$nhost;$nport\n");



$def = $wormbase->database('camace')."/wquery/SCRIPT:camace_nameDB_comm.def" unless $def;
my $TABLE = $wormbase->table_maker_query($acedb, $def);


while( <$TABLE> ){
  next if (/>/ or /\/\// );
  s/\"//g;  # remove "
  my($gene, $cds, $transcript, $pseudo) = split(/\s/);  

  next unless ($cds or $transcript or $pseudo);
  $log->write_to("WARNING: tablemaker returned $_") unless $gene=~/WBGene/;
  
  #$ace_genes{"$gene"}->{'cds'}->{"$cds"}        = 1 if $cds;
  my $seq_name = ($cds or $transcript or $pseudo);
  $seq_name =~ s/[a-z]$//; #remove isoform indication
  if( $ace_genes{$gene}->{name} and ($ace_genes{$gene}->{name} ne $seq_name) ) {
    my $testname = $ace_genes{"$gene"}->{'name'};
    unless ($exceptions{$testname}) {
#    unless ($ace_genes{"$gene"}->{'name'} eq 'C05G5.6') {
      $log->write_to("$gene has multiple sequence names ".$ace_genes{"$gene"}->{'name'}." and $seq_name\n");
    }
    next;
  }
  else {
    $ace_genes{$gene}->{name} = $seq_name;
    $ace_genes{$gene}->{status} = ($cds or $transcript or $pseudo) ?  1 : 0; #live if it has one these nametypes
  }
}

#Retrieve the data from the NS dump.
unless ($old) {my ($dump_fh);
	       $log->write_to("using $dump file for checking\n");
               open($dump_fh, "<$dump") or $log->log_and_die("Can't open $dump\n");
               my $gene;
               my $name;
               while ( my $line = <$dump_fh> ) {
                   chomp $line;
                   my @f = split ",", $line;
                   unless ($f[0] =~ /WBGene/) {next;}
                   $gene = $f[0];  
                   if ($f[2] =~ /\S+/) {
                       $name = $f[2];
                       $server_genes{$gene}->{name} = $name;
                   }
                   if ($f[3] =~ /live/) {
                       $server_genes{$gene}->{'status'} = "1";
                   }
                   else {
                       $server_genes{$gene}->{'status'} = "0";
                   }    
               }
}

#connect to name server and set domain to 'Gene'
if ($old) {
my $db;
eval {
  $db = NameDB_handler->new("$ndb;$nhost;$nport",
                             $nuser,
                             $npass, 
                             $wormbase->wormpub . "/DATABASES/NameDB");
};
$@ and do {
  $log->log_and_die("Could not connect to the NameDB : $@\n");
};


my $dom_id = $db->getDomainId('Gene');

# get nameserver data
my $query = "SELECT pi.object_public_id, 
                    pi.object_live, 
                    si.object_name
             FROM primary_identifier pi, secondary_identifier si, name_type nt
             WHERE  pi.object_id = si.object_id 
             AND    si.name_type_id = nt.name_type_id
             AND    nt.name_type_name = 'Sequence'
             AND    pi.domain_id = $dom_id
             ORDER BY pi.object_public_id";

# results
#| object_public_id |  live       | name
#| WBGene00044331   |           1 | T20B3.15      | 
#| WBGene00044331   |           1 | T20B3.15      |

print "$query\n";

my $sth = $db->dbh->prepare($query) || $log->log_and_die($db->dbh->errstr);
$sth->execute()|| $log->log_and_die($sth->errstr);;
					 
while (my ( $gene, $live, $name ) = $sth->fetchrow_array){
  $server_genes{$gene}->{name} = $name;
  $server_genes{$gene}->{status} = $live;
}
}

#How much work is to be done
my $genecount = scalar keys %ace_genes;
$log->write_to("Checking $genecount genes\n");

# get a complete list of gene names from both Ace and server

my $err_count = 0;
foreach my $gene_id (sort keys %ace_genes) {
  $err_count += &check_gene($gene_id);
}

#finish up
$log->write_to("There were $err_count errors found\n");

$log->mail();
exit(0);


#########################
# Check_gene subroutine #
#########################

sub check_gene {
  my $gene = shift;

  if ($server_genes{$gene} and not $ace_genes{$gene}) {
    $log->error("ERROR: $gene missing from acedb\n");
    return 1;
  } elsif ($ace_genes{$gene} and not $server_genes{$gene}) {
    $log->error("ERROR: $gene missing from server\n");
    return 1;
  } else {
    my $err = 0;

    if (not $server_genes{$gene}->{name}) {
      $log->error("ERROR: no name for $gene in nameserver\n");
      $err = 1;
    } elsif (not $ace_genes{$gene}->{name}) {
      $log->error("ERROR: no name for $gene in acedb\n");
      $err = 1;
    } elsif ($server_genes{$gene}->{name} ne $ace_genes{$gene}->{name}) {
      $log->error(sprintf("ERROR: gene $gene has different names in Ace (%s) and Server (%s)\n", 
                          $ace_genes{$gene}->{name}, 
                          $server_genes{$gene}->{name}));
      $err = 1;      
    }

    if ($ace_genes{$gene}->{status} != $server_genes{$gene}->{status}){
      $log->error("ERROR: $gene live or dead ? ace".$ace_genes{$gene}->{status}." ns".$server_genes{$gene}->{status}."\n");
      $err = 1;
    }

    return $err;
  }
}
