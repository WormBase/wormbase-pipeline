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

my ($help, $debug, $test, $store, $database, $dump, $old);

GetOptions (	'help'       => \$help,     # help documentation
            	'debug=s'    => \$debug,    # only email specified user
	    	'test'       => \$test,     # use the test server credentials
  		'store:s'    => \$store,    # wormbase storable object
                'database=s' => \$database, # specify an alternative database to the standard geneace
                'dumpfile=s' => \$dump,     # file dumped from the datomic Name service check the confluence guide on how to produce this
                'old'        => \$old,      # option to use the old MySQL nameserver
           );

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

log->log_and_die("you need to specify -dumpfile and point it to a datomic nameservice dump file\n") unless ($dump);

# establish log file.
my $log = Log_files->make_build_log($wormbase);
my $errorcount = "0";
#connect to Geneace and read in data
my $acedb = ($database and -d $database) ? $database : $wormbase->database('geneace');
my $def = "$acedb/wquery/SCRIPT:geneace_nameDB_comm.def";
my $TABLE = $wormbase->table_maker_query($acedb, $def);

my %ace_genes;
while( <$TABLE> ){
	next if (/>/ or /\/\// );
	next unless (/\w/);
	s/\"//g;  # remove "
	my($gene, $cgc, $seq, $status) = split(/\s/);
	$log->error("$gene has no status:$_\n") unless $status;
	$errorcount++ unless $status;
	$ace_genes{"$gene"}->{'cgc'} = $cgc if $cgc;
	$ace_genes{"$gene"}->{'seq'} = $seq if $seq;
	$ace_genes{"$gene"}->{'sts'} = ($status eq 'Live' ? 1 : 0);
}
           


my %name_types = ( '1' => 'cgc', '3' => 'seq');
my %server_genes;					 

if ($old) {

#connect to name server and set domain to 'Gene'
my $DB      = 'nameserver_live;web-wwwdb-core-02;3449';
my $PASS    = 'wormpub';
my $USER    = 'wormpub';
my $DOMAIN  = 'Gene';
my $db = NameDB_handler->new($DB,$USER,$PASS);

# get nameserver data
my $query = 'SELECT pi.object_public_id, si.name_type_id, si.object_name
             FROM   primary_identifier pi,secondary_identifier si
             WHERE  pi.object_id = si.object_id 
             AND   (si.name_type_id = 3 OR si.name_type_id = 1) 
             AND    pi.domain_id = 1';
# results
#| object_public_id | name_type_id | name
#| WBGene00044331   |            1 | clec-25       |
#| WBGene00044331   |            3 | T20B3.15      |
                                 
my $sth = $db->dbh->prepare($query);
$sth->execute();

while (my ( $gene, $name_type, $name ) = $sth->fetchrow_array){
        $server_genes{$gene}->{ $name_types{$name_type} } = $name;
}

#This is to get dead and briggsae genes that dont have names.
$query = 'SELECT object_public_id, object_live
          FROM   primary_identifier 
          WHERE domain_id = 1';

$sth = $db->dbh->prepare($query);
$sth->execute();

while (my ( $gene, $live ) = $sth->fetchrow_array){
        $server_genes{$gene}->{'sts'} = $live;
}


}

# Datomic world section
    unless ($old) {
	$log->write_to("Using $dump file for checking\n");
        my ( $gene, $name_type, $name);
        my ($dump_fh);
        open($dump_fh, "<$dump") or $log->log_and_die("Can't open $dump\n"); 
        while ( my $line = <$dump_fh> ) {
            chomp $line;
            my @f = split ",", $line;
            unless ($f[0] =~ /WBGene/) {next;}
            $gene = $f[0];            
            #Store the cgc and seq names
            if ($f[1] =~ /\S+/) {
                $server_genes{$gene}->{'cgc'} = $f[1]; 
            }
            if ($f[2] =~ /\S+/) {
                $server_genes{$gene}->{'seq'} = $f[2]; 
            }
            #Status 1 or 0
            if ($f[3] =~ /live/) {
                $server_genes{$gene}->{'sts'} = "1";
            }
            else {
                $server_genes{$gene}->{'sts'} = "0";
            }
        }		 
    }
    
# loop through the genes in the nameserver and check them
# delete from the acedb list as it goes as there's no need to check twice
foreach my $gene (keys %server_genes) {
	&check_gene($gene);
	delete($ace_genes{$gene}) if $ace_genes{$gene};
}

# any genes left in the acedb list are absent from the nameserver
foreach (keys %ace_genes ){
    if (($_ =~ /ENS/) || ($_ =~ /YMR133W/)) {
	print "Skipping ENSEMBL ID:$_\n" if $debug;
    }
    else {
	$log->error("ERROR: $_ absent from nameserver\n");
	$errorcount++;
    }
}
$log->write_to("INFO: $errorcount errors found\n") if ($log->report_errors > 0);
$log->write_to("No errors found\n") if ($log->report_errors == 0);
$log->mail;
exit(0);

sub check_gene {
	my $gene = shift;
	#check presence in both sets
	if( $ace_genes{"$gene"} ){
	  if( $server_genes{"$gene"} ){
	      # check Live 
	      if ($ace_genes{"$gene"}->{'sts'} != $server_genes{"$gene"}->{'sts'}){
		$log->error("ERROR: $gene live or dead ? ns".$server_genes{"$gene"}->{'sts'}." g".$ace_genes{"$gene"}->{'sts'}."\n");
		$errorcount++;
		return;
	      }
	      return if ($ace_genes{"$gene"}->{'sts'} == 0) ; #dont check dead genes.
	      # for each name type eg cgc, seq
	      my ($i,$type);	
	      while( ($i,$type) = each(%name_types) ){
		#confirm both have them
		if ($ace_genes{"$gene"}->{$type} or $server_genes{"$gene"}->{$type}) {
                  #determine if either doesn't
                  unless ($ace_genes{"$gene"}->{$type}) {
                    $log->error("ERROR: Geneace - $gene missing $type name ".$server_genes{"$gene"}->{$type}."\n");
                    $errorcount++;
		    next;
                  }			
                  unless ($server_genes{"$gene"}->{$type}) {
		      unless (($gene =~ /WBGene00006756/) ||($gene =~ /WBGene00044083/) || ($gene =~ /WBGene00000481/) || ($gene =~ /WBGene00001147/))  {    # B0564.1 Y105E8A.7 ZC416.8)
			  $log->error("ERROR: Server - $gene missing $type name ".$ace_genes{"$gene"}->{$type}."\n");
			  $errorcount++;
		      }
		      next;
                  }
                  # confirm they are the same if both have name type	
                  if ($ace_genes{"$gene"}->{$type} ne $server_genes{"$gene"}->{$type}) {
                    $log->error("ERROR: $gene different $type names for $gene G:".$ace_genes{"$gene"}->{$type}." N:".$server_genes{"$gene"}->{$type}."\n");
                    $errorcount++;
		    next;
                  }
		}
	      }
	      #print STDERR "$gene ok\n" if $debug;
	      return;
	}
	else {
	  # nameserver query will not return genes imported without a valid name.
	  $log->error("ERROR: $gene missing from server\n");
	  $errorcount++;
	  return;
	}	
      }
      else {
	$log->error("ERROR: $gene missing from acedb\n");
	$errorcount++;
	return;
      }
}
