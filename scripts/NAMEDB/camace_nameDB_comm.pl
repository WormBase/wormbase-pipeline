#!/nfs/disk100/wormpub/bin/perl -w

use lib '/nfs/WWWdev/SANGER_docs/cgi-bin/Projects/C_elegans/lib';
use lib $ENV{'CVS_DIR'};

use strict;
use NameDB_handler;
use Wormbase;
use Log_files;
use Ace;
use Carp;
use Getopt::Long;
use Storable;

my ($help, $debug, $test, $store, $database, $def);

GetOptions (	"help"       => \$help,
            	"debug=s"    => \$debug,
		"test"       => \$test,
		"store:s"    => \$store,
		"database:s" => \$database,
		"def:s"      => \$def
           );

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}


# establish log file.
my $log = Log_files->make_build_log($wormbase);

#connect to database and read in data
my $acedb = $database or $wormbase->database('camace');
$def = "$acedb/wquery/SCRIPT:camace_nameDB_comm.def" unless $def;
my $TABLE = $wormbase->table_maker_query($acedb, $def);

my %ace_genes;
while( <$TABLE> ){
	next if (/>/ or /\/\// );
	s/\"//g;  # remove "
	my($gene, $cds, $transcript, $pseudo) = split(/\s/);
	next unless ($cds or $transcript or $pseudo);
	#$ace_genes{"$gene"}->{'cds'}->{"$cds"}        = 1 if $cds;
	$ace_genes{"$gene"}->{'tra'}->{"$transcript"} = 1 if $transcript;
	$ace_genes{"$gene"}->{'pse'}->{"$pseudo"}     = 1 if $pseudo;
	$ace_genes{"$gene"}->{'sts'} = ($cds or $transcript or $pseudo) ?  1 : 0; #live if it has one these nametypes
}
           
#connect to name server and set domain to 'Gene'
my $DB    	= 'wbgene_id;mcs2a';
my $PASS 	= "wormpub";
my $USER 	= "wormpub";
my $DOMAIN  = 'Gene';
my $db = NameDB_handler->new($DB,$USER,$PASS);

# get nameserver data
my $query = "SELECT primary_identifier.object_public_id, 
                    primary_identifier.object_live, 
                    secondary_identifier.object_name,
                    secondary_identifier.name_type_id 
             FROM primary_identifier,secondary_identifier 
             WHERE  primary_identifier.object_id = secondary_identifier.object_id 
             AND   (secondary_identifier.name_type_id = 3) 
             ORDER BY object_public_id";

# results
#| object_public_id |  live       | name_type_id | name
#| WBGene00044331   |           1 |            2 | T20B3.15      |  no longer used
#| WBGene00044331   |           1 |            3 | T20B3.15      |
				 
my $sth = $db->dbh->prepare($query);
$sth->execute();

my %name_types = (
		  '2' => 'cds',
		  '3' => 'seq'
		  );

my %server_genes;					 
while (my ( $gene, $live, $name, $name_type ) = $sth->fetchrow_array){
	$server_genes{$gene}->{ $name_types{$name_type} }->{"$name"} = 1;
	$server_genes{$gene}->{ 'sts' } = $live;
}


foreach my $gene (keys %ace_genes) {
	&check_gene($gene);
	delete($ace_genes{$gene}) if $ace_genes{$gene};
}

$log->mail;
exit(0);


sub check_gene {
    my $gene = shift;
    #check presence in both sets
    if( $ace_genes{"$gene"} ){
	if( $server_genes{"$gene"} ){
	    # check Live 
	    if ($ace_genes{"$gene"}->{'sts'} != $server_genes{"$gene"}->{'sts'}){
		$log->write_to("ERROR: $gene live or dead ? ace".$ace_genes{"$gene"}->{'sts'}." ns".$server_genes{"$gene"}->{'sts'}."\n");
		return;
	    }
	    # for each name type eg cgc, seq
	    foreach my $type ('seq') {
		foreach my $cds (keys %{$ace_genes{"$gene"}->{"$type"} } ) {
		    unless ( $server_genes{"$gene"}->{"$type"}->{"$cds"} ) {
			$log->write_to("ERROR: $gene - ${type}_name missing in server\n");
			next;
		    }
		}
	    }
	    #server will have seq names for all CDS, acedb wont. 
#	    foreach my $seq (keys %{$server_genes{"$gene"}->{'seq'}}) {
	#	next if ($ace_genes{"$gene"}->{'cds'} );
		#unless ($ace_genes{"$gene"}->{'seq'}->{"$seq"} ) {
	 #   $log->write_to("ERROR: $gene - only Sequence_name in server\n");
#		}
#	    }
	    
	    #print STDERR "$gene ok\n" if $debug;
	    return;
	}
	else {
	    $log->write_to("ERROR: $gene missing from server\n");
	    return;
	}	
    }
    else {
	$log->write_to("ERROR: $gene missing from acedb\n");
	return;
    }
}
