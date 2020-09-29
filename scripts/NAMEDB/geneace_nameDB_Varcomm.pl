#!/software/bin/perl -w

use lib $ENV{'CVS_DIR'};
#use lib "$ENV{CVS_DIR}/NAMEDB/lib";

use strict;
use NameDB_handler;
use Wormbase;
use Log_files;
use Ace;
use Carp;
use Getopt::Long;
use Storable;

my ($help, $debug, $test, $store, $database, $outfile, $old, $dump);

GetOptions (	'help'       => \$help,
            	'debug=s'    => \$debug,
	    	'test'       => \$test,
  		'store:s'    => \$store,
                'database=s' => \$database,
		'outfile=s'  => \$outfile,
                'old'        => \$old,
                'dumpfile=s' => \$dump
           );

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

if (!defined $outfile) {
    $outfile = "VariationNSdiscrepancies.ace";
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);
my $errorcount = "0";

open (OUT, ">$outfile");

#connect to Geneace and read in data
my $acedb = ($database and -d $database) ? $database : $wormbase->database('geneace');

if ($acedb =~ (/\/$/)) {
    print "Removing final / from database path\n";
    $acedb =~ s/\/$//;
}

my $def = "$acedb/wquery/SCRIPT:geneace_nameDB_Varcomm.def";
my $TABLE = $wormbase->table_maker_query($acedb, $def);

my %ace_vars;
while( <$TABLE> ){
    #WBVar name status
	next if (/>/ or /\/\// );
	next unless (/\w/);
	s/\"//g;  # remove "
	my($var, $pub, $status) = split(/\s/);
	#Treat Suppressed genes the same as Dead as the nameserver has no concept of this 3rd status
	if ($status =~ /Supressed/){$status = "Dead";}
	$log->error("$var has no status:$_\n") unless $status;
	unless ($status) {
	    $errorcount++;
	    print "error1 $var\n";
	}
	$ace_vars{"$var"}->{'pub'} = $pub if $pub;
	$ace_vars{"$var"}->{'sts'} = ($status eq 'Live' ? 1 : 0);
}
           
#connect to name server and set domain to 'Gene'
my %server_vars;	
if ($old) {
my $DB      = 'nameserver_live;web-wwwdb-core-02;3449';
my $PASS    = 'wormpub';
my $USER    = 'wormpub';
my $DOMAIN  = 'Variation';
my $db = NameDB_handler->new($DB,$USER,$PASS);

# get nameserver data
my $query = 'SELECT pi.object_public_id, pi.object_live, si.object_name FROM   primary_identifier pi,secondary_identifier si WHERE  pi.object_id = si.object_id AND pi.domain_id = 3';


#SELECT pi.object_public_id, si.object_name FROM   primary_identifier pi,secondary_identifier si WHERE  pi.object_id = si.object_id AND pi.domain_id = 3'; 

#SELECT pi.object_public_id, si.name_type_id, si.object_name
#             FROM   primary_identifier pi,secondary_identifier si
#             WHERE  pi.object_id = si.object_id 
#             AND   (si.name_type_id = 3 OR si.name_type_id = 1) 
#             AND    pi.domain_id = 3';
# results
#| WBVar00000487    | 1 | bn85        |
#| WBVar00000488    | 1 | bn86        |
#| WBVar00000489    | 1 | bn87        |
#| WBVar00000490    | 1 | bn88        |

    my $sth = $db->dbh->prepare($query);
    $sth->execute();					 
    while (my ( $var, $live, $name ) = $sth->fetchrow_array){
        $server_vars{$var}->{'sts'} = $live;
        $server_vars{$var}->{'pub'} = $name;
    }
}
else {
    print "Running on the new Datomic dump file\n";
    #WBVar00094646,live,or6
    #WBVar00094643,live,op509
    my ($dump_fh);
    open($dump_fh, "<$dump") or $log->log_and_die("Can't open $dump\n");
    while ( my $line = <$dump_fh> ) {
        my $var;
        my $name;
        my $status;
        chomp $line;
        my @f = split ",", $line;
        unless ($f[0] =~ /WBVar/) {next;}
        $var = $f[0];  
        if ($f[2] =~ /\S+/) {
            $name = $f[2];
            $server_vars{$var}->{'pub'} = $name;
        }
        if ($f[1] =~ /\S+/) {
            $status = $f[1];
            if ($status =~ /live/){
            $server_vars{$var}->{'sts'} = '1';
            }
            if ($status =~ /dead/){
            $server_vars{$var}->{'sts'} = '0';
            }
        }
    }
}

# loop through the genes in the nameserver and check them
# delete from the acedb list as it goes as there's no need to check twice
foreach my $var (keys %server_vars) {
	&check_var($var);
	delete($ace_vars{$var}) if $ace_vars{$var};
}

# any genes left in the acedb list are absent from the nameserver
foreach (keys %ace_vars ){
    if (($_ =~ /ENS/) || ($_ =~ /YMR133W/)) {
	print "Skipping ENSEMBL ID:$_\n" if $debug;
    }
    else {
	# We don't really care about WBVariations not in the nameserver so changing to just flag non WBVar{8} IDs
	$log->error("ERROR: $_ absent from nameserver\n") unless ($_ =~ /WBVar\d{8}/);
	$errorcount++ unless ($_ =~ /WBVar\d{8}/);
    }
}
$log->write_to("INFO: $errorcount errors found\n") if ($log->report_errors > 0);
$log->write_to("No errors found\n") if ($log->report_errors == 0);
$log->mail;
exit(0);

sub check_var {
    my $var = shift;
    #check presence in both sets

#is it in geneace
    if( $ace_vars{"$var"} ){
#is it also in the nameserver
	if( $server_vars{"$var"} ){
	    # check Live 
	    if ($ace_vars{"$var"}->{'sts'} != $server_vars{"$var"}->{'sts'}){
		my $nsstatus;
		my $acestatus;
		if ($server_vars{"$var"}->{'sts'} eq "0") {$nsstatus = "Dead";}
		elsif ($server_vars{"$var"}->{'sts'} eq "1") {$nsstatus = "Live";}
		if ($ace_vars{"$var"}->{'sts'} eq "0") {$acestatus = "Dead";}
		elsif ($ace_vars{"$var"}->{'sts'} eq "1") {$acestatus = "Live";}
		$log->error("ERROR: $var live or dead ? ns:$nsstatus g:$acestatus\n");
		$errorcount++;
		print "error 3 $var\n";
		return;
	    }
	    return if ($ace_vars{"$var"}->{'sts'} == 0) ; #dont check dead genes.
	    # for each name type eg cgc, seq
	    
# check name
#confirm both have them
	    if ($ace_vars{"$var"}->{'pub'} or $server_vars{"$var"}->{'pub'}) {
		#determine if either doesn't
		unless ($ace_vars{"$var"}->{'pub'}) {
                    $log->error("ERROR: Geneace - $var missing 'pub' name ".$server_vars{"$var"}->{'pub'}."\n");
                    $errorcount++;
		    print "error 4 $var\n";
		    return;
		}			
		unless ($server_vars{"$var"}->{'pub'}) {
                    $log->error("ERROR: Server - $var missing 'pub' name ".$ace_vars{"$var"}->{'pub'}."\n");
                    $errorcount++;
		    print "error 5 $var\n";
		    return;
		}
		# confirm they are the same if both have name type	
		if ($ace_vars{"$var"}->{'pub'} ne $server_vars{"$var"}->{'pub'}) {
		    $log->error("ERROR: $var different 'pub' names for $var G: ".$ace_vars{"$var"}->{'pub'}." N: ".$server_vars{"$var"}->{'pub'}."\n");
		    print OUT "\nVariation $var\nPublic_name ".$server_vars{"$var"}->{'pub'}."\n\n"; 
		    $errorcount++;
		    print "error 6 $var\n";
		    return;
		}
	    }
	    #print STDERR "$var ok\n" if $debug;
	    return;
	}
	else {
	    # nameserver query will not return genes imported without a valid name.
	    $log->error("ERROR: $var missing from server\n");
	    $errorcount++;
	    print "error 7 $var\n";
	    return;
	}	
    }
    else {
	$log->error("ERROR: $var missing from acedb\n");
	if ($server_vars{"$var"}->{'pub'}) {
	    print OUT "Variation $var\nRemark \"Variation stub created from the Nameserver report\"\nPublic_name ".$server_vars{"$var"}->{'pub'}."\n\n";
	}
	else {
	    print OUT "\n\/\/Variation $var has no name data, skipping\n";	
	}
	$errorcount++;
	print "error 8 $var\n";
	return;
    }
}
