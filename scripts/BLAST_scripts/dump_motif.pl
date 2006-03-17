#!/usr/local/ensembl/bin/perl -w
#
# Originally written by Marc Sohrmann (ms2@sanger.ac.uk)
#
# Dumps protein motifs from ensembl mysql (protein) database to an ace file
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2006-03-17 11:43:18 $

use lib $ENV{'CVS_DIR'};

use strict;
use DBI;
use Getopt::Long;
use Wormbase;
use Storable;
use Log_files;

my ($WPver, $database, $mysql, $method);
my ($store, $test, $debug);

GetOptions(
	   "database:s" => \$database,
	   "mysql"      => \$mysql,
	   "method=s"   => \$method,
	   "store:s"   => \$store,
	   "test"      => \$test,
	   "debug:s"   => \$debug
	  );

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

my $log = Log_files->make_build_log($wormbase);
my $dump_dir = '/acari/work2a/wormpipe/dumps';

# define the names of the methods to be dumped
my @methods;
if ($method ) {
  push(@methods,$method)
}else{
  @methods= qw(ncoils seg signalp tmhmm hmmpfam);
}
$log->write_to("Dumping methods".@methods."\n");

# mysql database parameters
my $dbhost = "ecs1f";
my $dbuser = "wormro";
my $dbname = "worm_pep";
$dbname = $database if $database;
print "Dumping motifs from $dbname\n";
my $dbpass = "";

# to get the current time...
sub now {
    return sprintf ("%04d-%02d-%02d %02d:%02d:%02d",
                     sub {($_[5]+1900, $_[4]+1, $_[3], $_[2], $_[1], $_[0])}->(localtime));
}

# create output files
open(ACE,">$dump_dir/".$dbname."_motif_info.ace") || die "cannot create ace file:$!\n";
open(LOG,">$dump_dir/".$dbname."_motif_info.log") || die "cannot create log file:$!\n";

# make the LOG filehandle line-buffered
my $old_fh = select(LOG);
$| = 1;
select($old_fh);

$old_fh = select(ACE);
$| = 1;
select($old_fh);


$log->write_to("DUMPing protein motif data from ".$dbname." to ace\n---------------------------------------------------------------\n\n");

# connect to the mysql database
print LOG "connect to the mysql database $dbname on $dbhost as $dbuser [".&now."]\n\n";
my $dbh = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
    || $log->log_and_die("cannot connect to db, $DBI::errstr\n");

# get the mapping of method 2 analysis id
my %method2analysis;
print LOG "get mapping of method to analysis id [".&now."]:\n";
my $sth = $dbh->prepare ( q{ SELECT analysis_id
                               FROM analysis
                              WHERE program = ?
                           } );

foreach my $meth (@methods) {
    $sth->execute ($meth);
    (my $anal) = $sth->fetchrow_array;
    $method2analysis{$meth} = $anal;
    $log->write_to("$meth  $anal\n");
}

# prepare the sql querie
my $sth_f = $dbh->prepare ( q{ SELECT protein_id, seq_start, seq_end, hit_id, hit_start, hit_end, score
                                 FROM protein_feature
                                WHERE analysis_id = ?
                             } );

# get the motifs
my %motifs;
my %pfams;
foreach my $meth (@methods) {
  $log->write_to("processing $meth\n");
  $sth_f->execute ($method2analysis{$meth});
  my $ref = $sth_f->fetchall_arrayref;
  foreach my $aref (@$ref) {
    my ($prot, $start, $end, $hid, $hstart, $hend, $score) = @$aref;
    my $line;
    if ($meth eq "hmmpfam") {
      if( $hid =~ /(\w+)\.\d+/ ) {
	$hid = $1;
      }
       $line = "Motif_homol \"PFAM:$hid\" \"pfam\" $score $start $end $hstart $hend";
      push (@{$motifs{$prot}} , $line);
    }
    else {
      $line = "Feature \"$meth\" $start $end $score";
      push (@{$motifs{$prot}} , $line);
    }
  }
}

# print ace file
my $prefix = "WP";
if( "$database" eq "worm_brigpep") {
  $prefix = "BP";
}
foreach my $prot (sort {$a cmp $b} keys %motifs) {
    print ACE "\n";
    print ACE "Protein : \"$prefix:$prot\"\n";
    foreach my $line (@{$motifs{$prot}}) {
        print ACE "$line\n";
    }
}

    
$sth->finish;
$sth_f->finish;
$dbh->disconnect;

close ACE;

$log->write_to("\nEnd of Motif dump\n");
print "\nEnd of Motif dump\n";

$log->mail;
exit(0);
