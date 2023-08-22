#!/usr/bin/env perl
#
# Originally written by Marc Sohrmann (ms2@sanger.ac.uk)
#
# Dumps protein motifs from ensembl mysql (protein) database to an ace file
#
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2014-08-18 08:58:32 $

use lib $ENV{'CVS_DIR'};

use strict;
use DBI;
use Getopt::Long;
use Wormbase;
use Storable;
use Log_files;

my ($WPver, @methods);
my ($store, $test, $debug,$dump_dir,$dbname);

GetOptions(
	   "database:s" => \$dbname,
	   "methods=s"   => \@methods,
	   "store:s"   => \$store,
	   "test"      => \$test,
	   "debug:s"   => \$debug,
	   "dumpdir=s" => \$dump_dir,
	   "dbname=s"  => \$dbname,
	  );

$dump_dir ||= $ENV{'PIPELINE'}.'/dumps';

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

my $log = Log_files->make_build_log($wormbase);

# define the names of the methods to be dumped
@methods = qw(hmmpanther pfscan pirsf prints scanprosite smart ncbifam ncoils seg signalp tmhmm pfam superfamily mobidblite hamap) unless @methods;

$log->write_to("Dumping methods".@methods."\n");

# mysql database parameters
my $dbhost = $ENV{'WORM_DBHOST'};;
my $dbuser = "wormro";
my $dbport = $ENV{'WORM_DBPORT'};
$dbname ||= "worm_ensembl_elegans";
print "Dumping motifs from $dbname\n";
my $dbpass = "";

# to get the current time...
sub now {
    return sprintf ("%04d-%02d-%02d %02d:%02d:%02d",
                     sub {($_[5]+1900, $_[4]+1, $_[3], $_[2], $_[1], $_[0])}->(localtime));
}

# create output files
open(ACE,">$dump_dir/".$dbname."_motif_info.ace") or $log->log_and_die("cannot create ace file:$!\n");
open(LOG,">$dump_dir/".$dbname."_motif_info.log") or $log->log_and_die( "cannot create log file:$!\n");

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
my $dbh = DBI -> connect("DBI:mysql:$dbname:$dbhost:$dbport", $dbuser, $dbpass, {RaiseError => 1})
    || $log->log_and_die("cannot connect to db, $DBI::errstr\n");

# get the mapping of method 2 analysis id
my %method2analysis;
print LOG "get mapping of method to analysis id [".&now."]:\n";
my $sth = $dbh->prepare ( q{ SELECT analysis_id
                               FROM analysis
                              WHERE logic_name = ?
                           } );

foreach my $meth (@methods) {
    $sth->execute ($meth);
    (my $anal) = $sth->fetchrow_array;
    $method2analysis{$meth} = $anal;
    $log->write_to("method: $meth => analyis_id: $anal\n");
}


# prepare the sql query
my $sth_f;
if (-e "/software/worm") { # running on Sanger
  $sth_f = $dbh->prepare ( q{ SELECT stable_id, seq_start, seq_end, hit_name, hit_start, hit_end, score
                                 FROM protein_feature,translation_stable_id
                                WHERE analysis_id = ? AND translation_stable_id.translation_id = protein_feature.translation_id
                             } );
} else {
  $sth_f = $dbh->prepare ( q{ SELECT stable_id, p.seq_start, p.seq_end, hit_name, hit_start, hit_end, score
                                 FROM protein_feature p, translation t
                                WHERE analysis_id = ? AND t.translation_id = p.translation_id
                             } );
}  


# get the motifs
my %motifs;
my %pfams;
my %panther;
my %cds2wormpep;
my %feature2remark; # store remark lines, so we can populate them in the Motif objects

$wormbase->FetchData('cds2wormpep',\%cds2wormpep);

foreach my $meth (@methods) {
  $log->write_to("processing $meth\n");
  $sth_f->execute ($method2analysis{$meth});
  my $ref = $sth_f->fetchall_arrayref;
  foreach my $aref (@$ref) {
    my ($_prot, $start, $end, $hid, $hstart, $hend, $score) = @$aref;
    my $prot=($cds2wormpep{$_prot}||$_prot);

    $score||=0; # for NULL scores in some analysis, as else it will break Ace

    my $line;
    if ($meth eq "pfam") {
       if( $hid =~ /(\w+)\.\d+/ ) {
	$hid = $1;
       }
       $line = "Motif_homol \"PFAM:$hid\" \"pfam\" $score $start $end $hstart $hend";
    }elsif($meth=~/hmmpanther|pfscan|pirsf|prints|scanprosite|smart|ncbifam|superfamily/){
       my $prefix = uc($meth);
       my $motif = "$prefix:$hid"; # shorthand for the acedb id
       $line = "Motif_homol \"$motif\" \"$meth\" $score $start $end $hstart $hend";

       $feature2remark{"$motif"}="Remark \"$meth motif $hid\" From_analysis $meth\nDatabase $meth ${meth}_ID \"$hid\"";

       $panther{"$motif"}=$hid if $meth eq 'hmmpanther';
    }else { # tmhmm seg signalp ncoils
       $line = "Feature \"$meth\" $start $end $score";

       $feature2remark{"$meth"}="Remark \"$meth motif\" From_analysis $meth";
    }
    push (@{$motifs{$prot}} , $line) if $line;
  }
}


foreach my $prot (sort {$a cmp $b} keys %motifs) {
    print ACE "\n";
    # cds2wormpep conversion
    print ACE "Protein : \"$prot\"\n";
    foreach my $line (@{$motifs{$prot}}) {
        print ACE "$line\n";
    }
}

print ACE "\n";

#################################
while(my($k,$v)= each %panther){
   print ACE "Motif : \"$k\"\n";
   print ACE "Database hmmpanther PantherID \"$v\"\n\n"
}

while(my($k,$v)=each %feature2remark){
   print ACE "Motif : \"$k\"\n";
   print ACE "$v\n\n";
}
#################################

$sth->finish;
$sth_f->finish;
$dbh->disconnect;

close ACE;

$log->write_to("\nEnd of Motif dump\n");

$log->mail;
exit(0);

