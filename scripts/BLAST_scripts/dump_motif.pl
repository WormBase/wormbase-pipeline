#!/usr/local/bin/perl5.6.1 -w
#
# Originally written by Marc Sohrmann (ms2@sanger.ac.uk)
#
# Dumps protein motifs from ensembl mysql (protein) database to an ace file
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2003-06-06 09:45:16 $


use strict;
use DBI;
use Getopt::Long;

my ($debug, $WPver, $database);

GetOptions("debug" => \$debug,
	   "verison" => \$WPver,
	   "database" => \$database
	  );

# define the names of the methods to be dumped
my @methods = qw(ncoils seg signalp tmhmm hmmpfam);

# mysql database parameters
my $dbhost = "ecs1f";
my $dbuser = "wormro";
my $dbname = "wormprot";
$dbname = $database if $database;
my $dbpass = "";

# to get the current time...
sub now {
    return sprintf ("%04d-%02d-%02d %02d:%02d:%02d",
                     sub {($_[5]+1900, $_[4]+1, $_[3], $_[2], $_[1], $_[0])}->(localtime));
}

# create output files
my $dump_dir = "/acari/work2a/wormpipe/dumps";
open(ACE,">$dump_dir/".$dbname."_motif_info.ace") || die "cannot create ace file";

open(LOG,">$dump_dir/".$dbname."_motif_info.log") || die "cannot create log file";

# make the LOG filehandle line-buffered
my $old_fh = select(LOG);
$| = 1;
select($old_fh);

$old_fh = select(ACE);
$| = 1;
select($old_fh);


print LOG "DUMPing protein motif data from ".$dbname." to ace [".&now."]\n";
print LOG "---------------------------------------------------------------\n\n";

# connect to the mysql database
print LOG "connect to the mysql database $dbname on $dbhost as $dbuser [".&now."]\n\n";
my $dbh = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, $dbpass, {RaiseError => 1})
    || die "cannot connect to db, $DBI::errstr";

# get the mapping of method 2 analysis id
my %method2analysis;
print LOG "get mapping of method to analysis id [".&now."]:\n";
my $sth = $dbh->prepare ( q{ SELECT analysisId
                               FROM analysisprocess
                              WHERE program = ?
                           } );

foreach my $method (@methods) {
    $sth->execute ($method);
    (my $anal) = $sth->fetchrow_array;
    $method2analysis{$method} = $anal;
    print LOG "$method  $anal\n";
}

# prepare the sql querie
my $sth_f = $dbh->prepare ( q{ SELECT proteinId, start, end, hid, hstart, hend, score
                                 FROM protein_feature
                                WHERE analysis = ?
                             } );

# get the motifs
my %motifs;
my %pfams;
foreach my $method (@methods) {
    print LOG "processing $method\n";
    $sth_f->execute ($method2analysis{$method});
    my $ref = $sth_f->fetchall_arrayref;
    foreach my $aref (@$ref) {
        my ($prot, $start, $end, $hid, $hstart, $hend, $score) = @$aref;
        if ($method eq "hmmpfam") {
            my $line = "Motif_homol \"PFAM:$hid\" \"pfam\" $score $start $end $hstart $hend";
            push (@{$motifs{$prot}} , $line);
	}
        else {
            my $line = "Feature \"$method\" $start $end $score";
            push (@{$motifs{$prot}} , $line);
	}
    }
}

# print ace file
my $prefix = "WP";
if( "$database" eq "worm_brigprot") {
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

print LOG "\nEnd of Motif dump\n";
print "\nEnd of Motif dump\n";
close LOG;
exit(0);
