#!/usr/local/bin/perl5.6.1 -w
#
# Originally written by Marc Sohrmann (ms2@sanger.ac.uk)
#
# Dumps protein motifs from ensembl mysql (protein) database to an ace file
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2002-11-05 13:53:39 $


use strict;
use DBI;
use Getopt::Std;

# define the names of the methods to be dumped
my @methods = qw(ncoils seg signalp tmhmm hmmpfam);

# mysql database parameters
my $dbhost = "ecs1f";
my $dbuser = "wormro";
my $dbname = "wormprot";
my $dbpass = "";

# to get the current time...
sub now {
    return sprintf ("%04d-%02d-%02d %02d:%02d:%02d",
                     sub {($_[5]+1900, $_[4]+1, $_[3], $_[2], $_[1], $_[0])}->(localtime));
}

# create output files
open (LOG, ">motifs_dump.log") || die "cannot create log file";
open (ACE, ">motifs_dump.ace") || die "cannot create ace file";

# make the LOG filehandle line-buffered
my $old_fh = select(LOG);
$| = 1;
select($old_fh);

$old_fh = select(ACE);
$| = 1;
select($old_fh);


print LOG "DUMPing protein motif data from mysql to ace [".&now."]\n";
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
           # my $line = "Pep_homol \"WP:prot\" \"pfam\" $score $hstart $hend $start $end";
           # push (@{$pfams{$hid}} , $line);
	}
        else {
            my $line = "Feature \"$method\" $start $end $score";
            push (@{$motifs{$prot}} , $line);
	}
    }
}

# print ace file
foreach my $prot (sort {$a cmp $b} keys %motifs) {
    print ACE "\n";
    print ACE "Protein : \"WP:$prot\"\n";
    foreach my $line (@{$motifs{$prot}}) {
        print ACE "$line\n";
    }
}

#foreach my $pfam (sort {$a cmp $b} keys %pfams) {
#    print ACE "\n";
#    print ACE "Motif : \"WP:$pfam\"\n";
#    foreach my $line (@{$motifs{$pfam}}) {
#        print ACE "$line\n";
#    }
#}
    
$sth->finish;
$sth_f->finish;
$dbh->disconnect;

print LOG "\nEnd of Motif dump\n";
print "\nEnd of Motif dump\n";
