#!/usr/local/bin/perl

# Marc Sohrmann (ms2@sanger.ac.uk)

# takes as input a swissprot or trembl .fasta file,
# and deletes all worm, fly, human and yesast entries

use strict;
use Getopt::Std;
use DB_File;
use vars qw($opt_s $opt_t $opt_d);

getopts ("std:");

my %ORG;
my %DES;

if ($opt_s && $opt_t) {
    die "usage -s for swissprot,-t for trembl, -d directory name where files are < file with ids\n";
}
elsif ($opt_s) {
    unless (-s "$opt_d/swissprot2org") {
        die "swissprot2org not found or empty";
    }
    dbmopen %ORG, "$opt_d/swissprot2org", 0666 or die "cannot open swissprot2org DBM file";
    unless (-s "$opt_d/swissprot2des") {
        die "swissprot2des not found or empty";
    }
    dbmopen %DES, "$opt_d/swissprot2des", 0666 or die "cannot open swissprot2des DBM file";
}
elsif ($opt_t) {
    unless (-s "$opt_d/trembl2org") {
        die "trembl2org not found or empty";
    }
    dbmopen %ORG, "$opt_d/trembl2org", 0666 or die "cannot open trembl2org DBM file";
    unless (-s "$opt_d/trembl2des") {
        die "trembl2des not found or empty";
    }
    dbmopen %DES, "$opt_d/trembl2des", 0666 or die "cannot open trembl2des DBM file";
}
else {
    die "usage -s for swissprot,-t for trembl, -d directory name where files are < file with ids\n";
}

print "\n";
while (<>) {
    chomp;
    /^(\S+)/;
    my $db;
    my $fetch_db;
    my $prefix;
    my $id = $1;
    if ($opt_s) {
        $db = "SwissProt";
        $fetch_db = "slimswissprot";
        $prefix = "SW";
    }
    elsif ($opt_t) {
        $db = "TrEMBL";
        $fetch_db = "slimtrembl";
        $prefix = "TR";
    }
    else {
        die "wrong options";
    }
    my $entry = `/nfs/acari/wormpipe/Pipeline/fetch.pl -g $opt_d/$fetch_db.gsi -i $id`;

    chomp $entry;
    my @ary = split (/\n/, $entry);
    my $header = shift @ary;
    $header =~ /^\S+\s+(\S+)/;
    my $accession = $1;

    my $a = `pfetch $id`;

    $a =~ /^>\S+\s+(\S+)/;
    $accession = $1;

    my $seq = join ("", @ary);
    if ($seq) {
        print "Protein : \"$prefix:$id\"\n";
        print "Peptide \"$prefix:$id\"\n";
        print "Title \"$DES{$id}\"\n";
        print "Species \"$ORG{$id}\"\n";
        print "Database \"$db\" \"$id\" \"$accession\"\n";
        print "\n";
        print "Peptide : \"$prefix:$id\"\n";
        print "$seq\n";
        print "\n";
    }
    else {
        die "couldn't fetch sequence for $id";
    }
}

dbmclose %ORG;
dbmclose %DES;



