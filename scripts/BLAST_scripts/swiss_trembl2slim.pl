#!/usr/local/bin/perl

# Marc Sohrmann (ms2@sanger.ac.uk)

use strict;
use Getopt::Std;
use DB_File;
use vars qw($opt_s $opt_t);

getopts ("st");

my $usage = "cat swissprot/trembl .fasta file | swiss_trembl2slim.pl\n";
$usage .= "-s for swissprot\n";
$usage .= "-t for trembl\n";


my %exclude;
$exclude{'Caenorhabditis elegans'} = 1;
$exclude{'Drosophila melanogaster'} = 1;
$exclude{'Saccharomyces cerevisiae'} = 1;
$exclude{'Homo sapiens'} = 1;

my %HASH;

if ($opt_s && $opt_t) {
    die "$usage";
}
elsif ($opt_s) {
    unless (-s "swissprot2org") {
        die "swiss2org not found or empty";
    }
    dbmopen %HASH, "swissprot2org", 0666 or die "cannot open DBM file";
}
elsif ($opt_t) {
    unless (-s "trembl2org") {
        die "trembl2org not found or empty";
    }
    dbmopen %HASH, "trembl2org", 0666 or die "cannot open DBM file";
}
else {
    die "$usage";
}

read_fasta (\*STDIN);
dbmclose %HASH;

sub read_fasta {
    local (*FILE) = @_;
    my ($id, $acc, $seq);
    while (<FILE>) {
        chomp;
        if (/^>(\S+)\s+(\S+)/) {
            my $new_id = $1;
	    my $new_acc = $2;
            if ($id) {
                $seq =~ tr/a-z/A-Z/;
                my $org;
	 	if (exists($HASH{$id})){
	            $org = $HASH{$id};
                    if (exists $exclude{$org}) {
                    	$id = $new_id; $acc = $new_acc; $seq = "" ;
                    	next;
		    }
                    my $count = 0;
                    $seq = reverse $seq;
                    print ">$id $acc ($org)";
                    while (my $base = chop $seq) {
                    	if ($count++ % 50 == 0) {
                        	print "\n";
                    	}
                    	print $base;
                    }
                    print "\n";
            	}
            }
	    $id = $new_id ;$acc = $new_acc; $seq = "" ;
	} 
        elsif (eof) {
            if ($id) {
                $seq .= $_ ;
                $seq =~ tr/a-z/A-Z/;
                my $org;
	 	if (exists($HASH{$id})){
	            $org = $HASH{$id};
                    if (exists $exclude{$org}) {
                    	next;
		    }
                    my $count = 0;
                    $seq = reverse $seq;
                    print ">$id $acc ($org)";
                    while (my $base = chop $seq) {
                    	if ($count++ % 50 == 0) {
                        	print "\n";
                    	}
                    	print $base;
                    }
                    print "\n";
            	}
            }
        }
        else {
            $seq .= $_ ;
        }
    }
}
