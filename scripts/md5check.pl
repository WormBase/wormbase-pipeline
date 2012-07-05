#!/usr/bin/env perl

use Ace;
use Digest::MD5 qw(md5_hex);
use Getopt::Long;
use strict;

my $species = 'Caenorhabditis japonica';
my ($database,$write);

GetOptions('species=s'  => \$species,
	   'database=s' => \$database,
           'overwrite'  => \$write,	
);
my $db = Ace->connect(-path => $database);
my $seqIt = $db->fetch_many(-query => "Find Sequence WHERE Species=\"$species\" AND Method=\"Genomic_canonical\"");
while(my $seq=$seqIt->next){
        my $ace_md5=$seq->MD5;
        my $seq_md5=md5_hex(uc($seq->asDNA));
        if ($seq_md5 eq $ace_md5){
		print "$seq MD5 ok\n" unless $write;
	}else{
		print "ERROR: $seq has a wrong MD5 checksum\n" unless $write;
		# $db->parse("Sequence : $seq\nMD5 $seq_md5\n\n")||print Ace->error if $write;
		print "Sequence : $seq\nMD5 $seq_md5\n\n" if $write;
	}
}

