#!/usr/local/bin/perl

# uses the file provided by Peter Sterk (EBI) to write an ace file
# containing containing: the parent object, the protein_id without the
# version number, the version number


use strict;
$| = 1;

while (<>) {
    chomp;
    if (/\-/) {
        next;
    }
    elsif (/^(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)$/) {
        undef my ($acc , $pid , $ver , $gene , $id);
	$acc = $1;
	$pid = $2;
	$ver = $3;
	$gene = $4;
        open (ID , "getz5 \'\[embl\-acc:$acc\]\' \-f id |");
        while (<ID>) {
           chomp;
           /ID\s+CE(\w+)\s+/;
           $id = $1;
           }
	close ID;
        if (($id eq "") || (!defined $id)) {
            open (ID , "getz5 \'\[emblnew\-acc:$acc\]\' \-f id |");
            while (<ID>) {
               chomp;
               /ID\s+CE(\w+)\s+/;
               $id = $1;
	    }
	    close ID;
            if (($id eq "") || (!defined $id)) {
                warn "no parent for gene $gene , accession $acc\n";
	        next;
            }

	}
        print "Sequence $gene\nProtein_id $id $pid $ver\n\n";
    }
}

