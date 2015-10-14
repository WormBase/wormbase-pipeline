#!/usr/bin/env perl
# fix the genespan based on it's connected CDSes
#
use Ace;
use Ace::Sequence;
use strict;

my $db = Ace->connect(-path => shift()) ||die(Ace->Error);
my $genes = $db->fetch_many(-query => 'Find Gene Species="Brugia malayi"');
while (my $gene = $genes->next){
    my ($start,$stop,$seq);
    next unless $gene->Species eq 'Brugia malayi';

    my @children;
    if ($gene->Corresponding_Transcript){
      @children = $gene->Corresponding_Transcript;
    }elsif($gene->Corresponding_CDS){
      @children = $gene->Corresponding_CDS;
    }elsif($gene->Corresponding_Pseudogene){
      @children = $gene->Corresponding_Pseudogene;
    }else{
      next;
    }

    foreach my $cds (@children) {
        my $c = Ace::Sequence->new($cds);
        die("cannot find parent of $cds\n") unless $c;
        $c->refseq($c->parent);
        $start||=$c->start;
        $stop||=$c->end;
        if ($c->start > $c->end){ #reverse
           $start = $c->start if $c->start > $start;
           $stop = $c->end if $c->end < $stop;
        }else{
           $start = $c->start if $c->start < $start;
           $stop = $c->end if $c->end > $stop;
        }
        $seq = $c->parent();
    }
    print<<HERE;
Sequence : $seq
Gene_child $gene $start $stop

Gene : $gene
Sequence $seq

HERE
}
