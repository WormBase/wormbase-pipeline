#!/usr/bin/env perl
# rename CDSes, so the longest one is a

use Ace;
use Ace::Sequence;

my $db=Ace->connect(-path => shift);

my @suffix=qw(a b c d e f g h i j k l m n o p q r s t u v w x y z);

my $it = $db->fetch_many(Gene => '*');
while (my $gene = $it->next){
   next unless $gene->Corresponding_CDS;
   my @CDSes = $gene->Corresponding_CDS;
   my $n=0;
   foreach my $cds (sort {Ace::Sequence->new($b)->length() <=> Ace::Sequence->new($a)->length()} @CDSes){
    my $cname;
    if (scalar(@CDSes) == 1){
       $cname = $gene->Sequence_name;
    }else{
       $cname = $gene->Sequence_name . $suffix[$n];
       $n++;
    }
    print "-R CDS : $cds _$cname\n" unless $cname eq $cds;
   }
}
