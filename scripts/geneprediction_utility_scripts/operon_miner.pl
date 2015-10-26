#!/bin/env perl

use Getopt::Long;

my ($infile,$gff,$ace);
GetOptions(
    'infile=s' => \$infile,
    'gff'   => \$gff,
    'ace' => \$ace,
)||die($!);

open INF,$infile;

# $genes{seq}=@{%{start => 1, stop => 1, strand => '+', name => 'abc'},}
my %genes;
my %genes2operon;
my %operons;

# read genes
while(<INF>){
 if (/gene\tgene/){
   chomp;
   my @F=split;
   push(@{$genes{$F[0]}},{start => $F[3],stop => $F[4],strand => $F[6],name => $F[9]});
 }
}

seek(INF,0,0);

my $count = 1;
while(<INF>){
  if (/cufflinks\s+Transcript/){
    chomp;
    my @F=split;
    my $operon = "BMOP$count";
    my @genesHit = grep {$_->{start} <= $F[4] && $_->{stop} >= $F[3] && $_->{strand} eq $F[6] } @{$genes{$F[0]}};
    next unless scalar @genesHit > 1;
    my @genesHit = sort {$a->{start} <=> $b->{start}} @genesHit;
    map {$operon = $genes2operon{$_->{name}} if $genes2operon{$_->{name}}} @genesHit;
    map {$genes2operon{$_->{name}}=$operon} @genesHit;
    if ($operons{$operon}){
        ${$operons{$operon}}->{start} = $genesHit[0]->{start} if $genesHit[0]->{start} < ${$operons{$operon}}->{start};
        ${$operons{$operon}}->{stop} = $genesHit[-1]->{stop} if $genesHit[-1]->{stop} > ${$operons{$operon}}->{stop};
    }
    else {
        ${$operons{$operon}}->{start} = $genesHit[0]->{start};
        ${$operons{$operon}}->{stop} = $genesHit[-1]->{stop};
        ${$operons{$operon}}->{strand} = $genesHit[0]->{strand};
        ${$operons{$operon}}->{seq} = $F[0];
        $count++;
    }
    ${$operons{$operon}}->{members} = {} unless ${$operons{$operon}}->{members};
    map { ${$operons{$operon}}->{members}->{$_->{name}} = 1} @genesHit;
  }
}

while (my($k,$v)=each %operons){
   if ($gff){
    print "$$v->{seq}\toperon\toperon\t$$v->{start}\t$$v->{stop}\t.$vv->{strand}\t\.\tOperon \"$k\"\n";
   }
   elsif($ace){
    my ($start,$stop) = $$v->{strand} eq '+'?($$v->{start},$$v->{stop}):($$v->{stop},$$v->{start});
    print "Sequence : $$v->{seq}\nOperon $k $start $stop\n\n";
    print "Operon : $k\nCanonical_parent $$v->{seq}\nSpecies \"Brugia malayi\"\n";
    map {print "Contains_Gene $_\n"} keys %{$$v->{members}};
    print "Description \"Predicted operon from cufflinks transcripts\"\nMethod operon\n\n";
   }
   else{
    print "$k\t$$v->{seq}\t$$v->{start}\t$$v->{stop}\t$$v->{strand}\t";
    print join "\t",keys %{$$v->{members}};
    print "\n"
   }
 
}
