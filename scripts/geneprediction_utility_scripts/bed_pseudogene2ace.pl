#!/usr/bin/env perl

my %block;

while(<>){
  chomp;
  my @F=split;
  my $gstart=$F[5] eq '+'?$F[1]+1:$F[2];
  my $gstop =$F[5] eq '+'?$F[2]:$F[1]+1;

  print "Sequence : $F[0]\n";
  print "Pseudogene $F[3] ",$gstart," ",$gstop,"\n";
  print "\n";

  print "Pseudogene : $F[3]\n";
  print "Sequence $F[0]\n";

  my @starts = split(',',$F[11]);
  my @sizes  = split(',',$F[10]);

#  printf "Source_exons 1 %d\n", abs($gstart-$gstop)+1;

  # based on CStarts in BED, based on SMap first value in ACe :-(

  for(my $n=0;$n < scalar(@starts);$n++){
   if ($F[5] eq '+'){
      printf "Source_exons %d %d\n",$starts[$n]+1,$starts[$n]+$sizes[$n];     
   }else{
      my $end = $F[2] -$F[1];
      printf "Source_exons %d %d\n", $end - ($starts[$n]+$sizes[$n])+1, $end - $starts[$n] +1; 
   }
 }
  print "Method Pseudogene\n";# remove that later
  print "\n";
}
