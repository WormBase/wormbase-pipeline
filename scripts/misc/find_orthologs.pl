#!/usr/local/bin/perl5.6.1 -w

$b_e = shift;
$e_b = shift;

open (BE,"<$b_e");
while (<BE>) {
  chomp;
  @data = split;
  $brig_el{$data[0]} = $data[1];
}


open (EB,"<$e_b");
while (<EB>) {
  chomp;
  @data = split;
  $el_brig{$data[0]} = $data[1];
  
  $ort = $brig_el{$data[1]};
  if( $data[0] eq  $ort) {
    print "$data[0] ortholog =  $data[1] \n";
  }
}
close BE;
close EB;
