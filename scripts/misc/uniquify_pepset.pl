#!/usr/local/bin/perl5.8.0 -w  

my $pep;
while (<>) {
  if (/\>(\w+)/ ) {
    $pep = $1;
    $peps{$pep} = "";
  }
  else {
    $peps{$pep} .= $_;
  }
}

foreach (sort keys %peps) {
  print "\>$_\n$peps{$_}";
}
