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


=pod

=head1 NAME

 uniquify_fasta.pl

 Very simple script to make a version of a fasta file unique ie no duplicate entries.

 Use like this

 WS123.pep | perl uniquify_fasta.pl > WS123.pep.unique

=cut
