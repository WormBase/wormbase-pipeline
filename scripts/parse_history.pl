#!/usr/local/bin/perl



$history = shift;
%genes = ();


open (HISTORY, "<$history");
while (<HISTORY>) {
    
    undef ($id);
    undef ($acc);
    undef ($start);
    undef ($end);
    undef ($changes);
    
    chomp;
    
    ($id,$acc,$start,$end) = split (/\t/,$_);
    
    $current = $genes{$id};
    $current++;
    $genes{$id} = $current;
	


}
close HISTORY;


foreach $key (keys %genes) {
    
    print "$key  \t: $genes{$key}\n";


}

exit(0);
