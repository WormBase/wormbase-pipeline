
use strict;

my (%probes, %wormbase_names, $probe_name, %sequences, %to_keep);

my $wormbase_probenames = shift(@ARGV);
open(my $nfh, $wormbase_probenames) or die "Could not open $wormbase_probenames for reading\n";

while(<$nfh>) {
  /^(\S+)/ and $wormbase_names{$1} = 1;
}


while(<>) {
  /^\>(\S+)/ and do {
    $probe_name = $1;
    # sometimes you get duplicates in a file. Not sure why this is, but only take the first one
    if (exists $probes{$probe_name}) {
      $probe_name = "";
    }
    next;
  };

  /^(\S+)/ and do {
    if ($probe_name) {
      $probes{$probe_name} .= uc($1);
    }
  };
}

foreach my $probe_id (keys %probes) {
  $sequences{$probes{$probe_id}}->{$probe_id} = 1;
}



foreach my $seq (keys %sequences) {
  my @pids = sort keys %{$sequences{$seq}};

  if (scalar(@pids) > 1) {
    # multiple probes with the same sequence; look for the one in WormBase
    my %found_in_wb;

    foreach my $pid (@pids) {
      my ($array, $pname) = split(/:/, $pid); 

      if (exists $wormbase_names{$pname} or
          exists $wormbase_names{"${array}_${pname}"}) {
        $found_in_wb{$pid} = 1;
      }
    }
    my ($first, @others) = keys %found_in_wb;

    if (@others) {
      print STDERR "Skipping @pids : multiple were found in WormBase\n";
      next;
    } elsif (not defined $first) {
      # no wormbase ids found, probably because this seq is redundant w.r.t. another seq 
      # on another array, and the probe name on that array was chosen by Wen as the representative 
      # for the seq. But we will keep one representative for this array
      print STDERR "No wormbase names for @pids, so choosing first as representative\n";
      $to_keep{$pids[0]} = 1; 
    } else {
      print STDERR "Choosing wormbase name $first as representative\n";
      $to_keep{$first} = 1;
    }
  } else {
    $to_keep{$pids[0]} = 1; 
  }  
}

foreach my $pid (sort keys %probes) {
  if (exists $to_keep{$pid}) {
    print ">$pid\n";
    print "$probes{$pid}\n";
  }
}
