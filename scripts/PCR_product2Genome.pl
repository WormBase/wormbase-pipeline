#!/software/bin/perl -w

# Version: $Version: $
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-07-15 10:47:01 $

use strict;
use warnings;
use Getopt::Long;

use lib $ENV{'CVS_DIR'};

use Wormbase;
use Log_files;

use Ace;
use Coords_converter;

################################
# command-line options         #
################################

my ($debug,
    $test,
    $store, 
    $species,
    $load,
    $database,
    $no_load,
    $target,
    $acefile,
    $prior_mappings,
    $EPCR, 
    $IPCRESS,
    $MAX_PROD_SIZE,
    $MIN_PROD_SIZE,
    $only_unmapped,
    @prod_test,
    %prod_test,
    );

GetOptions(
  "debug=s"         => \$debug,
  "test"            => \$test,
  "species:s"       => \$species,
  "store=s"         => \$store,
  "database=s"      => \$database,
  "noload"          => \$no_load,
  "target=s"        => \$target,
  "acefile=s"       => \$acefile,
  "ipcress=s"       => \$IPCRESS,
  "epcr=s"          => \$EPCR,
  "minsize=s"       => \$MIN_PROD_SIZE,
  "maxsize=s"       => \$MAX_PROD_SIZE,
  "onlyunmapped"    => \$only_unmapped,
  "prodtest=s@"     => \@prod_test,
    );



############################
# recreate configuration   #
############################
my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("cant restore wormbase from $store\n"); 
}
else { 
  $wormbase = Wormbase->new( -debug => $debug, 
                             -test => $test, 
                             -organism => $species); 
}

my $log = Log_files->make_build_log($wormbase);
###########################################

$IPCRESS = "ipcress" if not defined $IPCRESS;
$EPCR = "e-PCR" if not defined $EPCR;
$MIN_PROD_SIZE = 50 if not defined $MIN_PROD_SIZE;
$MAX_PROD_SIZE = 50000 if not defined $MAX_PROD_SIZE;
$database = $wormbase->autoace if not $database;
$acefile  = $wormbase->acefiles . "/PCR_products.mappings.ace" if not $acefile;

$target   = $wormbase->genome_seq if not $target;
$log->log_and_die("Could not find $database\n") if not -d $database;
$log->log_and_die("Could not find $target\n") if not -e $target;

my $coords = Coords_converter->invoke($database, 0, $wormbase);

my (%mapped, $missed_by_epcr, $total_count, $fail_count);

my $prods_to_map = &get_pcr_products();
$total_count = scalar(keys %$prods_to_map);

$log->write_to(sprintf("%d unmapped products to map\n\n", scalar(keys %$prods_to_map)));

#
# ipcress round 1
#
&map_with_ipcress($prods_to_map, \%mapped, 1);

foreach my $id (keys %$prods_to_map) {
  if (exists $mapped{$id}) {
    delete $prods_to_map->{$id};
  }
}

#
# e-PCR
#
my $missed = scalar(keys %$prods_to_map);
$log->write_to("$missed products were missed by ipcress round 1.\n");
if ($missed) {
  $log->write_to("Trying e-PCR\n");

  &map_with_epcr($prods_to_map, \%mapped);

  foreach my $id (keys %$prods_to_map) {
    if (exists $mapped{$id}) {
      delete $prods_to_map->{$id};
    }
  }
}


#
# icress round 2 and 3
#
foreach my $mm (2, 3) {
  $missed =  scalar(keys %$prods_to_map);
  $log->write_to("$missed  products remain missed.\n");
  if ($missed) {
    $log->write_to(sprintf("Trying ipcress with %d mismatches\n", $mm)); 

    &map_with_ipcress($prods_to_map, \%mapped, $mm);

    foreach my $id (keys %$prods_to_map) {
      if (exists $mapped{$id}) {
        delete $prods_to_map->{$id};
      }
    }
  }
}


#
# log the outstanding unmapped products
#
$fail_count = scalar(keys %$prods_to_map);
$log->write_to(sprintf("\nFailed to map %d products (%.1f percent)\n\n", 
                       $fail_count, 100 * ($fail_count / $total_count)));
foreach my $id (keys %$prods_to_map) {
  $log->write_to(sprintf("Could not map %s %s %s\n", 
                         $id, 
                         $prods_to_map->{$id}->{left_primer}, 
                         $prods_to_map->{$id}->{right_primer}));
}

# Write and load an ace file for the mapped products
#
my (%prods_by_loc, %non_primary_hits);

foreach my $mp (sort keys %mapped) {

  my (@best);
  foreach my $hit (sort { $a->{edits} <=> $b->{edits} } @{$mapped{$mp}}) {
    if (not @best or ($best[-1]->{edits} == $hit->{edits})) {
      push @best, $hit;
    } else {
      last;
    }
  }

  if (scalar(@best) > 1) {

    # all @best now have same number of total mismatches. If one of the hits can 
    # is on the "right" clone (as determined by the clone name being part of the
    # product name), choose it as the primary hit; otherwise, choose the one with 
    # the shortest implied product length
    
    my (@good_looking_seq, @others);
    foreach my $hit (@best) {
      my ($seq, $st, $en) = $coords->LocateSpan($hit->{chr}, $hit->{start}, $hit->{end});

      if ($mp =~ /$seq/) {
        push @good_looking_seq, $hit;
      } else {
        push @others, $hit;
      }
    }
    if (scalar(@good_looking_seq) == 1) {
      @best = (@good_looking_seq, @others)
    } else {
      @best = sort { $a->{len} <=> $b->{len} } @best;    
    }


    # keep all matches iff they are perfect; otherwise just take the primary one
    if ($best[0]->{edits} > 0) {
      @best = ($best[0]);
    }

    # keep all matches that have length within a sensible range of the primary
    # (this weeds out absurd non-primary hits)
    @best = grep { $_->{len} < ($best[0]->{len} * 2) } @best;
  
    # Finally, mark all but the first hit as non-primary
    for(my $idx = 1; $idx < @best; $idx++) {
      $best[$idx]->{product} .= "_${idx}";
      
      push @{$non_primary_hits{$mp}}, {
        product => $best[$idx]->{product},
        method  => $best[$idx]->{method},
      };
    }

  }


  foreach my $hit (@best) {
    my ($seq, $st, $en) = $coords->LocateSpan($hit->{chr}, $hit->{start}, $hit->{end});
    
    $hit->{chr} = $seq;
    $hit->{start} = $st;
    $hit->{end} = $en;
    
    push @{$prods_by_loc{$seq}}, $hit;
  }
}


open (my $acefh, ">$acefile") or $log->log_and_die("Could not open $acefile for writing\n");

foreach my $prod (sort keys %non_primary_hits) {
  foreach my $obj (@{$non_primary_hits{$prod}}) {
    my $name = $obj->{product};
    my $method = $obj->{method};
    print $acefh "\nPCR_product : \"$name\"\n";
    print $acefh "Remark \"Created to represent non-primary genome mapping of $prod\"\n";
    print $acefh "Method \"$method\"\n";
  }
}

foreach my $chr (sort keys %prods_by_loc) {
  print $acefh "\nSequence : \"$chr\"\n";
  foreach my $mp (sort { $a->{product} cmp $b->{product} } @{$prods_by_loc{$chr}}) {

    printf $acefh "PCR_product\t%s\t%d\t%d\n", $mp->{product}, $mp->{start}, $mp->{end};
  }
}
close($acefh);

if (not $no_load) {
  $wormbase->load_to_database($database, $acefile, "PCR_product_mapping", $log);
}


$log->mail;
exit;


#############################
sub map_with_ipcress {
  my ($products, $mapped, $mismatches) = @_;

  my $ipfile = "/tmp/PCR_products.WB.$$.ipcress.txt";
  open(my $ipfh, ">$ipfile") or $log->log_and_die("Could not open ipcress file for reading\n");
  foreach my $prod_name (sort keys %$products) {
    printf($ipfh "%s\t%s\t%s\t%d\t%d\n", 
           $prod_name, 
           $products->{$prod_name}->{left_primer},
           $products->{$prod_name}->{right_primer},
           $MIN_PROD_SIZE, 
           $MAX_PROD_SIZE);
  }
  close($ipfh);
  
  open(my $iprunfh, "$IPCRESS --memory 512 --mismatch $mismatches --pretty FALSE $ipfile $target |") 
      or $log->log_and_die("Could not open ipcress command\n");
  while(<$iprunfh>) {
    /ipcress:\s+(\S+)\s+(\S+)\s+(\d+)\s+\S\s+(\d+)\s+(\d+)\s+\S\s+(\d+)\s+(\d+)\s+(\S+)/ and do {
      my ($tseq, $prod_id, $prod_len, $start, $mm1, $start2, $mm2, $ori) = ($1, $2, $3, $4, $5, $6, $7, $8);

      next if $ori =~ /single/;

      $start++;
      my $end = $start + $prod_len - 1;

      #if ($ori =~ /revcomp/) {
      #  ($start, $end) = ($end, $start);
      #}

      $tseq =~ s/:filter\S+//; 
      my $mm = $mm1 + $mm2;

      my $obj = {
        product => $prod_id,
        method => $products->{$prod_id}->{method},
        chr => $tseq,
        start => $start,
        end => $end, 
        len => $prod_len,
        edits => $mm,
      };

      push @{$mapped->{$prod_id}}, $obj;
    }
  }

  unlink $ipfile;
}




#############################
sub map_with_epcr {
  my ($products, $mapped) = @_;

  my $epcrfile = "/tmp/PCR_products.WB.$$.epcr.txt";
  open(my $epcrfh, ">$epcrfile") or $log->log_and_die("Could not open epcr file for reading\n");
  foreach my $prod_name (sort keys %$products) {
    printf($epcrfh "%s\t%s\t%s\t%d-%d\n", 
           $prod_name, 
           $products->{$prod_name}->{left_primer},
           $products->{$prod_name}->{right_primer},
           $MIN_PROD_SIZE, 
           $MAX_PROD_SIZE);
  }
  close($epcrfh);

  open(my $epcrrunfh, "$EPCR -n 1 -g 1 -t 3 $epcrfile $target |")
      or $log->log_and_die("Could not open e-PCR commmand\n");
  while(<$epcrrunfh>) {
    /^(\S+)\s+(\S+)\s+(\S)\s+(\d+)\s+(\d+)\s+(\d+)\/\d+\-\d+\s+(\d+)\s+(\d+)/ and do {

      my $obj = {
        product => $2,
        method => $products->{$2}->{method},
        chr => $1,
        start => $4, 
        end   => $5,
        len => $6, 
        edits => $7 + $8,
      };

      push @{$mapped->{$2}}, $obj;
    }
  }

  unlink $epcrfile;
}

#############################
sub get_pcr_products {

  # General strategy
  # 1. Fetch products that have Left_mapping_primer and Right_mapping_primer set
  #    (these are cases where we have added our own primer sequences to the object)
  # 2. Fetch the rest

  if (@prod_test) {
    map { $prod_test{$_} = 1 } @prod_test;
  }

  my (%good_prods);
  my $tace = $wormbase->tace;

  my $tb_def1 = &get_initial_table_def();
  my $tb_cmd1 = "Table-maker -p \"$tb_def1\"\nquit\n";
  open(my $tace_fh1, "echo '$tb_cmd1' | $tace $database |");
  while(<$tace_fh1>) {
    if (/^\"(\S+)\"/) {
      s/\"//g;

      my @l = split;

      next if @prod_test and not exists $prod_test{$l[0]};

      $good_prods{$l[0]} = {
        method       => $l[1],
        left_primer  => uc($l[2]),
        right_primer => uc($l[3]),
      };
    }
  }
  close($tace_fh1);

  my (%prods, %meths);
  
  my $tb_def2 = &get_secondary_table_def();
  my $tb_cmd2 = "Table-maker -p \"$tb_def2\"\nquit\n";
  open(my $tace_fh2, "echo '$tb_cmd2' | $tace $database |");
  while(<$tace_fh2>) {
    if (/^\"(\S+)\"/) {
      s/\"//g;

      my @l = split;

      next if @prod_test and not exists $prod_test{$l[0]};

      $meths{$l[0]} = $l[1];

      if ($l[2] and $l[3]) {
        push @{$prods{$l[0]}}, [$l[2], $l[3]];
      } elsif ($l[2]) {
        push @{$prods{$l[0]}}, [$l[2], ""];
      } else {
        $prods{$l[0]} = [];
      }
    }
  }
  close($tace_fh2);

  foreach my $prod (sort keys %prods) {
    if (scalar(@{$prods{$prod}}) != 2) {
      $log->write_to(sprintf("Product $prod did not have 2 oligos ( %d ), so skipping\n", scalar(@{$prods{$prod}})));
      next;
    } elsif (not $prods{$prod}->[0]->[1] or not $prods{$prod}->[1]->[1]) {
      $log->write_to(sprintf("Product $prod has empty oligo sequence for one or both of its oligos; skipping\n"));
      next;
    }

    $good_prods{$prod} = {
      left_primer  => uc($prods{$prod}->[0]->[1]),
      right_primer => uc($prods{$prod}->[1]->[1]),
      method       => $meths{$prod},
    };
  }

  unlink $tb_def1, $tb_def2;
  return \%good_prods;
}


#############################
sub get_secondary_table_def {

  my $tmdef = "/tmp/pcr_prod_2.$$.def";

  open(my $qfh, ">$tmdef") or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $condition = 
      "NOT Mapping_primers AND Method AND Method != \"Variation_PCR\"";
  if ($only_unmapped) {
    $condition .= " AND NOT Canonical_parent";
  }

  my $query = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class PCR_product 
From 1
Condition $condition
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class Method 
From 1 
Tag Method 

Colonne 3 
Width 12 
Mandatory 
Visible 
Class 
Class Oligo 
From 1 
Tag Oligo 
 
Colonne 4 
Width 12 
Mandatory 
Visible 
Text 
From 3 
Tag Sequence 



EOF

  print $qfh $query;
  close($qfh);

  return $tmdef;
}

################################
sub get_initial_table_def {

  my $tmdef = "/tmp/pcr_prod_1.$$.def";

  open(my $qfh, ">$tmdef") or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $condition = ($only_unmapped) ? "Condition NOT Canonical_parent" : "";

  my $query = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class PCR_product 
From 1 
$condition
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class
Class Method 
From 1 
Tag Method
 
Colonne 3 
Width 12 
Mandatory 
Visible 
Text 
From 1 
Tag Left_mapping_primer 
 
Colonne 4 
Width 12 
Mandatory 
Visible 
Text 
From 1 
Tag Right_mapping_primer 


EOF

  print $qfh $query;
  close($qfh);

  return $tmdef;
}
