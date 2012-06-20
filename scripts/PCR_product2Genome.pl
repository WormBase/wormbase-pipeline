#!/software/bin/perl -w

# Version: $Version: $
# Last updated by: $Author: klh $
# Last updated on: $Date: 2012-06-20 13:00:14 $

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
    $previous_mapping_dir,
    $load,
    $database,
    $no_load,
    $target,
    $acefile,
    $EPCR, 
    $IPCRESS,
    $MAX_PROD_SIZE,
    $MIN_PROD_SIZE,
    );

GetOptions(
  "debug=s"        => \$debug,
  "test"           => \$test,
  "species:s"      => \$species,
  "store=s"        => \$store,
  "previousdir=s"  => \$previous_mapping_dir,
  "database=s"     => \$database,
  "noload"         => \$no_load,
  "target=s"       => \$target,
  "acefile=s"      => \$acefile,
  "ipcress=s"      => \$IPCRESS,
  "epcr=s"         => \$EPCR,
  "minsize=s"      => \$MIN_PROD_SIZE,
  "maxsize=s"      => \$MAX_PROD_SIZE,
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

my (%mapped, $missed_by_ipcress, $missed_by_epcr, $total_count, $fail_count);

my $prods_to_map = &get_pcr_products();
$total_count = scalar(keys %$prods_to_map);
&map_with_ipcress($prods_to_map, \%mapped);

foreach my $id (keys %$prods_to_map) {
  if (exists $mapped{$id}) {
    delete $prods_to_map->{$id};
  } else {
    $missed_by_ipcress++;
  }
}

$log->write_to("$missed_by_ipcress products were not found be ipcress - Trying e-PCR\n");

&map_with_epcr($prods_to_map, \%mapped);

foreach my $id (keys %$prods_to_map) {
  if (exists $mapped{$id}) {
    delete $prods_to_map->{$id};
  }
}

#
# log the outstanding unmapped products
#
my $fail_count = scalar(keys %$prods_to_map);
$log->write_to("Failed to map %d products (%d percent)\n", $fail_count, 100 * ($fail_count / $total_count));
foreach my $id (keys %$prods_to_map) {
  $log->write_to(sprintf("Could not map %s %s %s\n", $id, @{$prods_to_map->{$id}}));
}

# Write and load an ace file for the mapped products
#
my %prods_by_loc;

foreach my $mp (sort keys %mapped) {

  my (@best, @best2);
  foreach my $hit (sort { $a->{edits} <=> $b->{edits} } @{$mapped{$mp}}) {
    if (not @best or ($best[-1]->{edits} == $hit->{edits})) {
      push @best, $hit;
    } else {
      last;
    }
  }

  foreach my $hit (sort { $a->{len} <=> $b->{len} } @best) {
    my $overlap = 0;
    foreach my $ohit (@best2) {
      if ($ohit->{chr} eq $hit->{chr} and 
          $hit->{start} <= $ohit->{end} and 
          $hit->{end} >= $ohit->{start}) {
        $overlap = 1;
        last;
      } 
    }
    if (not $overlap) {
      push @best2, $hit;
    }
  }

  foreach my $hit (@best2) {
    my ($seq, $st, $en) = $coords->LocateSpan($hit->{chr}, $hit->{start}, $hit->{end});

    $hit->{chr} = $seq;
    $hit->{start} = $st;
    $hit->{end} = $en;

    push @{$prods_by_loc{$seq}}, $hit;
  }
}


open (my $acefh, ">$acefile") or $log->log_and_die("Could not open $acefile for writing\n");
foreach my $chr (sort keys %prods_by_loc) {
  print $acefh "\nSequence : \"$chr\"\n";
  foreach my $mp (sort { $a->{product} cmp $b->{product} } @{$prods_by_loc{$chr}}) {
    printf $acefh "PCR_product\t%s\t%d\t%d\n", $mp->{product}, $mp->{start}, $mp->{end};
  }
}
close($acefh);

$log->mail;
exit;


#############################
sub map_with_ipcress {
  my ($products, $mapped) = @_;

  my $ipfile = "/tmp/PCR_products.WB.$$.ipcress.txt";
  open(my $ipfh, ">$ipfile") or $log->log_and_die("Could not open ipcress file for reading\n");
  foreach my $prod_name (sort keys %$products) {
    printf($ipfh "%s\t%s\t%s\t%d\t%d\n", 
           $prod_name, 
           @{$products->{$prod_name}}, 
           $MIN_PROD_SIZE, 
           $MAX_PROD_SIZE);
  }
  close($ipfh);
  
  open(my $iprunfh, "$IPCRESS --memory 512 --mismatch 1 --pretty FALSE $ipfile $target |") 
      or $log->log_and_die("Could not open ipcress command\n");
  while(<$iprunfh>) {
    /ipcress:\s+(\S+)\s+(\S+)\s+(\d+)\s+\S\s+(\d+)\s+(\d+)\s+\S\s+(\d+)\s+(\d+)\s+(\S+)/ and do {
      my ($tseq, $prod_id, $prod_len, $start, $mm1, $end, $mm2, $ori) = ($1, $2, $3, $4, $5, $6, $7, $8);

      next if $ori =~ /single/;

      $start++;
      $end = $start + $prod_len + 1;

      if ($ori =~ /revcomp/) {
        ($start, $end) = ($end, $start);
      }

      $tseq =~ s/:filter\S+//; 
      my $mm = $mm1 + $mm2;


      my $obj = {
        product => $prod_id,
        chr => $tseq,
        start => $start,
        end => $end, 
        len => $prod_len,
        edits => $mm2,
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
           $prod_name, @{$products->{$prod_name}}, 
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
        chr => $1,
        start => ($3 eq '-') ? $5 : $4, 
        end   => ($3 eq '-') ? $4 : $5, 
        len => $6, 
        edits => $7 + $8,
      };

      push @{$mapped->{$2}}, $obj;
    }
  }
}

#############################
sub get_pcr_products {

  #return {
  #  # found by ipcress
  #  'sjj_AC8.3'    => ['CAATTGTGCGAAATAAGTGTCAA', 'CGCTTTCGAGAGTACGTCTATGT'],
  #  'sjj_AC8.4'    => ['TATTACTGGCACTTTGGCTCATT', 'TGGTTATGGAGAAGAGCACGTAT'],
  #  'sjj_B0041.5'  => ['CTGTCTCCAGTCTCAATCTCGTT', 'CGCTGTATTAATTTTCACTTCCG'],
  #  'sjj_B0280.1'  => ['CTGCATCCGAGAAACTGTCA',    'TCAACGCGATGGATTTATCA'],
    # missed by ipcress, found by epcr
  #  'cenix:300-c5' => ['GATTTTCCTTACCCTTCATGTCC', 'AGACACCTGGTGGTAATAGGGAT'], 
    # missed by both
  #  'mv_B0284.4'   => ['AAGTGGTATTTCCTCCATTTCCATTG', 'TGGGCGCTAATGGAAGCACAGA'],
  #};

  my (%prods, %good_prods);

  my $tb_def = &get_table_def();
  my $tace = $wormbase->tace;

  my $tb_cmd = "Table-maker -p \"$tb_def\"\nquit\n";
  open(my $tace_fh, "echo '$tb_cmd' | $tace $database |");
  while(<$tace_fh>) {
    if (/^\"(\S+)\"/) {
      s/\"//g;
      my @l = split;

      if ($l[1] and $l[2]) {
        push @{$prods{$l[0]}}, [$l[1], $l[2]];
      } elsif ($l[1]) {
        push @{$prods{$l[0]}}, [$l[1], ""];
      } else {
        $prods{$l[0]} = [];
      }
    }
  }

  foreach my $prod (sort keys %prods) {
    if (scalar(@{$prods{$prod}}) != 2) {
      $log->write_to(sprintf("Product $prod did not have 2 oligos ( %d ), so skipping\n", scalar(@{$prods{$prod}})));
      next;
    } elsif (not $prods{$prod}->[0]->[1] or not $prods{$prod}->[1]->[1]) {
      $log->write_to(sprintf("Product $prod has empty oligo sequence for one or both of its oligos; skipping\n"));
      next;
    }

    $good_prods{$prod} = [uc($prods{$prod}->[0]->[1]), uc($prods{$prod}->[1]->[1])];
  }

  unlink $tb_def;
  return \%good_prods;
}


#############################
sub get_table_def {
  my ($only_unmapped) = @_;

  my $tmdef = "/tmp/pcr_prod.$$.def";

  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $condition = ($only_unmapped) 
      ? "Condition Method AND Method != \"Variation_PCR\" AND NOT Sequence" 
      : "Condition Method AND Method != \"Variation_PCR\"";

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
Class Oligo 
From 1 
Tag Oligo 
 
Colonne 3 
Width 12 
Mandatory 
Visible 
Text 
From 2 
Tag Sequence 

EOF

  print $qfh $query;
  close($qfh);

  return $tmdef;
}

