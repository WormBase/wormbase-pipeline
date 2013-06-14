#!/usr/bin/env perl

# Version: $Version: $
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-06-14 10:38:57 $

use strict;
use warnings;
use Getopt::Long;

use lib $ENV{'CVS_DIR'};

use Wormbase;
use Log_files;

use strict;

my ($debug,
    $test,
    $store, 
    @gff_current, 
    @gff_previous,
    );

GetOptions(
  "debug=s"        => \$debug,
  "test"           => \$test,
  "store=s"        => \$store,
  "currentgff=s@"   => \@gff_current,
  "previousgff=s@"  => \@gff_previous,
    );


my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("cant restore wormbase from $store\n"); 
}
else { 
  $wormbase = Wormbase->new( -debug => $debug, 
                             -test => $test );
}

my $log = Log_files->make_build_log($wormbase);

if (not @gff_current or grep { not -e $_ } @gff_current) {
  $log->log_and_die("You must supply at least one valid current GFF file with -currentgff\n");
}
if (not @gff_previous or grep { not -e $_ } @gff_previous ) {
  $log->log_and_die("You must supply at least one valid previous GFF file with -previousgff\n");
}

my (%previous_data, %current_data);

foreach my $pair ([\@gff_current, \%current_data], 
                  [\@gff_previous, \%previous_data]) {
  my ($files_arr, $data) = @$pair;

  foreach my $file (@$files_arr) {
    my $fh;
    if ($file =~ /\.gz$/) {
      open( $fh, "gunzip -c $file |") or $log->log_and_die("Could not open gunzip stream to $file\n");
    } else {
      open( $fh, $file) or $log->log_and_die("Could not open $file for reading\n");
    }
    while(<$fh>) {
      /^\#/ and next;
      
      /^\S+\s+(\S+)\s+(\S+)/ and do {
        my ($src, $type) = ($1, $2);
        
        $data->{$src . " " . $type}++;
      }
    }
  }
}

my (%all_keys, %only_in_current, %only_in_previous, %data_loss, %data_gain, %identical);
foreach my $key (keys %previous_data, keys %current_data) {
  $all_keys{$key} = 1;
}

foreach my $key (sort keys %all_keys) {
  if (exists $previous_data{$key} and not exists $current_data{$key}) {
    $only_in_previous{$key} = 1;
  } elsif (exists $current_data{$key} and not exists $previous_data{$key}) {
    $only_in_current{$key} = 1;
  } elsif ($previous_data{$key} > $current_data{$key}) {
    $data_loss{$key} = 1;
  } elsif ($previous_data{$key} < $current_data{$key}) {
    $data_gain{$key} = 1;
  } else {
    $identical{$key} = 1;
  }
}

$log->write_to("\nGFF report: comparing gff_current (@gff_current) with gff_previous (@gff_previous)\n\n");
$log->write_to("Source/Types present previously but now missing (*INVESTIGATE*):\n\n");
foreach my $k (sort keys %only_in_previous) {
  $log->write_to(sprintf(" %-20s (%d)\n", $k, $previous_data{$k}));
}
$log->write_to("\nSource/Types present in higher numbers previously than now (*INVESTIGATE*):\n\n");
foreach my $k (sort keys %data_loss) {
  $log->write_to(sprintf(" %-20s %d (%d)\n", $k, $current_data{$k}, $previous_data{$k}));
}
$log->write_to("\nSource/Types present now but not present previously (probably okay):\n\n");
foreach my $k (sort keys %only_in_current) {
  $log->write_to(sprintf(" %-20s %d\n", $k, $current_data{$k}));
}
$log->write_to("\nSource/Types present in lower numbers previously than now (probably okay):\n\n");
foreach my $k (sort keys %data_gain) {
  $log->write_to(sprintf(" %-20s %d (%d)\n", $k, $current_data{$k}, $previous_data{$k}));
}
$log->write_to(sprintf("\nAll other Source/Types (%d) have identical counts between previous and current\n", scalar(keys %identical)));

$log->mail();
exit(0);
