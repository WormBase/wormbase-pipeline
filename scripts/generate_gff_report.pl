#!/usr/bin/env perl

# Version: $Version: $
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-12-06 16:47:19 $

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
    $current_label,
    @gff_previous,
    $previous_label,
    $report_file,
    $checksum,
    $final,
    );

GetOptions(
  "debug=s"        => \$debug,
  "test"           => \$test,
  "store=s"        => \$store,
  "currentgff=s@"   => \@gff_current,
  "currentlabel=s"  => \$current_label,
  "previousgff=s@"  => \@gff_previous,
  "previouslabel=s" => \$previous_label,
  "report=s"        => \$report_file,
  "checksum"        => \$checksum,
  "final"           => \$final,
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

$current_label = "@gff_current" if not defined $current_label;
$previous_label = "@gff_previous" if not defined $previous_label;

my $report  = "\nCOMPARISON BETWEEN GFF CURRENT ($current_label) AND GFF PREVIOUS ($previous_label)\n\n";
my $identical = 0;

if ($checksum) {
  # can only do this id called with a single source and single dest
  if (scalar(@gff_previous) == 1 and
      scalar(@gff_current) == 1) {
    my ($prev) = @gff_previous;
    my ($cur)  = @gff_current;

    my ($md5_prev, $md5_cur, $md5fh);

    open($md5fh, "md5sum $prev |") or $log->log_and_die("Could not open md5sum command for $prev\n");
    while(<$md5fh>) {
      /^(\S+)/ and $md5_prev = $1;
    }
    open($md5fh, "md5sum $cur |") or $log->log_and_die("Could not open md5sum command for $cur\n");
    while(<$md5fh>) {
      /^(\S+)/ and $md5_cur = $1;
    }

    if ($md5_prev eq $md5_cur) {
      $report .= "Current and previous GFF files have identical checksums\n\n";
      $identical = 1;
    } else {
      $report .= "Current and previous GFF files have differing checksums\n\n";
    }
  }
}

if (not $identical) {
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
  
  if (keys %only_in_previous) {
    $report .= "\nWARNING: There are source/type present previously but now missing (*INVESTIGATE*):\n\n";
  
    foreach my $k (sort keys %only_in_previous) {
      $report .= sprintf(" %-20s (%d)\n", $k, $previous_data{$k});
    }
  } else {
    $report .= "OK: No source/type present previously but now missing\n";
  }

  if (not $final) {
    # Not interested in the count differences in final model
    if (keys %data_loss) {  
      $report .= "\nWARNING: There are source/types present in higher numbers previously than now (*INVESTIGATE*):\n\n";
      foreach my $k (sort keys %data_loss) {
        $report .= sprintf(" %-20s %d (%d)\n", $k, $current_data{$k}, $previous_data{$k});
      }
    } else {
      $report .= "OK: No source/type present in lower numbers now than previous\n";
    }
  }

  if (keys %only_in_current) {
    $report .= "\nWARNING: There are source/type present now but not present previously (probably okay):\n\n";
    foreach my $k (sort keys %only_in_current) {
      $report .= sprintf(" %-20s %d\n", $k, $current_data{$k});
    }
  } else {
    $report .= "OK: No source/type present now but not present previously\n";
  }

  if (not $final) {
    if (keys %data_gain) {
      $report .= "\nWARNING: There are source/type present in lower numbers previously than now (probably okay):\n\n";
      foreach my $k (sort keys %data_gain) {
        $report .= sprintf(" %-20s %d (%d)\n", $k, $current_data{$k}, $previous_data{$k});
      }
    } else {
      $report .= "OK: No source/type present in higher numbers now than previously\n";
    }
  }
    
  if (not $final) {
    $report .= sprintf("\nOK: Other source/types (%d) have identical counts between previous and current\n", scalar(keys %identical));
  } 

}

if ($report_file) {
  $log->write_to("Writing results to $report_file\n");
  open(my $fh, ">$report_file") or $log->log_and_die("Could not open report file for writing\n");
  print $fh $report;
  close($fh);
} else {
  print $report;
  $log->write_to($report);
}

$log->mail();
exit(0);
