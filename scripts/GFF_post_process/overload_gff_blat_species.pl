#!/usr/bin/env perl
#
# overload_gff_blat_species.pl
#
# Adds species to each non-native BLAT match
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-07-21 11:07:59 $

use strict;
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

######################################
# variables and command-line options #
######################################

my (  $debug, $test, $store, $wormbase, $species );
my ( $gff3, $infile, $outfile, $changed_lines );

GetOptions(
    'debug=s'   => \$debug,
    'test'      => \$test,
    'store:s'   => \$store,
    'species:s' => \$species,
    'infile:s'  => \$infile,
    'outfile:s' => \$outfile,
    'gff3'      => \$gff3,
);

if ($store) {
    $wormbase = retrieve($store) or croak("Can't restore wormbase from $store\n");
} 
else {
    $wormbase = Wormbase->new(
        -debug => $debug,
        -test  => $test,
	-organism => $species,
    );
}

my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n");
}

open(my $gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
open(my $gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  

my $speciesh = &get_species_data();

while (<$gff_in_fh> ) {
  if ( /^#/ ) {
    print $gff_out_fh $_;
    next;
  }
  chomp;
  my @f = split(/\t/, $_);

  if (grep { $f[1] =~ /^$_/ } ('BLAT_WASHU', 'BLAT_NEMBASE', 'BLAT_NEMATODE', 'BLAT_Caen_EST_') or 
      grep { $f[1] eq $_ } ('EMBL_nematode_cDNAs-BLAT', 'NEMATODE.NET_cDNAs-BLAT', 'NEMBASE_cDNAs-BLAT')) {
    
    # It is possible that this script is being run on the same
    # input file multiple times (e.g. when sorting out problems)
    # in which case we do not want to add 'Species' multiple
    # times to the same line.      
    
    if (defined $f[1] && $f[8] !~ /;\sSpecies/ && $f[8] !~ /Species\=/) {
      my $id;
      
      if ($gff3) {
        ($id) = $f[8] =~ /Target\=(\S+)/;
      } else {
        ($id) = ( $f[8] =~ /Target \"Sequence:(\S+)\"/ );
      }
      
      my $hkey = $f[1];
      # the {'BLAT_NEMATODE'} hash holds the EMBL data which BLAT_Caen_EST_* uses as well
      if ($hkey =~ /BLAT_Caen_/ or $hkey =~ /^EMBL_nematode_cDNAs/) {
        $hkey = "BLAT_NEMATODE";
      } elsif ($f[1] =~ /NEMATODE\.NET_cDNAs/) {
        $hkey = "BLAT_WASHU";
      } elsif ($f[1] =~ /^NEMBASE_cDNAs/) {
        $hkey = "BLAT_NEMBASE";
      }
      
      if ( exists $speciesh->{$hkey}->{$id} ) {
        my $suffix;
        if ($gff3) { 
          $f[8] .= ";Species=" . $speciesh->{$hkey}->{$id};
        } else {
          $f[8] .= " ; Species \"" . $speciesh->{$hkey}->{$id} . "\"";
        }
        $changed_lines++;
      } else {
        $log->write_to("Cannot find species info for $id\n");
      }
    }
  }
  
  print $gff_out_fh join("\t", @f), "\n";
}
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n");
$log->write_to("Finished processing : $changed_lines lines modified\n");
$log->mail();
exit(0);

#####################################
sub get_species_data {

  my $ace_dir = $wormbase->autoace;
  my $tace = $wormbase->tace;

  ###################################
  # get the species of the sequences
  ###################################

  my $cmd1 = "Query Find Sequence Where Database = \"NEMATODE_NET\"\nshow -a Species\nquit";
  my $cmd2 = "Query Find Sequence Where Database = \"NEMBASE\"\nshow -a Species\nquit";
  my $cmd3 = "Query Find Sequence Where Database = \"EMBL\"\nshow -a Species\nquit";

  my (%species, $id, $db );
  
  print "Finding BLAT_WASHU data\n";
  open( TACE, "echo '$cmd1' | $tace $ace_dir |" );
  while (<TACE>) {
    chomp;
    next if (/acedb\>/);
    next if (/\/\//);
    if (/Sequence\s+:\s+\"(\S+)\"/) {
      $id = $1;
    }
    elsif (/Species\s+\"(.+)\"/) {
      $species{BLAT_WASHU}->{$id} = $1;
    }
  }
  close TACE;

  print "Finding BLAT_NEMBASE data\n";
  open( TACE, "echo '$cmd2' | $tace $ace_dir |" );
  while (<TACE>) {
    chomp;
    next if (/acedb\>/);
    next if (/\/\//);
    if (/Sequence\s+:\s+\"(\S+)\"/) {
      $id = $1;
    }
    elsif (/Species\s+\"(.+)\"/) {
      $species{BLAT_NEMBASE}->{$id} = $1;
    }
  }
  close TACE;
  
  print "Finding BLAT_NEMATODE data\n";
  open( TACE, "echo '$cmd3' | $tace $ace_dir |" );
  while (<TACE>) {
    chomp;
    next if (/acedb\>/);
    next if (/\/\//);
    if (/Sequence\s+:\s+\"(\S+)\"/) {
      $id = $1;
    }
    elsif (/Species\s+\"(.+)\"/) {
      $species{BLAT_NEMATODE}->{$id} = $1;
    }
  }
  close TACE;

  # Not all EMBL Sequences will be in the database when this script
  # is run, because the objects for the other single-species sets
  # (briggase, remanei etc) do not make their way into autoace
  # until the merge. We therefore fill in the gaps by picking these
  # up directly. 
  my %accessors = $wormbase->all_species_accessors;
  foreach my $owb (values %accessors) {
    my $cdnadir = $owb->maskedcdna;
    my $full_sp = $owb->full_name;
    
    foreach my $file (glob("$cdnadir/*.*")) {
      open(my $fh, $file);
      while(<$fh>) {
        /^\>(\S+)/ and do {
          if (not exists $species{BLAT_NEMATODE}->{$1}) {
            $species{BLAT_NEMATODE}->{$1} = $full_sp;
          }
        }
      }
    }
  }

  return \%species;
}
  


