#!/usr/bin/env perl
#
# make_UTR_GFF.pl
#
# Creates three_prime_UTR, five_prime_UTR and coding_exon lines
# from the Coding_transcript and curated GFF split files
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-12-02 16:14:54 $


use strict;
use warnings;
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;

my ( $debug, $test, $store, $wormbase, $species, $verbose, $gff3);

GetOptions(
  "debug=s"      => \$debug,
  "test"         => \$test,
  "store:s"      => \$store,
  "species:s"    => \$species,
  "verbose"      => \$verbose,
  "gff3"         => \$gff3,
    );


if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n");
} else { 
  $wormbase = Wormbase->new( -debug => $debug, 
                             -test => $test,
                             -organism => $species);
}

my $log = Log_files->make_build_log($wormbase);

my @bits_of_work;
if ($wormbase->assembly_type eq 'contig') {
  push @bits_of_work, {
    curated => $wormbase->GFF_file_name(undef, 'curated'),
    coding_transcript => $wormbase->GFF_file_name(undef, 'Coding_transcript'),
    output => ($gff3) ? $wormbase->GFF3_file_name(undef, 'UTR') : $wormbase->GFF_file_name(undef, 'UTR'),
  };
} else {
  foreach my $chr ($wormbase->get_chromosome_names(-prefix => 1, -mito => 1)) {
    push @bits_of_work, {
      curated => $wormbase->GFF_file_name($chr, 'curated'),
      coding_transcript => $wormbase->GFF_file_name($chr, 'Coding_transcript'),
      output => ($gff3) ? $wormbase->GFF3_file_name($chr, 'UTR') :$wormbase->GFF_file_name($chr, 'UTR'),
    };
  }
}

foreach my $bit (@bits_of_work) {
  my $curated = $bit->{curated};
  my $coding_trans = $bit->{coding_transcript};
  my $outfile = $bit->{output};

  my ( %cds, %transcripts);

  open(my $outfh, ">$outfile") 
      or $log->log_and_die("Could not open $outfile for writing\n");

  $log->write_to("Processing $curated\n") if $verbose;

  open(my $cfh, $curated) or $log->log_and_die("Could not open $curated for reading\n");
  while(<$cfh>) {
    /^\#/ and next;
    chomp;
    my @l = split(/\t/, $_);
    next if $l[1] ne 'curated';

    my ($cds) = $l[8] =~ /CDS\s+\"(\S+)\"/;

    if ($l[2] eq 'CDS') {
      $cds{$cds} =  {
        chr    => $l[0],
        start  => $l[3],
        end    => $l[4],
        strand => $l[6],
      };
    } elsif ($l[2] eq 'coding_exon') {
      push @{$cds{$cds}->{coding_exons}}, {
        start  => $l[3],
        end    => $l[4],
        phase  => $l[7],
      };
    }
  }

  $log->write_to("Processing $coding_trans\n") if $verbose;

  open(my $codfh, $coding_trans) or $log->log_and_die("Could not open $coding_trans for reading\n");
  while(<$codfh>) {
    /^\#/ and next;
    chomp;
    my @l = split(/\t/, $_);
    next if $l[1] ne 'Coding_transcript' or $l[2] ne 'exon';

    # each exon will end up as one of coding_exon, three_prime_utr, five_prime_utr, or a comination
    my ($trans) = $l[8] =~ /Transcript\s+\"(\S+)\"/;
    my $cds = &short_name($trans);

    $transcripts{$cds}->{$trans} = 1;
    
    $log->and_die("Could not find correspondind CDS for $trans\n")
        if not exists $cds{$cds};

    my @utrs;
    if ($l[3] < $cds{$cds}->{start}) {
      # UTR left
      push @utrs, {
        type => ($l[6] eq '+') ? 'five' : 'three',
        start =>  $l[3],
        end => ($l[4] < $cds{$cds}->{start}) ? $l[4] : $cds{$cds}->{start} - 1,
      };
    }
    if ($l[4] > $cds{$cds}->{end}) {
      push @utrs, {
        type => ($l[6] eq '+') ? 'three' : 'five',
        end => $l[4],
        start => ($l[3] > $cds{$cds}->{end}) ? $l[3] : $cds{$cds}->{end} + 1,
      };
    }

    foreach my $utr (@utrs) {
      print $outfh join("\t", 
                        $cds{$cds}->{chr}, 
                        "Coding_transcript",
                        "$utr->{type}_prime_UTR",
                        $utr->{start},
                        $utr->{end},
                        ".",
                        $cds{$cds}->{strand},
                        ".",
                        ($gff3) ? "Parent=Transcript:$trans" : "Transcript \"$trans\""), "\n";
    }
  }

  # Do not need the coding_exons for GFF3
  if (not $gff3) {
    foreach my $cds (sort keys %cds) {
      foreach my $trans (sort keys %{$transcripts{$cds}}) {
        foreach my $cex (@{$cds{$cds}->{coding_exons}}) {
          print $outfh join("\t", 
                            $cds{$cds}->{chr},
                            "Coding_transcript",
                            "coding_exon", 
                            $cex->{start},
                            $cex->{end},
                            ".",
                            $cds{$cds}->{strand},
                            #$cex->{phase},
                            ".",
                            "Transcript \"$trans\" ; CDS \"$cds\""), "\n";
        }
      }
    }
  }
  close($outfh) or $log->log_and_die("Could not close $outfile after writing\n");

  $wormbase->check_file($outfile, 
                        $log,
                        lines => ['^##', 
                                  "^\\S+\\s+Coding_transcript\\s+(three_prime_UTR|coding_exon|five_prime_UTR)\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+\\s+\\S+"],
                        gff => 1,
      );  
}


$log->mail();
exit(0);

##########################


# get CDS name (lifted from original version)
sub short_name {
  my ($name) = @_;
  my $cds_regex = $wormbase->cds_regex_noend;
  my ($cdsname) = $name =~ /($cds_regex)/;
  if (! defined $cdsname) {
    $log->log_and_die("There is a problem with extracting the CDS name from the Transcript name: $name\n");
  }
  return $cdsname;
}

1;
