#!/usr/local/bin/perl5.8.0 -w
#
# agp2dna.pl
#
# by Dan Lawson
#
# Reconstructs chromosome DNA consensus file from the agp file.
# Each clone segment of sequence is checked from the EMBL
# entry and the acedb derived DNA file.
#
# Usage : agp2dna.pl [-options]
#
# Last edited by: $Author: klh $
# Last edited on: $Date: 2012-06-26 10:43:18 $

use strict;
use lib $ENV{'CVS_DIR'};

use Bio::PrimarySeq;
use Bio::SeqIO;

use Wormbase;
use Getopt::Long;
use Carp;
use Storable;
use Log_files;

use vars qw ($seq_len $sv_acc $sv_ver );


########################################
# command-line options                 #
########################################

my $debug;              # Verbose debug mode
my $help;               # Help/Usage page
my $chrom;              # Single chromosome mode - provide chromosome name
my $test;               # use test environment in ~wormpub/TEST_BUILD
my $quicktest;          # same as -test but only uses one chromosome
my $verbose;
my $store;
my $wormbase;
my $genome_agp_file;
my $output_seq_file;
my $species;

GetOptions (
  "debug=s"      => \$debug,
  "help"         => \$help,
  "chrom=s"      => \$chrom,
  "test"         => \$test,
  "verbose"      => \$verbose,
  "store:s"      => \$store,
  "agpfile:s"    => \$genome_agp_file,
  "outseqfile:s" => \$output_seq_file,
  "species:s"    => \$species,
    );
($test = 1) if ($quicktest);


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
      );
}

$species = $wormbase->species;

# establish log file.
my $log = Log_files->make_build_log($wormbase);


##############################
# Script variables (run)     #
##############################

# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir   = $wormbase->basedir;
my $dnadir    = $wormbase->chromosomes;

my (%agp, @chromosomes);

if ($genome_agp_file) {
  &read_agp_file($genome_agp_file, \%agp);
} else {
  my $agpdir    = $wormbase->autoace."/yellow_brick_road";
  foreach my $file (glob("$agpdir/*.agp")) {
    &read_agp_file($file, \%agp);
  }
}

if ($chrom) {
  push @chromosomes, $chrom;
} else {
  @chromosomes = $wormbase->get_chromosome_names(-prefix => 1, -mito => 0 );
}

########################################
# Loop over each chromosome            # 
########################################
my @EMBL_seqs;

foreach my $chromosome (@chromosomes) {
    
  my ($wormbase_seq, $EMBL_con_seq) = ("", "");

  if (not exists $agp{$chromosome}) {
    $log->log_and_die("Could not find AGP entries for $chromosome\n");
  }
  
  open (DNA, "<$dnadir/${chromosome}.dna") or $log->log_and_die("Could not open DNA file for $chromosome\n");
  while (<DNA>) {
    chomp;
    next if (/^>/);
    $wormbase_seq .= uc($_);
  }
  close DNA;
  

  foreach my $agp_line (@{$agp{$chromosome}}) {
    my @f = split (/\s+/, $agp_line);

    # do block if not a gap
    unless ($f[4] eq "N") {
      my ($acc,$sv) = split (/\./,$f[5]);

      my ($clone_start, $clone_end) = ($f[6], $f[7]);
      my ($chrom_start, $chrom_end) = ($f[1], $f[2]);

      my $chrom_seg_len = $chrom_end - $chrom_start + 1;
      my $clone_seg_len = $clone_end - $clone_start + 1;

      my ($EMBL_seq, $EMBL_sv) = &sequence_fetch($acc);

      if (not $EMBL_seq or not $EMBL_sv) {
        $log->error("ERROR: Could not fetch the sequence for $acc [ACEDB:$sv <=> EMBL:$EMBL_sv\n");
      }

      if ($f[8] eq '-') {
        $EMBL_seq = &rev_comp($EMBL_seq);
      }      

      my $EMBL_slice = substr($EMBL_seq,$clone_start - 1, $clone_seg_len);
            
      # check against WormBase sequence
      my $wormbase_slice = uc (substr($wormbase_seq, $chrom_start - 1, $chrom_seg_len));
      
      # print log line
      $log->write_to(sprintf("[%8d : %8d] %8s => Adding %6d bp from position %8d to %8d from version $EMBL_sv\n", 
                             $chrom_start,
                             $chrom_end,
                             $acc,
                             $clone_seg_len,
                             $clone_start,
                             $clone_end));
      
      # add to consensus sequence
      $EMBL_con_seq .= $EMBL_slice;
      
      # Sequence_version difference
      if ($sv =! $EMBL_sv) {
	$log->error("ERROR: Discrepent sequence version for $acc [ACEDB:$sv <=> EMBL:$EMBL_sv\n");
      }
            
      # Sequence difference
      if ($EMBL_slice ne uc$wormbase_slice) {
        my $message = "ERROR: you are not adding the same sequence for $acc\n";

        my ($count_a,$count_c,$count_g,$count_t,$count_n) = $wormbase->DNA_string_composition($wormbase_slice);
        $message .= "  WormBase [$acc] : A=$count_a C=$count_c G=$count_g T=$count_t N=$count_n\n";

        ($count_a,$count_c,$count_g,$count_t,$count_n) = $wormbase->DNA_string_composition($EMBL_slice);
        $message .= "  EMBL     [$acc] : A=$count_a C=$count_c G=$count_g T=$count_t N=$count_n\n";
	$log->error($message);
      }
    }
    # else insert gap 
    else {
      # print log line
      $log->write_to("[$f[1] : $f[2]] Adding $f[5] bp of padding characters (sequence gap)\n");
      
      # add to consensus sequence
      $EMBL_con_seq .= 'N' x $f[5];
    }
  }
  
  push @EMBL_seqs, Bio::PrimarySeq->new(-id => "$chromosome", 
                                        -seq => $EMBL_con_seq);
}

if ($output_seq_file) {
  my $seqio = Bio::SeqIO->new(-format => 'fasta',
                              -file => ">$output_seq_file");
  foreach my $seqo (@EMBL_seqs) {
    $seqio->write_seq($seqo);
  }
  $seqio->close;
}


########################################
# hasta luego                          #
########################################
$log->mail;
print "Finished.\n" if ($verbose);
exit(0);


#####################
sub sequence_fetch {
  my ($acc) = @_;

  open(ENTRY, "wget -O - 'http://www.ebi.ac.uk/ena/data/view/$acc&display=fasta' |")
      or $log->log_and_die("Could not open EBI sequence fetch command for $acc\n");

  my ($seq, $ver);

  while (<ENTRY>) {
    chomp;

    /^\>ENA\|\S+\|(\S+)\.(\d+)/ and do { 
      $ver = $2;
      next;
    };

    /^(\S+)$/ and $seq .= uc($1);
  }
  close(ENTRY) or $log->log_and_die("Could not successfully close EBI sequence command for $acc\n");

  return ($seq, $ver);

}


####################
sub rev_comp {
  my ($dna) = @_;

  $dna = scalar(reverse($dna));
  $dna =~ tr/ACGT/TGCA/;
  
  return $dna;
}


####################
sub read_agp_file {
  my ($file, $agp_ref) = @_;
  
  open(F, $file) or $log->log_and_die("Could not open AGP $file for reading\n");
  while(<F>) {
    next if /^\#/;

    my ($chrom) = /^(\S+)\s+/;
    if ($species eq 'elegans') {
      $chrom = "CHROMOSOME_${chrom}";
    }

    push @{$agp_ref->{$chrom}}, $_;
  }
}

__END__

