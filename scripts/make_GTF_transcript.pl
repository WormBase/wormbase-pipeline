#!/software/bin/perl -w
#
# make_GTF_transcript.pl
#
# Ceates a GTF file for reading into Cufflinks of the transcripts of all expressed genes, pseudogene, non-coding transcripts, transposons and their isoforms
#
#
# by Gary Williams
#
# Last updated by: $Author: gw3 $                      
# Last updated on: $Date: 2015-05-14 09:38:09 $        

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use IO::Handle;
use Cwd;
use Ace;
use Sequence_extract;
use Coords_converter;


######################################
# variables and command-line options #
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $species, $noprefix);
my ($output, $dbdir);

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "test"        => \$test,
	    "verbose"     => \$verbose,
	    "store:s"     => \$store,
	    "output=s"    => \$output,
	    "database=s"  => \$dbdir,
	    "species=s"   => \$species,
	    "noprefix"    => \$noprefix # remove the CHROMSOME_ prefix from the chromosome name
	    );

#$test = 0;
#$debug = "gw3";

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species,
			     );
}



# Display help if required
&usage("Help") if ($help);

$dbdir          = $wormbase->test ? $wormbase->database("current") : $wormbase->autoace unless $dbdir; # Database path
my $common_dir      = "$dbdir/COMMON_DATA";       # GFF directory
my $gffdir          = "$dbdir/GFF_SPLITS";
my $prefix          = $wormbase->chromosome_prefix;

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


if (!defined($output)){
  die "No -output file\n";
}

my %dummy;
my %worm_gene2geneID_name = $wormbase->FetchData('worm_gene2geneID_name', \%dummy, $common_dir);

my $old_output_file_size = -s $output;


open (OUTDAT, ">$output") || die "Can't open $output\n";

my @file_types = qw(
	Coding_transcript
	Non_coding_transcript
	miRNA_primary_transcript
        pre_miRNA
        miRNA
        7kncRNA
	ncRNA
	rRNA
        piRNA
        asRNA
        lincRNA
	snoRNA
	scRNA
	snRNA
	stRNA
	tRNA
        circRNA
);


my @files;
if ($wormbase->assembly_type eq 'contig') {
  @files = map { "$gffdir/$_.gff" } @file_types;
} else {
  my @chroms = $wormbase->get_chromosome_names(-prefix => 1, -mito => 1);
  foreach my $chrom (@chroms) {
    foreach my $pattern (map { "$gffdir/${chrom}_${_}.gff" } @file_types) {
      push @files, glob($pattern);
    }
  }
}


foreach my $file (@files) {
  open (IN, "<$file") || $log->log_and_die("Can't open $file\n");
  my $count=0;
  while (my $line = <IN>) {
    $count++;
    chomp $line;
    next if ($line =~ /^#/);
    next if ($line =~ /^\s*$/);

    my @line = split /\t+/, $line;
    next if $line[2] eq 'intron';

    if ($line[8] =~ /Transcript\s+(\S+)/) {
      my $transcript_id = $1;
      $transcript_id =~ s/"//g;

      if ($noprefix) {
	$line[0] =~ s/${prefix}(.+)/$1/;
      }

      if ($line[2] ne 'exon') {
	$line[2] = 'transcript';
      }
      $line[1] = 'WormBase';
      my $cds_regex = $wormbase->cds_regex_noend;
      my ($sequence_name) = ($transcript_id =~ /($cds_regex)/);

#      if ($transcript_id eq "F53G12.15") {print "linecount=$count, have F53G12.15, file = $file line=$line\n";}
#      if ($transcript_id eq "Y65B4BL.10") {print "linecount=$count, have Y65B4BL.10, file = $file line=$line\n";}

      if (!defined $sequence_name) {$sequence_name = $transcript_id}# the tRNAs (e.g. C06G1.t2) are not changed correctly
      my $gene_id = $worm_gene2geneID_name{$sequence_name};
      if (!defined $gene_id) {$log->log_and_die ("Can't find gene_id $sequence_name for: $line\n")}
      $line[8] = "gene_id \"$gene_id\"; transcript_id \"$transcript_id\"";
      $line = join "\t", @line;
      print OUTDAT "$line\n";
    }
  }
  close (IN);
}

#print "Now doing Pseudogenes\n";
# now do the Pseudogenes - these don't have Transcripts
    
@files = glob("$gffdir/*Pseudogene.gff");

foreach my $file (@files) {
  open (IN, "<$file") || $log->log_and_die("Can't open $file\n");
  while (my $line = <IN>) {
    chomp $line;
    next if ($line =~ /^#/);
    next if ($line =~ /^\s*$/);

    my @line = split /\t+/, $line;
    if ($line[8] =~ /Pseudogene\s+(\S+)/) {
      if ($line[2] eq 'intron') {next}
      if ($line[1] eq 'Transposon_Pseudogene') {next} # these don't have gene_ids
      my $transcript_id = $1;
      $transcript_id =~ s/"//g;

      if ($noprefix) {
	$line[0] =~ s/${prefix}(.+)/$1/;
      }

      if ($line[2] ne 'exon') {
	$line[2] = 'transcript';
      }
      $line[1] = 'WormBase';
      my $cds_regex = $wormbase->cds_regex_noend;
      my ($sequence_name) = ($transcript_id =~ /($cds_regex)/);
      if (!defined $sequence_name) {$sequence_name = $transcript_id}# the tRNAs (e.g. C06G1.t2) are not changed correctly
      my $gene_id = $worm_gene2geneID_name{$sequence_name};
      if (!defined $gene_id) {$log->log_and_die ("Can't find Pseudogene gene_id $sequence_name for: $line\n")}
      $line[8] = "gene_id \"$gene_id\"; transcript_id \"$transcript_id\"";
      $line = join "\t", @line;
      print OUTDAT "$line\n";

    }
  }
  close (IN);
}

if (defined $old_output_file_size) {
  my $output_file_size = -s $output;
  if ($old_output_file_size < $output_file_size * 0.9 || $old_output_file_size > $output_file_size * 1.1) {
    $log->error("WARNING: old output file size: $old_output_file_size, new output file size: $output_file_size\n");
  }
}




close(OUTDAT);

$log->mail();
print "Finished" if ($verbose);
exit(0);




###############################
# Prints help and disappears  #
###############################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}


