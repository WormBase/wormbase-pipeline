#!/usr/bin/env perl

use strict;
use Getopt::Long;

use Storable;
use Wormbase;
use Bio::SeqIO;

my ($wormbase,$test,$debug, $store, $acefile, $noload, $database);

GetOptions( "debug=s"     => \$debug,
            "test"        => \$test,
            "acefile=s"   => \$acefile,
            "database=s"  => \$database,
            "noload"      => \$noload,
            "store:s"     => \$store);


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( 
    -debug    => $debug,
    -test     => $test,
      );
}

my $log = Log_files->make_build_log($wormbase);

$acefile = $wormbase->acefiles . "/tier3_stubs.ace" if not defined $acefile;
$database = $wormbase->autoace if not defined $database;



open(my $acefh, ">$acefile") or $log->log_and_die("Could not open $acefile for writing\n");

my %accessors = $wormbase->tier3_species_accessors;

while (my($species,$wb)=each %accessors){
  my $seqdir = $wb->sequences;
  my $gff3 = sprintf("%s/%s.gff3", $wb->sequences, $species);
  my $prot = sprintf("%s/%s.prot.fa", $wb->sequences, $species);

  $log->write_to("Doing $species from $prot and $gff3\n");

  my ($prot_seqio, $gff3_fh);
  if (-e $prot) {
    $prot_seqio = Bio::SeqIO->new(-format => 'fasta', -file => $prot);
  } elsif (-e "$prot.gz") {
    $prot_seqio = Bio::SeqIO->new(-format => 'fasta', -fh => "gunzip -c $prot.gz |");
  } else {
    $log->write_to("Could not find $prot or $prot.gz for $species\n");
    next;
  }

  while (my $seq = $prot_seqio->next_seq) {
    my $id = sprintf("%s_%s", $wb->ncbi_bioproject, $seq->id);
    printf($acefh "\nProtein : \"%s\"\n", $id);
    printf($acefh "Species \"%s\"\n", $wb->full_name);
    printf($acefh "Corresponding_CDS \"%s\"\n", $id);
    print $acefh "Live\n";

    printf($acefh "\nPeptide : \"%s\"\n", $id);
    printf($acefh "%s\n", $seq->seq);
  }

  if (-e $gff3) {
    open($gff3_fh, $gff3);
  } elsif (-e "$gff3.gz") {
    open($gff3_fh, "gunzip -c $gff3.gz |");
  } else {
    $log->log_and_die("Could not find $gff3 or $gff3.gz");
  }

  my (%geneid2name, %cdsname2geneid);

  while(<$gff3_fh>) {
    /^\#/ and next;
    chomp;

    my @l = split(/\t/, $_);

    if ($l[2] eq 'gene' or $l[2] eq 'mRNA') {
      my @attr = split(/;/, $l[8]);
      my ($name_attr) = grep { /^Name=/ } @attr;
      my ($id_attr) = grep { /^ID=/ } @attr;

      my ($name) = $name_attr =~ /Name=(\S+)/; 
      my ($id) = $id_attr =~ /ID=(\S+)/;

      next unless($name && $id); # both need to be defined to produce a valid ace-file

      my $ace_name = sprintf("%s_%s", $wb->ncbi_bioproject, $name);

      if ($l[2] eq 'gene') {
        $geneid2name{$id} = $ace_name;

        printf($acefh "\nGene : \"%s\"\n", $ace_name);
        printf($acefh "Public_name \"%s\"\n", $name);
        printf($acefh "Species \"%s\"\n", $wb->full_name);
        print $acefh "Method \"Gene\"\n";
        print $acefh "Remark \"Tier3 species gene\" Inferred_automatically \"tier3_stubs.pl\"\n";
        print $acefh "Live\n";
      } else {
        my ($parent_attr) = grep { /^Parent=/ } @attr;
        my ($parent) = $parent_attr =~ /Parent=(\S+)/;
        $cdsname2geneid{$ace_name} = $parent; 
 
        printf($acefh "\nCDS : \"%s\"\n", $ace_name);
        printf($acefh "Species \"%s\"\n", $wb->full_name);
        print $acefh "Method \"curated\"\n";
        printf($acefh "Database WormBase TierIII \"%s\"\n", $name);
        print $acefh "Remark \"Tier3 species gene\" Inferred_automatically \"tier3_stubs.pl\"\n";
        print $acefh "Coding CDS\n";
        print $acefh "From_laboratory \"HX\"\n";
      }
    }
  }

  foreach my $cds_name (sort keys %cdsname2geneid) {
    my $gene_id = $cdsname2geneid{$cds_name};
    $log->log_and_die("Could not find gene name for $cds_name\n") 
        if not exists $geneid2name{$gene_id};
    my $gene_name = $geneid2name{$gene_id};

    printf($acefh "\nCDS \"%s\"\n", $cds_name);
    printf($acefh "Gene \"%s\"\n", $gene_name);
  }
}

close($acefh) or $log->log_and_die("Could not successfully close $acefile after writing\n");

if (not $noload) {
  $log->write_to("Loading...\n|");
  $wormbase->load_to_database($database, $acefile, 'tier3_stubs', $log);
}

$log->mail();
exit();
