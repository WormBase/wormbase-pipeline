#!/software/bin/perl -w
#
# WBGene_span.pl
#
# by Anthony Rogers
#
# Creates SMapped Gene spans for Gene objects
#
# Last edited by: $Author: klh $
# Last edited on: $Date: 2015-04-23 16:33:33 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Coords_converter;
use Getopt::Long;
use Log_files;
use Storable;

my (
    $database, $species,  $test, $store, 
    $debug,  $chromosome, $no_load, $no_dump,
);

GetOptions(
    'database=s'   => \$database,
    'test'         => \$test,
    "store:s"      => \$store,
    'debug=s'      => \$debug,
    'chromosome=s' => \$chromosome,
    'species:s'    => \$species,
    'noload'       => \$no_load,
    'nodump'       => \$no_dump,
);

my $wormbase;
if ($store) {
  $wormbase = retrieve($store)
    or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug => $debug, -test => $test, -organism => $species );
}

my $log = Log_files->make_build_log($wormbase);

$database = $wormbase->autoace unless $database;
$log->write_to("Generating WBGene spans from database $database\n");

my $coords = Coords_converter->invoke( $database, undef, $wormbase );
my %worm_gene2geneID_name = $wormbase->FetchData('worm_gene2geneID_name');

my (%gene_coords, %gene_span);

my  @methods = qw(Coding_transcript Non_coding_transcript Pseudogene ncRNA 7kncRNA tRNA miRNA pre_miRNA miRNA_primary_transcript snRNA snlRNA snoRNA rRNA scRNA stRNA lincRNA asRNA piRNA Transposon_CDS Transposon_Pseudogene);
  
foreach my $method (@methods) {
  print "checking $method \n" if $debug;
  my $gff_fh = $wormbase->open_GFF_file( undef, $method, $log );
  while (<$gff_fh>) {
    my @data = split;
    if (   ( $data[1] eq 'Coding_transcript' )
           or ( $data[1] eq 'Non_coding_transcript' )
           or ( $data[1] eq 'Pseudogene' )
           or ( $data[1] eq 'curated' )
           or ( $data[1] eq 'tRNA' )
           or ( $data[1] eq 'snRNA' )
	   or ( $data[1] eq 'miRNA_mature' )
           or ( $data[1] eq 'miRNA_primary_transcript' )
	   or ( $data[1] eq 'pre_miRNA' )
           or ( $data[1] eq 'rRNA' )
           or ( $data[1] eq 'scRNA' )
           or ( $data[1] eq 'stRNA' )
           or ( $data[1] eq 'snRNA' )
           or ( $data[1] eq 'ncRNA' )
	   or ( $data[1] eq '7kncRNA' )
           or ( $data[1] eq 'snoRNA' )
           or ( $data[1] eq 'snlRNA' )
	   or ( $data[1] eq 'asRNA' )
	   or ( $data[1] eq 'piRNA' )
	   or ( $data[1] eq 'lincRNA' )
           or ( $data[1] eq 'Transposon_CDS' ) 
           or ( $data[1] eq 'Transposon_Pseudogene') ) {
      next if ( $data[2] eq 'exon'
                or $data[2] eq 'coding_exon'
                or $data[2] eq 'intron' );
      my ($gene) = $data[9] =~ /\"(\S+)\"$/;

      push @{$gene_coords{$gene}}, {
        chr     =>  $data[0],
        start   =>  $data[3],
        end     =>  $data[4],
        strand  =>  $data[6],
      };
    }
  }
  close($gff_fh);
}
    
foreach my $tran_id ( keys %gene_coords ) {
  my $WBgene = $worm_gene2geneID_name{$tran_id};
  if ( !defined $WBgene ) {
    my $cds = $tran_id;
    my $cdsregex = $wormbase->cds_regex_noend;
    $cds =~ s/($cdsregex)\.\d+/$1/;# convert transcript ID to sequence name to get WBGene ID
    $WBgene = $worm_gene2geneID_name{$cds};
    if ( !defined $WBgene ) {
      $log->write_to("*** $tran_id is not a key of worm_gene2geneID_name\n");
      next;
    }
  }
  foreach my $tran (@{$gene_coords{$tran_id}}) {
    if (not exists $gene_span{$WBgene}) {
      $gene_span{$WBgene} = {
        chr    => $tran->{chr},
        strand => $tran->{strand},
        min    => $tran->{start},
        max    => $tran->{end},
      }
    } else {
      if ($tran->{chr}    ne $gene_span{$WBgene}->{chr} or
          $tran->{strand} ne $gene_span{$WBgene}->{strand}) {
        $log->log_and_die("Multiple transcripts for $WBgene on different chromosomes/strands! Bailing\n");
      }
      $gene_span{$WBgene}->{min} = $tran->{start} if $tran->{start} < $gene_span{$WBgene}->{min};
      $gene_span{$WBgene}->{max} = $tran->{end} if $tran->{end} > $gene_span{$WBgene}->{max};
    }
  }
}


my $acefile = $wormbase->acefiles . "/WBgene_spans.ace";
open(my $ace_fh, ">$acefile" ) or $log->log_and_die->write_to("cant open output $acefile:\t$!\n");

foreach my $gene ( keys %gene_span ) {
  # write S-map details back to database
  my @coords = ($gene_span{$gene}->{strand} eq "+") 
      ? $coords->LocateSpan($gene_span{$gene}->{chr},
                            $gene_span{$gene}->{min},
                            $gene_span{$gene}->{max})
      : $coords->LocateSpan($gene_span{$gene}->{chr},
                            $gene_span{$gene}->{max},
                            $gene_span{$gene}->{min});
  
  print $ace_fh "\nSequence : \"$coords[0]\"\n";
  print $ace_fh "Gene_child $gene $coords[1] $coords[2]\n";
}
close($ace_fh);


$wormbase->load_to_database( "$database", "$acefile", "WBGene_span", $log )
    unless $no_load;

unless ($no_dump) {
  $wormbase->run_script("dump_gff_batch.pl -database ". $wormbase->autoace
			. " -methods gene -dump_dir ". $wormbase->gff_splits,$log
		       ) unless $no_dump;
}

$log->mail;

exit(0);

=pod

=Title

    WBGene_span.pl

=head1 Overview

    parses GFF files to create a gene span for each WBGene from the start of the 5' most point of all transcripts to the 3'.

=head1 Output

    writes acefile output to $database/acefiles and gff to $database/CHROMOSOMES
    if acefile is created it will be loaded in to the database

=head1 Options

    -database  default is /wormsrv2/autoace
    -test      to use the test database
    -gff       to write a GFF file of the genespans as they are created rather than dumping them after loading to database
    -no_ace    dont write an acefile

=head1 Example

    perl WBGene_span.pl -database ~wormpub/DATABASES/current_DB -gff

=cut

