#!/usr/local/bin/perl5.8.0 -w
#
# processGFF.pl
#
# by Anthony Rogers
#
# handles post processing of GFF files
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2006-02-14 09:50:41 $


use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use strict;
use Storable;
use Log_files;

my ($help, $debug, $test, $quicktest, $database, $store );
my ($utr, $clone_acc, $gene_acc, $nematode);
my ($input, $output);
GetOptions (
	    "help"          => \$help,
	    "debug=s"       => \$debug,
	    "test"          => \$test,
	    "store:s"       => \$store,
	    "quicktest"     => \$quicktest,
	    "database:s"    => \$database,

	    "utr"           => \$utr,
	    "clone_acc"     => \$clone_acc,
	    "gene_acc"      => \$gene_acc,
	    "nematode"      => \$nematode,

	    "input:s"       => \$input,
	    "output:s"      => \$output,
	   );

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}

my $log = Log_files->make_build_log($wormbase);

$database = $wormbase->autoace unless $database;
my @chroms = qw( I II III IV V X MtDNA);
@chroms = qw(III) if $quicktest;

if( $input and $output ){
  &GFF_with_UTR($input, $output)               if $utr;
  &GFF_clones_with_accessions($input, $output) if $clone_acc;
  &GFF_genes_with_accessions($input, $output)  if $gene_acc;
  $wormbase->run_script('add_species_to_BLAT_GFF.pl', $log) if $nematode;
}
else {
  if( $utr ) {
    foreach my $chrom ( @chroms ) {
      my $in = $wormbase->gff_splits."/CHROMOSOME_${chrom}_UTR.gff";
      my $out = $wormbase->gff_splits."/CHROMOSOME_${chrom}_UTR_CDS.gff";
      &GFF_with_UTR( $in, $out );
    }
  }
  if( $clone_acc ) {
    my $db = Ace->connect(-path=>$database) || do { print "Connection failure to $database: ",Ace->error; die();};
    foreach my $chrom ( @chroms ) {
      my $in = $wormbase->gff_splits."/CHROMOSOME_${chrom}_Genomic_canonical.gff";
      my $out = $wormbase->gff_splits."/CHROMOSOME_${chrom}_clone_acc.gff";
      &GFF_clones_with_accessions( $in, $out, $db );
    }
    $db->close;
  }
  if( $gene_acc ) {
    my $db = Ace->connect(-path=>$database) || do { print "Connection failure to $database: ",Ace->error; die();};
    foreach my $chrom ( @chroms ) {
      my $in = $wormbase->gff_splits."/CHROMOSOME_${chrom}_CDS.gff";
      my $out = $wormbase->gff_splits."/CHROMOSOME_${chrom}_CDS_acc.gff";
      &GFF_genes_with_accessions( $in, $out, $db );
    }
    $db->close;
  }
}

$log->mail;

exit(0);

sub GFF_with_UTR{
  my $utr = shift;
  my $utr_cds = shift;

  $log->write_to("writing UTR GFF file $utr_cds\n");
  open( UTR, "<$utr" );
  open( UTR_CDS, ">$utr_cds");
  while (<UTR>) {
    if (/^\#/){
      print UTR_CDS $_;
    }
    else{
      (/_UTR:(\S+)\"/);
      print UTR_CDS "$`" . "_UTR:$1\" CDS=\"$1\"\n";
    }
  }
  close UTR;
  close UTR_CDS;
  
  system ("mv -f $utr_cds $utr");
}


sub GFF_clones_with_accessions{
  my $infile   = shift;
  my $outfile = shift;
  my $db = shift;
  my $wormdb = "$database";

  unless( $db ){
    $db = Ace->connect(-path=>$database) || do { print "Connection failure to $wormdb: ",Ace->error; die();};
  }

  open (GFF, "<$infile")  || die "Can't open GFF file: $infile\n\n";
  open (OUT, ">$outfile") || die "Can't write to $outfile\n";

  $log->write_to("writing clone accession GFF file $outfile\n");
  while (<GFF>) {

    next if (/^\#/);
    chomp;

    my @gff = split (/\t/,$_);
    next unless ( $gff[1] eq "Genomic_canonical" ); #Genome_canonical dump includes assembly tags etc
    my ($gen_can) = $gff[8] =~ /Sequence \"(\S+)\"/; 
    my $obj = $db->fetch(Sequence=>$gen_can);
    if (!defined ($obj)) {
      print "Could not fetch sequence '$gen_can'\n";
      next;
    }

    # Grab accession and sequence version info from the database
    my $acc;
    my $sv;

    if($obj->at('DB_info.Database.EMBL.NDB_AC')){
      ($acc) = $obj->at('DB_info.Database.EMBL.NDB_AC');
      ($sv) = $obj->at('DB_info.Database.EMBL.NDB_SV');
    }
    elsif($obj->at('DB_info.Database.GenBank.NDB_AC')){
      ($acc) = $obj->at('DB_info.Database.GenBank.NDB_AC');
      ($sv) = $obj->at('DB_info.Database.GenBank.NDB_SV');
    }

    # now just want the numerical suffix of sequence version field
    $sv =~ s/.*\.//;

     print OUT "$_ acc=$acc ver=$sv\n";
     print "$_ acc=$acc ver=$sv\n" if $wormbase->debug;

    $obj->DESTROY();

  }
  close(GFF);
  close(OUT);
}

########################################

sub GFF_genes_with_accessions{
  my $infile  = shift;
  my $outfile = shift;

  $log->write_to("writing Genes with acc GFF file $outfile\n");
  my $db = Ace->connect(-path=>$database) || do { print "Connection failure to $database: ",Ace->error; die();};

  open (GFF, "<$infile") || die "Can't open GFF file: $infile\n\n";
  open (OUT, ">$outfile") || die "Can't write to $outfile\n";
  while (<GFF>) {

    next if (/^\#/);    
    chomp;

    my @gff = split (/\t/,$_);
    next if ($gff[1] ne "curated");
    (my $gen_can) = $gff[8] =~ /CDS \"(\S+)\"/; 
    my $obj = $db->fetch(CDS=>$gen_can);
    if (!defined ($obj)) {
      print "Could not fetch CDS '$gen_can'\n";
      next;
    }

    my $acc = $obj->Corresponding_protein(1);
    $acc =~ s/WP\://g;

    if ($acc ne "") {
      $acc =~ s/WP\://g;
      print OUT "$_ wp_acc=$acc\n";
    }
    else {  # could be a tRNA
      print OUT "$_\n";
    }
    undef ($acc); 

    $obj->DESTROY();

  }
  close (GFF);
  close (OUT);
}
