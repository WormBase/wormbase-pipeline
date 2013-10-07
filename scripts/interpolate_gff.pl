#!/usr/bin/env perl

use strict;
use lib $ENV{CVS_DIR};
use Modules::Physical_map;
use Wormbase;
use Getopt::Long;
use IO::File;
use File::Basename;

my (
  $store,  $test, $prep, $debug, $genes,
  $clones, $alleles, $all, $wormbase, $chromosome, $gff3,
    );

GetOptions(
  'test'         => \$test,
  'debug:s'      => \$debug,
  'genes'        => \$genes,
  'clones'       => \$clones,
  'alleles'      => \$alleles,
  'all'          => \$all,
  'store:s'      => \$store,
  'chromosome:s' => \$chromosome,
  'prepare'      => \$prep,
  'gff3'           => \$gff3,
) || die `perldoc $0`;

if ($store) {
    $wormbase = Storable::retrieve($store)
      or croak("Can't restore wormbase from $store\n");
}
else { 
  $wormbase = Wormbase->new( -debug => $debug, -test => $test ); 
}

my $log = Log_files->make_build_log($wormbase) ;# prewarning will be misused in a global way

my $acedb      = $wormbase->autoace;
my $chromdir   = $wormbase->gff_splits;
my $outdir     = $wormbase->acefiles;
##############################################################################################
#generate a new mapper based on the files (also needs to point to better gene/gene files)

my $mapping_store_file = "$acedb/logs/rev_physicals.yml";
my $mapping_log_file = "$acedb/logs/rev_physicals.log";
my $fixes_file = "$outdir/genetic_map_fixes.ace";
my $chr_prefix = $wormbase->chromosome_prefix;

# check the mappings
if ($prep) {
  unlink $mapping_store_file if -e $mapping_store_file;

  my $mapper = Physical_mapper->new( 
    undef,
    $acedb, 
    glob("$chromdir/${chr_prefix}*_gene.gff") );
  
  $mapper->check_mapping( $acedb, $fixes_file, $mapping_log_file, $log );
  $mapper->freeze($mapping_store_file);

  if (not -s $fixes_file) {
    $log->write_to("There were no fixes, so probably no need to look!\n");
  }

  $log->mail();
  exit(0);
}

my $mapper = Physical_mapper->new($mapping_store_file);
my $rev_genes = $mapper->gmap();

$log->write_to("\n\ngenerating acefiles:\n-------------------\n");

my @chromosomes;
if ($chromosome) {
  @chromosomes = ( $chromosome );
} else {
  @chromosomes = $wormbase->get_chromosome_names(-prefix => 1);
}

# specifies the Allele Methods, that should get parsed/dumped/interpolated
my @allele_methods = ('Allele',
                      'Deletion_allele',
                      'Insertion_allele',
                      'Deletion_and_Insertion_allele',
                      'Substitution_allele',
                      'Transposon_insertion');
my @gene_methods = ('gene');
my @clone_methods = ('clone_acc');

foreach my $chrom (@chromosomes) {
  my @data;
  if ($alleles or $all) {
    &dump_data( $wormbase, $chrom, @allele_methods );
    push @data, @allele_methods;
  }
  push @data, @gene_methods   if $genes or $all;
  push @data, @clone_methods  if $clones or $all;

  foreach my $f (@data) {
    my $file = "$chromdir/${chrom}_$f.gff";
    my $outfile = "$outdir/interpolated_${f}_$chrom.ace";

    my $fh = IO::File->new( $file, "r" ) || ( $log->write_to("cannot find: $file\n") && next );

    my $of = IO::File->new("> $outfile");
    $log->write_to("writing to: $outfile\n");
    
    while (<$fh>) {
      next if /\#/;
      s/\"//g;
      my @fields = split(/\t+/, $_);
      
      my ( $chr, $source, $feature) = ( $fields[0], $fields[1], $fields[2]);
      $chr =~ s/$chr_prefix//;

      my ($id, $ctag);
      if ($gff3) {              
        my ($first) = split(/;/, $fields[8]);
        ($ctag, $id) = $first =~ /^ID:(\S+):(\S+)/;
      } else {
        ($ctag, $id) = $fields[8] =~ /^(\S+)\s+(\S+)/;
      }
      
      my $class;
      if ( $source eq 'Genomic_canonical' && $feature eq 'region' ) {
        $class = 'Sequence';
      }
      elsif ( $source eq 'Allele' && $ctag eq 'Variation' ) {
        $class = 'Variation';
      }
      elsif ( $source eq 'gene' && $feature eq 'gene' ) {
        $class = 'Gene';
        # do not interpolate genes that already have a defined genetic position
        next if exists $rev_genes->{$id};
      }
      else { next }
      
      my $pos = ( $fields[3] + $fields[4] ) / 2;    # average map position
      my $aceline = $mapper->x_to_ace( $id, $pos, $chr, $class );
      
      print $of $aceline if $aceline ; # mapper returns undef if it cannot be mapped (like on the telomers)
      $log->write_to( "cannot map $class : $id (might be on a telomer) - phys.pos $chr : $pos\n"
          ) if ( !$aceline );    #--
    }
    
    close $of;
    close $fh;
  }
}
###########################################################################################
$log->mail();
exit 0;

sub dump_data {
  my ( $wormbase, $chromosome, @methods ) = @_;
  
  my $giface = $wormbase->giface;
  my $meth=join(',',@methods);

  my $cmd = "GFF_method_dump.pl -database $acedb -method $meth -dump_dir $chromdir -chromosome $chromosome -giface $giface";
  $wormbase->run_script($cmd);
}
