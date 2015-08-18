#!/usr/bin/env perl

use strict;
use lib $ENV{CVS_DIR};
use Modules::Physical_map;
use Wormbase;
use Getopt::Long;
use IO::File;
use File::Basename;

my (
  $store,  $test, $prep_fix, $debug, $genes, $prep_nofix, $fix_acefile,
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
  'preparefix'   => \$prep_fix,
  'preparenofix' => \$prep_nofix,
  'fixacefile=s' => \$fix_acefile,
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
my $chr_prefix = $wormbase->chromosome_prefix;
$fix_acefile = "$outdir/genetic_map_fixes.ace" if not defined $fix_acefile;

# check the mappings
if ($prep_fix or $prep_nofix) {
  unlink $mapping_store_file if -e $mapping_store_file;

  my $mapper = Physical_mapper->new( 
    undef,
    $acedb, 
    glob("$chromdir/${chr_prefix}*_gene.gff") );

  my $errors = $mapper->check_and_fix_mapping( 
    $acedb, 
    $fix_acefile, 
    $mapping_log_file, 
    $prep_fix,
    $log );

  $mapper->freeze($mapping_store_file);
  
  if ($prep_fix) {
    $log->write_to("There were $errors errors after attempting to fix. Check out $mapping_log_file if this is > 0!.\n");
  } else {
    $log->write_to("There were $errors errors, did not attempt to fix. If >0, check out $mapping_log_file to resolve manually.\n");
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
my @allele_methods = ('Deletion_allele',
                      'Insertion_allele',
                      'Deletion_and_Insertion_allele',
                      'Substitution_allele',
                      'Transposon_insertion');
my @gene_methods = ('gene');
my @clone_methods = ('clone_acc');

my %check_results;
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
        if (exists $rev_genes->{$id}) {
          $check_results{$id} = {
            obj => $id,
            chr => $chr,
            ppos => ($fields[3] + $fields[4]) / 2,
            gpos => $rev_genes->{$id}->[0],
            interpolated => 0,
          };
          next;
        }
      }
      else { 
        next; 
      }
      
      my $pos = ( $fields[3] + $fields[4] ) / 2;    # average map position
      my $map_pos = $mapper->map($pos, $chr);
      if (defined $map_pos) {
        my $aceline = sprintf("$class : \"$id\"\nInterpolated_map_position \"$chr\" $map_pos\n\n");
        print $of $aceline if $aceline ; 
        $check_results{$id} = {
          obj => $id,
          chr => $chr,
          ppos => $pos,
          gpos => $map_pos,
          interpolated => 1,
        };
      } else {
        $log->write_to( "cannot map $class : $id (might be on a telomer) - phys.pos $chr : $pos\n" );
      }
    }
    
    close $of;
    close $fh;
  }
}

# check results; the interpolated genes should all have consistent genetic and physical positions
foreach my $chr ($wormbase->get_chromosome_names(-prefix => 0, -mito => 0)) {
  $log->write_to("Checking interpolations for $chr\n");
  my @entries = grep { $_->{chr} eq $chr } values %check_results;
  @entries = sort { $a->{ppos} <=> $b->{ppos} } @entries;
  for(my $i=1; $i < @entries; $i++) {
    if ($entries[$i]->{gpos} < $entries[$i-1]->{gpos}) {
      $log->write_to(sprintf("Inconsistent: %s %s (%s %s %d) %s (%s %s %d)\n", 
                             $chr, 
                             $entries[$i-1]->{obj}, 
                             $entries[$i-1]->{ppos}, 
                             $entries[$i-1]->{gpos}, 
                             $entries[$i-1]->{interpolated}, 
                             $entries[$i]->{obj}, 
                             $entries[$i]->{ppos}, 
                             $entries[$i]->{gpos}, 
                             $entries[$i]->{interpolated}));
    }
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
