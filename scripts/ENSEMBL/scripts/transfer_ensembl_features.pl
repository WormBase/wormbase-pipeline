#!/usr/bin/env perl

=head1 NAME

  transfer_ensembl_features.pl

=head1 DESCRIPTION

=head1 OPTIONS

=head1 EXAMPLE

=cut

use strict;
use warnings;
use Carp;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw ( throw warning ) ;
use Getopt::Long;

my $inport      = 3306;
my $outport     = 3306;

my ($inuser, $inhost, $indbname);
my ($outuser, $outpass, $outhost, $outdbname);

my ($repeat, $dnaalign, $proteinalign, $simple, $seq_region_map, $verbose);


&GetOptions( 'inhost|sourcehost:s'       => \$inhost,
             'inuser|sourceuser:s'       => \$inuser,
             'indbname|sourcedbname:s'   => \$indbname,
             'inport|sourceport:n'       => \$inport,
             'outhost|targethost:s'      => \$outhost,
             'outuser|targetuser:s'      => \$outuser,
             'outpass|targetpass:s'      => \$outpass,
             'outdbname|targetdbname:s'  => \$outdbname,
             'outport|targetport:n'      => \$outport,
             'simple'                    => \$simple,
             'repeat'                    => \$repeat,
             'dnaalign'                  => \$dnaalign,
             'proteinalign'              => \$proteinalign,     
             'seqregionmap'              => \$seq_region_map,
             'verbose'                   => \$verbose,
    );

die "You must give details of the source database with -indbhost, -indbport, -indbuser and -indbname\n"
    if not defined $inhost or not defined $inport or not defined $inuser or not defined $indbname;
die "You must give details of the target database with -outdbhost, -outdbport, -outuser, -outpass, and -outdbname\n"
    if not defined $outhost or not defined $outport or not defined $outuser or not defined $outpass or not defined $outdbname;

my %seq_region_map;
if (defined $seq_region_map) {
  open(my $fh, $seq_region_map) or die "Could not open $seq_region_map for reading\n";
  while(<$fh>) {
    /^(\S+)\s+(\S+)/ and do {
      $seq_region_map{$1} = $2;
    }
  }
}


my $sourcedb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $inhost,
                                        -user   => $inuser,
                                        -port   => $inport,
                                        -dbname => $indbname );

my $targetdb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $outhost,
                                        -user   => $outuser,
                                        -pass   => $outpass,
                                        -port   => $outport,
                                        -dbname => $outdbname );

my @s_slices = sort { $b->length <=> $a->length } @{$sourcedb->get_SliceAdaptor->fetch_all('toplevel')};


my (%target_analyses);

foreach my $s_sl (@s_slices) {
  my $name = $s_sl->seq_region_name;
  if (exists $seq_region_map{$name}) {
    $name = $seq_region_map{$name};
  }

  $verbose and print STDERR "Processing " . $s_sl->name . "\n";
  
  my $t_sl = $targetdb->get_SliceAdaptor->fetch_by_region('toplevel', $name);
  die "Could not fetch slice $name from target database\n" if not defined $t_sl;

  if ($repeat) {
    $verbose and print STDERR "   Transferring repeats...\n";
    &transfer($t_sl, $targetdb->get_RepeatFeatureAdaptor, $s_sl->get_all_RepeatFeatures);
  }

  if ($simple) {
    $verbose and print STDERR "   Transferring simple features...\n";
    &transfer($t_sl, $targetdb->get_SimpleFeatureAdaptor, $s_sl->get_all_SimpleFeatures);
  }

  if ($dnaalign) {
    $verbose and print STDERR "   Transferring dnaalign features...\n";
    &transfer($t_sl, $targetdb->get_DnaAlignFeatureAdaptor, $s_sl->get_all_DnaAlignFeatures);
  }

  if ($proteinalign) {
    $verbose and print STDERR "   Transferring proteinalign features...\n";
    &transfer($t_sl, $targetdb->get_ProteinAlignFeatureAdaptor, $s_sl->get_all_ProteinAlignFeatures);
  }
}


###########################################
sub transfer {
  my ($t_sl, $fadap, $feats) = @_;

  my @to_store;

  foreach my $f (@$feats) {
    $f->slice($t_sl);
    my $logic = $f->analysis->logic_name;
    if (not exists $target_analyses{$logic}) {
      my $ana = $targetdb->get_AnalysisAdaptor->fetch_by_logic_name($logic);
      if (not defined $ana) {
        die "Could not find analysis '$logic' in target database\n";
      } 
      $target_analyses{$logic} = $ana;
    }
    $f->analysis($target_analyses{$logic});

    push @to_store, $f;
  }

  $fadap->store(@to_store);
}

  

