#!/usr/local/bin/perl5.8.0 -w

# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-03-01 17:51:52 $

package Geneace;

use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Ace;
use Wormbase;
use Coords_converter;
use strict;

my $def_dir = "/wormsrv1/geneace/wquery/";
my $curr_db = "/nfs/disk100/wormpub/DATABASES/current_DB";
my $geneace_dir = "/wormsrv1/geneace";
my $tace = &tace;

sub init {
  my $class = shift;
  my $this={};
  bless($this, $class);
}

sub geneace {
  my $this = shift;
  my $geneace_dir = "/wormsrv1/geneace";
  return $geneace_dir;
}

sub curr_db {
  my $this = shift;
  my $curr_db_dir = "/nfs/disk100/wormpub/DATABASES/current_DB";
  return $curr_db_dir;
}

sub test_geneace {
  my $this = shift;
  my $test_dir = "/nfs/disk100/wormpub/DATABASES/TEST_DBs/CK1TEST";
  return $test_dir;
}

sub parse_inferred_multi_pt_obj {
  my ($this, $version) = @_;
  my ($multi_obj, %locus_multi, $locus);
  my $multi_file = glob("/wormsrv1/geneace/JAH_DATA/MULTI_PT_INFERRED/inferred_multi_pt_obj_WS$version");
  my @multi_file = `cat $multi_file`;

  foreach (@multi_file){
    chomp;
    if ($_ =~ /^Multi_pt_data : (\d+)/){
      $multi_obj = $1;
    }
    if ($_ =~ /^Combined Locus \"(.+)\" 1 Locus \"(.+)\" 1 Locus \"(.+)\"/){
      $locus = $2;
      $locus_multi{$locus} = $multi_obj;
    }
  }
  return %locus_multi;
}


sub other_name {
  my ($this, $db, $option) = @_; # $db is db handle
  my %main_other;
  push( my @result, $db->find("Find Gene_name * where Other_name_for AND !(Cb-* OR Cr-*)") );
  if ($option eq "hash"){
    foreach(@result){
      push(@{$main_other{$_ -> Other_name_for(1)}}, $_);
    }
    return %main_other;
  }
  if ($option eq "array"){
    return @result;
  }
  if (!$option){
    print "You need to specify datatype to resturn: array or hash\n";
  }
}
sub cgc_name_is_also_other_name {
  my ($this, $db) = @_; # $db is db handle
  push( my @exceptions, $db->find("Find Gene_name * where CGC_name_for; Other_name_for") );
  return @exceptions;
}

sub loci_have_multi_pt {
  my ($this, $db) = @_;
  push( my @loci_has_multi, $db->find("Find Locus * where Multi_point") );

  my %loci_2_multi;
  foreach (@loci_has_multi){
    push(@{$loci_2_multi{$_}}, $_ -> Multi_point(1) );
  }
  return %loci_2_multi;
}

sub clone_to_lab {
  my $this = shift;
  my %clone_lab;

  my $clone_to_lab=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/clone_to_lab.def" quit
EOF
  open (FH, "echo '$clone_to_lab' | $tace $curr_db |") || die "Couldn't access $curr_db\n";
  while (<FH>){
    chomp $_;
    if ($_ =~ /\"(.+)\"\s+\"(.+)\"/){
      $clone_lab{$1} = $2;
    }
  }
  return %clone_lab;
}

sub upload_database {
  my($this, $db_dir, $command, $tsuser, $log) = @_;
  open (Load_GA,"| $tace -tsuser \"$tsuser\" $db_dir >> $log") || die "Failed to upload to Geneace";
  print Load_GA $command;
  close Load_GA;
}

#sub other_name_exceptions {
#  # list of loci which is used as other_name and CGC_name
#  my $this = shift;
#  my $db = &database;
#  push(my @exceptions, $db->
#     }

sub array_comp {

  my ($this, $ary1_ref, $ary2_ref, $option)=@_;
  my(@union, @isect, @diff, %union, %isect, %count, $e);
  @union=@isect=@diff=();
  %union=%isect=();
  %count=();
  foreach $e(@$ary1_ref, @$ary2_ref){
    $count{$e}++;
  }
  foreach $e (keys %count){
    push (@union, $e);
    if ($count{$e}==2){push @isect, $e;}
    else {push @diff, $e;}
  } 
  return \@diff, \@isect, \@union if !$option;
  return @diff if $option eq "diff";
  return @isect if $option eq "same";
  return @union if $option eq "union";
}

sub get_clone_chrom_coords {

  my %clone_info;
  my $gff_version = get_wormbase_version() - 1;  # current_DB version
  my @clone_coords = `cut -f 1,4,5,9 /wormsrv2/autoace/GFF_SPLITS/WS$gff_version/CHROMOSOME_*.clone_acc.gff`;
  foreach (@clone_coords){
    chomp;
    my ($chrom, $start, $end, $ninth) = split(/\t+/, $_);
    $chrom =~ s/[A-Z]+_//;
    $ninth =~ /S.+\s+\"(.+)\".+/;
    my $clone = $1;
    push(@{$clone_info{$clone}}, $chrom, $start, $end);
  }
  return %clone_info;
}

sub get_unique_from_array {
  my ($this, @array) = @_;
  my %seen=();
  my @new_array_no_dup;
  foreach (@array){push(@new_array_no_dup, $_) unless $seen{$_}++}
  return @new_array_no_dup;
}

sub get_overlapped_clones {
  my $this = shift;
  my %overlapped_clone;

  my $get_overlapped_clones=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/get_overlapped_clones.def" quit
EOF
  open (FH, "echo '$get_overlapped_clones' | $tace $curr_db |") || die "Couldn't access $curr_db\n";
  while (<FH>){
    chomp $_;
    my @clones = split(/\s+/, $_) if $_ =~ /^\".+/;
    for (my $i = 0; $i < scalar @clones; $i++){
      $clones[$i] =~ s/\"//g;
      push(@{$overlapped_clone{$clones[0]}}, $clones[$i]) if $i != 0;
    }
  }
  return %overlapped_clone;
}


1;
