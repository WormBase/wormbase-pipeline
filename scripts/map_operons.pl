#!/usr/local/bin/perl5.8.0 -w
#
# map_operons.pl

# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2006-02-23 16:02:01 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Data::Dumper;
use Getopt::Long;
use File::Copy;
use Ace;
use Log_files;

my ($debug, $store, $test, $noload);

GetOptions(
	   'debug=s'  => \$debug,
	   'store=s'  => \$store,
	   'test'     => \$test,
	   'noload'   => \$noload,
);

############################
# recreate configuration   #
############################
my $wb;
if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug, -test => $test, ) }

my $log = Log_files->make_build_log($wb);
my $acefile = $wb->acefiles."/operon_coords.ace";

my @chromosomes = $test ? qw(III) : qw(I II III IV V X);
my %gene_span;
foreach (@chromosomes){
  open (GS,"<".$wb->gff_splits."/CHROMOSOME_${_}_gene.gff") or $log->log_and_die("Cant open ".$wb->gff_splits."/CHROMOSOME_${_}_gene.gff :$!\n");
  while (<GS>) {
    # CHROMOSOME_III  gene    gene    16180   17279   .       +       .       Gene "WBGene00019182"
    my @data = split;
    next unless ($data[2] eq 'gene');
    $data[9] =~ s/\"//g;
    my $gene = $data[9];
    $gene_span{$gene}->{'chrom'} = $data[0];
    $gene_span{$gene}->{'start'} = $data[3];
    $gene_span{$gene}->{'end'}   = $data[4];
    $gene_span{$gene}->{'strand'}= $data[6];
  }
}

open (OUT,">$acefile") or $log->log_and_die("cant open $acefile : $!\n");
my $db = Ace->connect(-path => $wb->autoace) or $log->log_and_die("cant connect to ".$wb->autoace." :".Ace->error."\n");
my @operons = $db->fetch('Operon' => '*');
foreach my $operon(@operons) {
  my @genes = map($_->name, $operon->Contains_gene);
  my ($op_start, $op_end, $op_strand);
  foreach my $gene (@genes) {
    $op_start =  $gene_span{$gene}->{'start'} if (!(defined $op_start) or $op_start > $gene_span{$gene}->{'start'});
    $op_end   =  $gene_span{$gene}->{'end'}   if (!(defined $op_end)   or $op_end   < $gene_span{$gene}->{'end'});
    $op_strand = $gene_span{$gene}->{'strand'}if(!(defined $op_strand) or $op_strand eq $gene_span{$gene}->{'strand'});
    $op_strand = $gene_span{$gene}->{'chrom'};
  }

  print OUT "\nSequence : $op_strand\nOperon $operon ";
  ($op_start < $op_end) ? print OUT "$op_start $op_end" : print OUT "$op_end $op_start";
  print OUT "\n";
}

$wb->load_to_database($wb->autoace,"$acefile",'operon_span') unless $noload;

$log->mail;
exit(0);
