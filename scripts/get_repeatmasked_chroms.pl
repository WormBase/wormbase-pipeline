#!/usr/local/ensembl/bin/perl -w
#
# agp2ensembl
#
# Cared for by Simon Potter
# (C) GRL/EBI 2001
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#
# modified for reading in .agp files for worm ensembl

=pod

=head1 NAME

agp2ensembl.pl

=head1 SYNOPSIS

get_repeatmasked_chroms.pl -agp ~wormpipe/Elegans/WSXXX.agp

=head1 DESCRIPTION

extracts RepeatMasked sequence from DNA database using specified AGP file

=head1 OPTIONS

    -agp     agp file

=head1 CONTACT

ar2@sanger.ac.uk

=cut

use lib '/nfs/farm/Worms/Ensembl/ensembl-pipeline/modules';
use lib '/nfs/farm/Worms/Ensembl/ensembl/modules';
use lib '/nfs/disk100/humpub/modules/PerlModules';
use lib (-e '/wormsrv2/scripts') ? '/wormsrv2/scripts' : $ENV{'CVS_DIR'};


use strict;
use Getopt::Long;
#use Bio::Root::RootI;
#use Bio::Seq;
#use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Clone;
use Bio::EnsEMBL::RawContig;
#use Bio::SeqIO;
#use Hum::NetFetch qw( wwwfetch );
#use Hum::EMBL;
use Log_files;
use Wormbase;


my $agp;
my $test;
my $out_dir;
my $debug;

&GetOptions(
	    "agp:s"     => \$agp,
	    "output:s"  => \$out_dir,
	    "test"      => \$test,
	    "debug:s"   => \$debug
	   );


my $log;
if ( $test ) {
  $log = Log_files->make_log("$out_dir/CHROMOSOMES/get_repeatmasked_chroms.log");
}
else {
  $log = Log_files->make_build_log( $debug );
}


unless ($agp) {
    print STDERR "Must specify agp file\n";
    exit 1;
}

$out_dir = "/wormsrv2/autoace/CHROMOSOMES" unless $out_dir;
die "cant write to $out_dir\t$!\n" unless (-w $out_dir );
die "cant read agp file $agp\t$!\n" unless (-r $agp );

# open connection to EnsEMBL DB
my $dbobj;

$log->write_to("Connecting to worm_dna\n");

$dbobj = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
						       '-host'   => 'ecs1f',
						       '-user'   => 'wormro',
						       '-dbname' => 'worm_dna'
						      )
  or die "Can't connect to Database worm_dna";


my $clone_adaptor  = $dbobj->get_CloneAdaptor();
my $contig_adaptor = $dbobj->get_RawContigAdaptor();

$log->write_to("Building chromosomes\n");

open (AGP,"<$agp") or die "cant open $agp\n";
my %chrom_seq;

while ( <AGP> ) {

  #  X       1       2278    1       F       AL031272.3      1       2278    +
  my @data = split;
  my $agp_chrom = $data[0];
  my ($acc) = $data[5] =~ /(\w+)\./;
  my $length = $data[7];

  my $clone = $clone_adaptor->fetch_by_accession($acc);
  die "no clone in database for $acc\n" unless $clone;
  my $maskedseq = lc $clone->get_all_Contigs->[0]->get_repeatmasked_seq->seq; # only ever one contig/clone.
  $chrom_seq{$agp_chrom} .= substr($maskedseq,0,$length);
}
close AGP;


foreach ( sort keys %chrom_seq ) {
  my $text = "\t$_\t".length $chrom_seq{$_};
  $log->write_to("$text");
  $log->write_to("\n");
}

print STDERR "Outputting     ";

foreach my $seq (sort keys %chrom_seq ) {

  my $outfile = "$out_dir"."/CHROMOSOME_${seq}_masked.dna";
  open (OUT,">$outfile") or die "cant write $outfile\t$!\n";

  $log->write_to("\twriting chromosome $seq\n");
  print OUT ">CHROMOSOME_${seq} 1 ",length($chrom_seq{$seq}),"\n";
  my $width = 50;
  my $start_point = 0;
  while ( $start_point + $width < length( $chrom_seq{$seq} ) ) {
    print OUT substr($chrom_seq{$seq}, $start_point, $width ),"\n";
    $start_point += $width;
  }

  print OUT substr($chrom_seq{$seq}, $start_point),"\n";
  close OUT;
  system("gzip $outfile");
}

$log->write_to("Done\n");
$log->mail;

exit(0);
