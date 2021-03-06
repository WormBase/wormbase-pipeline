#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  mlagan_test.pl
#
#        USAGE:  ./mlagan_test.pl
#
#  DESCRIPTION:  test the MLAGAN database
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:   (Michael Han), <mh6@sanger.ac.uk>
#      COMPANY:
#      VERSION:  1.0
#      CREATED:  10/08/06 10:32:12 BST
#     REVISION:  1
#===============================================================================

use strict;


use Bio::EnsEMBL::Registry;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::LocatableSeq;
Bio::EnsEMBL::Registry->load_all();

# adaptors setup
my $mlss_a = Bio::EnsEMBL::Registry->get_adaptor( 'Compara', 'compara', 'MethodLinkSpeciesSet' );
my $gdb_a  = Bio::EnsEMBL::Registry->get_adaptor( 'Compara', 'compara', 'GenomeDB' );
my $ga_a   = Bio::EnsEMBL::Registry->get_adaptor( 'Compara', 'compara', 'GenomicAlignBlock' );
my $dnaf_a = Bio::EnsEMBL::Registry->get_adaptor( 'Compara', 'compara', 'DnaFrag' );

# genome databases setup
my $edb = $gdb_a->fetch_by_name_assembly('Caenorhabditis elegans');
my $bdb = $gdb_a->fetch_by_name_assembly('Caenorhabditis briggsae');
my $rdb = $gdb_a->fetch_by_name_assembly('Caenorhabditis remanei');
my $ndb = $gdb_a->fetch_by_name_assembly('Caenorhabditis brenneri');
my $pdb = $gdb_a->fetch_by_name_assembly('Pristionchus pacificus');
my $jdb = $gdb_a->fetch_by_name_assembly('Caenorhabditis japonica');
my $brdb = $gdb_a->fetch_by_name_assembly('Brugia malayi');

my $mlss = $mlss_a->fetch_by_method_link_type_GenomeDBs( 'PECAN', [ $edb, $bdb, $rdb,$ndb,$brdb,$pdb,$jdb] );

my $alignIO = Bio::AlignIO->newFh(
    -interleaved => 0,
    -fh          => \*STDOUT,
    -format      => 'clustalw',
    -idlength    => 40
);

# get a DNA fragment
my $elegans_I = $dnaf_a->fetch_by_GenomeDB_and_name( $edb, $ARGV[0] );

## throws a warning because of my meta table setup

my $all_blocks = $ga_a->fetch_all_by_MethodLinkSpeciesSet_DnaFrag( $mlss, $elegans_I, $ARGV[1], $ARGV[2] );
my @all_aligns;

# for all aligned blocks
foreach my $align (@$all_blocks) {

    my $simple_align = Bio::SimpleAlign->new();
    $simple_align->id( "ID#" . $align->dbID );

    # get all alignments
    foreach my $this ( @{ $align->get_all_GenomicAligns() } ) {
        #		$this->dnafrag->name, '[', $this->dnafrag_start,'-', $this->dnafrag_end, '](', $this->dnafrag_strand, ")\n",
        my $seq_name = $this->dnafrag->genome_db->name;
        $seq_name =~ s/(.)\w* (...)\w*/$1$2/;
        $seq_name .= '-'.$this->dnafrag->name;
	my $strand=$this->dnafrag_strand >0 ?'+':'-';
	$seq_name .= "($strand)";
        my $aligned_sequence = $this->aligned_sequence;
        my $seq              = Bio::LocatableSeq->new(
            -SEQ    => $aligned_sequence,
            -START  => $this->dnafrag_start,
            -END    => $this->dnafrag_end,
            -ID     => $seq_name,
            -STRAND => $this->dnafrag_strand
        );
        $simple_align->add_seq($seq);
    }
    push( @all_aligns, $simple_align );
}

foreach my $aln (@all_aligns) {
    print $alignIO $aln;
}
