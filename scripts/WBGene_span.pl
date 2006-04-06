#!/usr/local/bin/perl5.8.0 -w
#
# WBGene_span.pl
#
# by Anthony Rogers
#
# Creates SMapped Gene spans for Gene objects
#
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2006-04-06 12:33:14 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Coords_converter;
use Getopt::Long;
use Log_files;
use Storable;

my ( $database, $test, $store, $gff, $no_ace, $debug, $gff_file, $chromosome, $prepare, $no_gff );

GetOptions(
    'database=s'   => \$database,
    'test'         => \$test,
    "store:s"      => \$store,
    'gff'          => \$gff,
    'no_ace'       => \$no_ace,
    'debug=s'      => \$debug,
    'gff_file=s'   => \$gff_file,
    'chromosome=s' => \$chromosome,
    'prepare'      => \$prepare,
    'no_gff'       => \$no_gff
);

my $wormbase;
if ($store) { $wormbase = retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wormbase = Wormbase->new( -debug => $debug, -test => $test, ) }

my $log = Log_files->make_build_log($wormbase);

$database = $wormbase->autoace unless $database;
$log->write_to("Generating WBGene spans from database $database\n");

# get the required GFF files in place
if ($prepare) {

    # gff dump Coding_transcript
    print "dumping CDSes\n" if $test;

    $wormbase->run_command( "perl $ENV{'CVS_DIR'}/GFF_method_dump.pl -database $database -method Coding_transcript -dump_dir " . $wormbase->gff_splits );

    # cat them to the main gff files
    # foreach my $chrom ( @chromosomes ) {
    # &run_command("cat $database/CHROMOSOMES/CHROMOSOME_${chrom}.Coding_transcript.gff >> $database/CHROMOSOMES/CHROMOSOME_${chrom}.gff");
}
else {
    my %CDS_data;
    my $cds;
    my %gene2seq;

    my $def_file = "$database/wquery/WBGene2seq.def";
    $log->log_and_die("Tablemaker definition file WBGene2seq.def is required. Not present in $database/wquery\n") unless ( -e "$def_file" );

    my $tace   = $wormbase->tace;
    my $coords = Coords_converter->invoke( $database, undef, $wormbase );

    my $command = "Table-maker -p $def_file\nquit\n";

    open( FH, "echo '$command'| $tace $database | " );

    while (<FH>) {
        s/{\t}*/\t/g;
        s/\"//g;
        my @data = split;
        next unless $data[1];
        my $gene = shift @data;
        push( @{ $gene2seq{$gene} }, @data );
    }

    close FH;

    my @chromosomes = qw( I II III IV V X MtDNA);
    @chromosomes = ("$chromosome") if $chromosome;
    foreach my $chrom (@chromosomes) {
        my %gene_coords;

        foreach my $method (qw(Coding_transcript Non_coding_transcript tRNAscan-SE-1.23 miRNA)) {
            print "checking $method \n" if $test;
            $gff_file = $wormbase->gff_splits . "/CHROMOSOME_${chrom}_$method.gff" unless $gff_file;
            open( GFF, "<$gff_file" ) or die "cant open GFF $gff_file\n";
            while (<GFF>) {
                my @data = split;
                if (   ( $data[1] eq "Coding_transcript" )
                    or ( $data[1] eq "Non_coding_transcript" )
                    or ( $data[1] eq "Pseudogene" )
                    or ( $data[1] eq "curated" )
                    or ( $data[1] eq "tRNAscan-SE-1.23" )
                    or ( $data[1] eq "tRNA" )
                    or ( $data[1] eq "snRNA" )
                    or ( $data[1] eq "miRNA" )
                    or ( $data[1] eq "rRNA" )
                    or ( $data[1] eq "scRNA" )
                    or ( $data[1] eq "snoRNA" )
                    or ( $data[1] eq "tRNA" )
                    or ( $data[1] eq "stRNA" )
                    or ( $data[1] eq "snRNA" ) )
                {
                    next if ( $data[2] eq "exon" or $data[2] eq "coding_exon" or $data[2] eq "intron" );
                    my ($gene) = $data[9] =~ /\"(\w+\.\w+)/;
                    push( @{ $gene_coords{$gene} }, [ ( $data[3], $data[4], $data[6] ) ] );
                }
            }
            close GFF;
            undef $gff_file;
        }
        my %gene_span;
        foreach my $WBgene ( keys %gene2seq ) {
            foreach my $CDS ( @{ $gene2seq{$WBgene} } ) {
                next unless $gene_coords{$CDS};
                foreach my $transcript ( @{ $gene_coords{$CDS} } ) {
                    if ( !( defined $gene_span{$WBgene} ) ) {
                        $gene_span{$WBgene}->{'min'}    = $transcript->[0];
                        $gene_span{$WBgene}->{'max'}    = $transcript->[1];
                        $gene_span{$WBgene}->{'strand'} = $transcript->[2];
                    }
                    else {
                        if ( $transcript->[0] < $gene_span{$WBgene}->{'min'} ) {
                            $gene_span{$WBgene}->{'min'} = $transcript->[0];
                        }
                        if ( $gene_span{$WBgene}->{'max'} < $transcript->[1] ) {
                            $gene_span{$WBgene}->{'max'} = $transcript->[1];
                        }
                    }
                }
            }
        }
        if ($gff) {
            open( OUTGFF, $wormbase->gff_splits . "/CHROMOSOME_${chrom}_WBgene.gff" ) or do { $log->write_to("cant open output\n"); die "cant open output\n"; }
        }
        my $acefile = $wormbase->acefiles . "/WBgene_spans_${chrom}.ace";
        unless ($no_ace) {
            open( ACE, ">$acefile" ) or do { $log->write_to("cant open output $acefile:\t$!\n"); die "cant open output $acefile:\t$!\n"; }
        }
        foreach my $gene ( keys %gene_span ) {
            if ($gff) {
                print OUTGFF "CHROMOSOME_${chrom}\tgene\tgene\t";
                print OUTGFF $gene_span{$gene}->{'min'}, "\t", $gene_span{$gene}->{'max'}, "\t";
                print OUTGFF ".\t", $gene_span{$gene}->{'strand'}, "\t.\tGene \"$gene\"\n";
            }

            # write S-map details back to database
            unless ($no_ace) {
                my @coords;
                @coords =
                    $gene_span{$gene}->{'strand'} eq "+"
                  ? $coords->LocateSpan( $chrom, $gene_span{$gene}->{'min'}, $gene_span{$gene}->{'max'} )
                  : $coords->LocateSpan( $chrom, $gene_span{$gene}->{'max'}, $gene_span{$gene}->{'min'} );
                print ACE "\nSequence : \"$coords[0]\"\n";
                print ACE "Gene_child $gene $coords[1] $coords[2]\n";
            }
        }
        close ACE;
        close OUTGFF;
        $wormbase->load_to_database( "$database", "$acefile", "WBGene_span" ) unless ($no_ace);    # I know this would be better done all together but . .
    }
}

$wormbase->run_script( "dump_gff_batch.pl -database " . $wormbase->autoace . " -methods gene -dump_dir " . $wormbase->gff_splits, $log ) unless $no_gff;
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

