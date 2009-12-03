#!/usr/bin/perl

use lib $ENV{CVS_DIR};
use Wormbase;
use strict;
use Ace;

# prop should get a -species option

my $wb = Wormbase->new(-autoace => glob('~wormpub/DATABASES/current_DB' ));
#my $db = Ace->connect( $wb->autoace ) or die "Can't open ace database:", Ace->error;
my $db = Ace->connect( glob('~wormpub/DATABASES/current_DB') ) or die "Can't open ace database:", Ace->error;

my $version = $db->version;


my (
    %NOTES,          %LOCUS, %GENBANK,      %CONFIRMED,
    %ORFEOME,        %GENES, %GENE_EXTENTS, %WORMPEP,
    %TRANSCRIPT2CDS, %genes_seen,         , %loci_seen
);

my $debug = 1 if $ENV{DEBUG};

# setting up things

print STDERR "getting confirmed genes\n" if $debug;
get_confirmed( $db, \%CONFIRMED );
print STDERR "getting genebank ids\n" if $debug;
get_genbank( $db, \%GENBANK );
print STDERR "getting transcript2cds\n" if $debug;
get_transcripts( $db, \%TRANSCRIPT2CDS );
print STDERR "getting wormpep\n" if $debug;
get_wormpep( $db, \%WORMPEP );
print STDERR "getting loci information\n" if $debug;
get_loci( $db, \%LOCUS );
print STDERR "getting genes\n" if $debug;
get_genes( $db, \%GENES );
print STDERR "getting notes\n" if $debug;
get_notes( $db, \%NOTES );
print STDERR "getting ORFEOME info\n" if $debug;
get_orfeome( $db, \%ORFEOME );

$db->close;

while (<>) {
    chomp;
    next if /^\#/;
    my ( $ref, $source, $method, $start, $stop, $score, $strand, $phase, $group ) = split /\t/;
    
    next if $source eq 'assembly_tag';    # don't want 'em, don't need 'em
    next if $method eq 'HOMOL_GAP';       # don't want that neither

    $ref   =~ s/^CHROMOSOME_//;
    $group =~ s/CHROMOSOME_//;

    $source = '' if $source eq '*UNKNOWN*';

    # Process top-level CDS and Transcript entries
    if (   $method =~ /Transcript|CDS|.*primary_transcript/
        && ( $source =~ /Coding_transcript|curated|.*RNA|miRNA/i )
        && $group =~ /[Transcript|CDS] "(\w+\.\d+[a-z]?\.?\d?)"/ ){
        ### Need to pick up transcript IDs, too

        my $match = $1;
        remember_gene_extents( $1, $ref, $start, $stop, $strand );
        my @notes;

        # Many of the hash lookups are keyed by CDS *only*
        my $lookup = ( $group =~ /Transcript/ ) ? $TRANSCRIPT2CDS{$match} : $match;

        # This probably needs to be more restrictive
        my $loci = $LOCUS{$lookup} || $LOCUS{$match};

        # WS133: Still required
        # ADD IN A SEPERATE ENTRY FOR EACH LOCUS
        if ($loci) {    # These are really genes...
            foreach (@$loci) {
                my $bestname = bestname($_);

                # Span for curated loci - only once per WBGene
                print join( "\t", $ref, 'curated', 'gene', $start, $stop, $score, $strand, $phase, "Locus $bestname" ), "\n"
                  if ( !$loci_seen{$bestname}++ );
            }
        }

        # WS133 - NOTES ARE PART OF CDS FEATURES, BUT NOT TRANSCRIPTS
        # Append some notes to top-level feature entries.
        unless ( $source eq 'curated' && $method eq 'CDS' ) {
            my $notes = $NOTES{$lookup} || $NOTES{$match};
            push @notes, map { qq(Note "$_") } @{$notes} if $notes;

            my $wormpep = $WORMPEP{$lookup} || $WORMPEP{$match};
#           push @notes, map { qq(WormPep "$_") } @{$wormpep} if $wormpep; # hmmm ....
            push @notes, qq(WormPep "$wormpep") if $wormpep;

            # This should be a translated WBGene ID
            my $locus_notes = $LOCUS{$lookup} || $LOCUS{$match};
            push @notes, map { qq(Note "${\&bestname($_)}") } @{$locus_notes} if $locus_notes;

            my $confirmed = $CONFIRMED{$lookup} || $CONFIRMED{$match};
            push @notes, qq(Prediction_status "$confirmed") if $confirmed;

            # Fudge factor: append WBGeneIDS to each transcript entry
            my $genes = $GENES{$lookup} || $GENES{$match};
            if ($genes) {
                push @notes, map {qq(Gene "$_")} @$genes;
            }
        }

        if ( $method =~ /.*primary_transcript/ ) {

            # Add a CDS attribute to Transcript entries
            # (Need to be able to associate CDSes with Transcripts)
            push @notes, qq(CDS "$lookup") if $TRANSCRIPT2CDS{$match};

            # Create a top level feature for Coding_transcript to accomodate
            # GBrowse modification on feature aggregation. Although Sanger
            # supplies Coding_transcript:protein_coding_primary_transcript
            # entries, we need individual entries of
            # Coding_transcript:Transcript to ensure proper aggregation.

            # What are the ramifications of this on things like RNAs?
            $group = join ' ; ', $group, @notes;
            if ( $method =~ /protein_coding_primary_transcript/ ) {
                print join( "\t",
                    $ref, 'Coding_transcript', 'Transcript', $start, $stop,
                    $score, $strand, $phase, $group ), "\n";
            }
            else {

                # This will be RNAs and such.
                # Should I also preserve the original entries? (This is basicallt the original entry)
                # Do I need a top level feature of method Transcript in order for aggregation to work correctly?
                
                print join( "\t", $ref, $source, $method, $start, $stop, $score, $strand, $phase, $group ), "\n";
            }
            next;
        }

        $group = join ' ; ', $group, @notes;
        print join( "\t", $ref, $source, $method, $start, $stop, $score, $strand, $phase, $group ), "\n";

        next;
    }

    # Skip Ant's fix for WBGenes
    next if $source eq 'Gene' && $method eq 'processed_transcript';

    if (   $method eq 'region'
        && $source eq 'Genomic_canonical'
        && $group =~ /Sequence "(\w+)"/ ){
        if ( my $accession = $GENBANK{$1} ) {
            $group .= qq( ; Note "Clone $1; Genbank $accession");
            print join( "\t",
                $ref, 'Genbank', $method, $start, $stop, $score, $strand,
                $phase, "Genbank \"$accession\"" ),"\n";
        }
    }

    next if ( $method eq 'intron' && $source =~ /^tRNAscan/ ) ; # messing up tRNA scanning

    if (   $method eq 'PCR_product'
        && $source eq 'Orfeome'
        && $group =~ /PCR_product "([^\"]+)"/ ) {
        my $amp = $ORFEOME{$1};
        $group .= qq( ; Amplified $amp) if defined $amp;
    }

    # fix variant fields: Variant "T" => Note "T"
    $group =~ s/(?:Variant|Insert) "(\w+)"/Note "$1"/;

    # fix reversed targets
    if ( $group =~ /Target \S+ (\d+) (\d+)/ ) {
        if ( $2 < $1 ) {
            $group =~ s/(\d+) (\d+)/$2 $1/;
            $strand = '-';
        }
    }
    print join( "\t",
        $ref,   $source, $method, $start, $stop,
        $score, $strand, $phase,  $group ),
        "\n";
}

# write out the extents of the genes
my %seenit;
for my $cds ( keys %GENES ) {
    ( my $base = $cds ) =~ s/[a-z]$//;
    next if $seenit{$base}++;
    next unless $GENE_EXTENTS{$base};
    my ( $seqid, $start, $stop, $strand ) = @{ $GENE_EXTENTS{$base} }{qw(seqid start stop strand)};
    for my $gene ( @{ $GENES{$cds} } ) {    # there *should* be only one gene per cds
        print join( "\t", $seqid, 'gene', 'processed_transcript', $start, $stop, '.', $strand, '.', qq(Gene "$gene") ), "\n";
    }
}

exit 0;

# remember the extreme left and right of transcripts
# for translating into gene extents.  This is called with
# CDS and UTR groups.
sub remember_gene_extents {
    my ( $group, $seqid, $start, $stop, $strand ) = @_;
    $group =~ s/[a-z]$//;    # get the base name of the CDS/UTR group
    
    $GENE_EXTENTS{$group}{start} = $start if !defined $GENE_EXTENTS{$group}{start} || $GENE_EXTENTS{$group}{start} > $start;
    $GENE_EXTENTS{$group}{stop}  = $stop if !defined $GENE_EXTENTS{$group}{stop} || $GENE_EXTENTS{$group}{stop} < $stop;
    $GENE_EXTENTS{$group}{strand} ||= $strand;
    $GENE_EXTENTS{$group}{seqid}  ||= $seqid;
}

# grab gene2cds.commondata
sub get_genes {
  my ($db,$hash) = @_;  # Keys are CDSes; values are WBGeneIDs
  #  my @genes = $db->fetch(-query=>'find Gene IS WBGene0* AND Corresponding_CDS');
  my @genes = $db->fetch(-query=>'find Gene');
  foreach my $obj (@genes) {
    my @cds = ($obj->Corresponding_CDS,$obj->Corresponding_Transcript);
    next unless @cds;
    foreach (@cds) {
      push @{$hash->{$_}},$obj;
    }
  }
}

sub get_transcripts {
    my ( $db, $hash ) = @_;    # Keys are transcript; values are CDSids
    my @transcripts = $db->fetch( -query => 'find Transcript' );
    foreach my $obj (@transcripts) {
        my $cds = $obj->Corresponding_CDS;
        $hash->{$obj} = $cds;
    }
}

# New gene model ( > WS123 approach)
# Fetch all genes that are also loci
sub get_loci {

    # hash keys are predicted gene names, values are one or more gene objects
    # (These gene objects will be translated into three-letter locus names prior to dumping)
    my ( $db, $hash ) = @_;

    # Fetch out all the cloned loci
    # This approach means that some genes will have duplicate entries (arising form the Other_names)
    my @loci = $db->fetch( -query => 'find Gene Molecular_info AND (CGC_name OR Other_name)' );

    foreach my $obj (@loci) {
        my @genomic = ( $obj->Corresponding_CDS, $obj->Corresponding_Transcript );
        @genomic >= 1 or next;
        foreach (@genomic) {
            push @{ $hash->{$_} }, $obj;
        }
    }
}

# allow lookup by wormpep id
sub get_wormpep {
    my ( $db, $hash ) = @_ ; # hash keys are predicted CDS names, values are one or more wormpep names
    $hash = $wb->FetchData( 'cds2wormpep', $hash );
}

# This could probably be shifted to the ?Gene class.
sub get_notes {
    my ( $db, $hash ) = @_ ; # hash keys are predicted gene names, values are one or more brief identifications
    my @genes = $db->fetch(
        -query   => 'find CDS Brief_identification',
        -filltag => 'Brief_identification'
    );

    # Should probably also look for notes attached to sequences, yes? As before...
    push(
        @genes,
        $db->fetch(
            -query   => 'find Sequence Brief_identification',
            -filltag => 'Brief_identification'
        )
    );
    
    push(
        @genes,
        $db->fetch(
            -query   => 'find Transcript Brief_identification',
            -filltag => 'Brief_identification'
        )
    );
    foreach my $obj (@genes) {
        my @notes = $obj->Brief_identification or next;

        next if ( $hash->{$obj} );
        $hash->{$obj} = \@notes;
    }
}

sub get_genbank {
    my ( $db, $hash ) = @_; # hash keys are cosmid names, values are genbank accessions (1 to 1)

    my @cosmids = $db->fetch(
        -query   => 'find Genome_Sequence Database',
        -filltag => 'Database'
    );

    for my $cosmid (@cosmids) {
        my @dbs = $cosmid->Database;
        foreach (@dbs) {
            foreach my $col ( $_->col ) {
                next unless $col eq 'NDB_AC';
                $hash->{$cosmid} = $col->right;
            }
        }
    }
}

# What is the confirmed status for CDSes?
# This has changed with WS133
sub get_confirmed {
    my ( $db, $hash ) = @_; # hash keys are predicted gene names, values are confirmation type
    my @confirmed = $db->fetch( -query => 'find CDS Prediction_status' );
    foreach my $obj (@confirmed) {
        my ($tag) = $obj->Prediction_status->row if ( $obj->Prediction_status );
        $hash->{$obj} = $tag;
    }
}

sub get_orfeome {
    my ( $db, $hash ) = @_;
    my @mv_primers = $db->fetch( -query => 'find PCR_Product mv*', -filltag => 'Amplified' );
    for my $obj (@mv_primers) {
        my $amplified = $obj->Amplified;
        $hash->{$obj} = $amplified;
    }
}

# Translate WB gene IDs into more useful molecular or three-letter names
sub bestname {
    my $gene = shift;
    # Public name is, oddly, never filled in.
    my $bestname =
         $gene->Public_name
      || $gene->CGC_name
      || $gene->Molecular_name
      || $gene->Sequence_name;
    return $bestname;
}
