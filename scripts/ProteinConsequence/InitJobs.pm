=head1 NAME

ProteinConsequence::InitJobs - set up VEP pipeline jobs

=cut

=head1 DESCRIPTION

Setup for WormBase Build variant mapping and VEP annotation pipeline.
Fetches all live variations for species and splits into batches.
Creates GFF files for VEP analysis from split GFFs.
Prepares VEP input files.
Dumps file containing feature positions/types.

=cut

package ProteinConsequence::InitJobs;

use strict;

use File::Path qw(make_path remove_tree);
use Data::Dumper;
use Wormbase;
use Storable 'store';

use base qw(ProteinConsequence::BaseProteinConsequence);


sub fetch_input {
    my $self = shift;
    
    my $output_dir = $self->required_param('output_dir');
    
    my $err;
    if (-d $output_dir) {
	remove_tree($output_dir, {error => \$err});
    }	
    make_path($output_dir, {error => \$err});
    die "Failed to create new output directory: " . Dumper($err) if $err && @$err;

    $self->dbc->disconnect_when_inactive(1);

    my $wb = Wormbase->new(
	-organism => $self->required_param('species'),
	-autoace  => $self->required_param('database'),
	-debug    => $self->required_param('debug'),
	-test     => $self->required_param('test'),
	);

    $self->construct_gff($wb);
    
    $self->prepare_input_files();

    my $variations = $self->get_all_allele_ids_table_maker($wb);

    $self->dbc->disconnect_when_inactive(0);


    my $batches = $self->split_vars_into_batches($variations);
    my @batch_ids = keys %$batches;
    
    $self->param('batch_info', [map { {batch_id => $_, variation_id_string => join(',', @{$batches->{$_}})} } @batch_ids]);
}


=head2 construct_gff

    Title:    construct_gff
    Function: creates GFF file for VEP input from split GFF files, dumps feature positions to file
    Args:     Wormbase object
    Returns:  n/a

=cut

sub construct_gff {
    my ($self, $wb) = @_;

    my $ace = Ace->connect(-path => $wb->autoace) or die "Could not create AcePerl connection";

    my @files;
    my $gff_dir = $self->required_param('gff_dir');
    my $out_file = $self->required_param('output_dir') . '/' . $self->required_param('species') . '.gff';
    my $fp_file = $self->required_param('output_dir') . '/' . $self->required_param('species') . '_feature_positions.dmp';
    my $pg_file = $self->required_param('output_dir') . '/' . $self->required_param('species') . '_pseudogenes.dmp';

    my $mRNA_details = $self->get_mRNA_details_and_genes($ace);
    my $noncoding_genes = $self->get_noncoding_genes($ace);

    for my $type ('gene',  'UTR', 'curated', 'miRNA', 'ncRNA', 'rRNA', 'scRNA', 'snoRNA', '7kncRNA',
		  'snRNA', 'snlRNA', 'tRNA', 'stRNA', 'piRNA', 'asRNA', 'lincRNA', 'Pseudogene', 'circRNA') {
	push @files, (glob("$gff_dir/*_${type}.gff"), "$gff_dir/${type}.gff");
    }

    my (%gene_positions, %cds_transcripts, %protein_ids);
    my $feature_positions = {};
    my $entries_printed = {};

    open (GFF_OUT, '>', $out_file);
    my ($intron, $exon);
    for my $file (@files) {
	open (GFF_IN, '<', $file);
	while (<GFF_IN>) {
	    next if $_ =~ /^#/;
	    chomp;
	    s/\"//g;
	    my @F = split;
	    if ($F[1] eq 'gene' and $F[2] eq 'gene') {
		$gene_positions{$F[9]}{'start'} = $F[3];
		$gene_positions{$F[9]}{'end'} = $F[4];
	    }
	    elsif ($F[1] eq 'Coding_transcript' and $F[2] eq 'coding_exon') {
		my $cds = $F[-1];
		my $transcript = $F[-4];
		unless (exists $protein_ids{$cds}) {
		    $protein_ids{$cds} = $self->get_protein_id($ace, $cds);
		}
		$cds_transcripts{$cds}{$transcript} = 1;
	    }
	    elsif ($F[2] eq 'Pseudogene') {
		$feature_positions->{$F[0]}{$F[3]}{$F[4]}{$F[9]} = 'pseudogene';
	    }

	    unless ($F[2] eq 'intron' or $F[1] eq 'gene' or ($F[1] eq 'Coding_transcript' and $F[2] eq 'coding_exon') or $F[2] eq 'region') {
		my $to_print;
		($to_print, $feature_positions, $entries_printed) = $self->convert_to_gff3(\%gene_positions, \%cds_transcripts, \%protein_ids,
											   $mRNA_details, $noncoding_genes,
											   $feature_positions, $entries_printed, @F);
		print GFF_OUT $to_print . "\n";
	    }
	}
	close (GFF_IN);
    }
    close (GFF_OUT);
	
    store $feature_positions, $fp_file;

    return;
}


=head2 convert_to_gff3

    Title:    convert_to_gff3
    Function: converts input from split GFF files into GFF3 format string to print
    Args:     hashref of gene positions ({$gene}{'start'} = $start, {$gene}{'end'} = $end),
              hashref mapping CDS to transcript IDs ({$cds}{$transcript} = 1),
              hashref mapping CDS to protein IDs ({$cds} = $protein),
              hashref of transcript details ({$transcript}{'start'} = $start,
                                             {$transcript}{'end'} = $end,
                                             {$transcript}{'gene'} = $gene),
              hashref mapping noncoding transcripts to genes ({$nc_transcript} = $nc_gene),
              hashref of feature positions ({$chr}{$start}{$end}{$id} = $type),
              hashref of GFF entries already printed ({$feature_id} = 1),
              array of columns from split GFF file
    Returns:  GFF3 format string to print, hashref of feature positions, hashref of GFF entries
              already printed

=cut

sub convert_to_gff3 {
    my ($self, $gene_positions, $cds_transcripts, $protein_ids, $mRNA_details, $noncoding_genes, $feature_positions, $entries_printed, @F) = @_;

    my $biotype = $F[2];
    if ($F[2] eq 'coding_exon') {
	$F[2] = 'CDS';
	$biotype = 'protein_coding';
    }
    elsif ($F[2] eq 'asRNA') {
	$F[2] = 'antisense_RNA';
    }
    elsif ($F[2] eq 'pre_miRNA') {
	$biotype = 'miRNA';
    }
    elsif ($F[2] eq 'Pseudogene') {
	if ($F[1] eq 'Pseudogene' or $F[1] eq 'tRNA_Pseudogene') {
	    $F[2] = 'pseudogenic_transcript';
	}
	elsif ($F[1] eq 'rRNA_Pseudogene') {
	    $F[2] = 'pseudogenic_rRNA';
	}
	$biotype = 'pseudogene';
    }
    elsif ($F[2] eq 'transposable_element_ncRNA') {
	$F[2] = 'ncRNA';
	$biotype = 'transposon';
    }
    elsif ($F[2] eq 'transposable_element_pseudogene') {
	$F[2] = 'pseudogenic_transcript';
	$biotype = 'transposon_pseudogene';
    }

    $F[1] = 'WormBase_VEP';

    if ($F[2] ne 'CDS' and $F[2] ne 'exon' and $F[2] ne 'miRNA' and $F[2] !~ /UTR/) {
	my $transcript = pop @F;
	pop @F;

	my $transcript_info = "ID=$transcript;Name=$transcript;biotype=$biotype";
	my $transcript_line = join("\t", @F, $transcript_info);
	my $gene_id = $noncoding_genes->{$transcript};
	if (defined $gene_id) {
	    my $gene_start = $gene_positions->{$gene_id}{'start'};
	    my $gene_end = $gene_positions->{$gene_id}{'end'};
	    $transcript_line = join("\t", $F[0], $F[1], 'gene', $gene_start, $gene_end, $F[5], $F[6], '.',
				    "ID=$gene_id;Name=$gene_id;gene_id=$gene_id;biotype=$biotype") . "\n$transcript_line" .
					";Parent=$gene_id";
	    $feature_positions->{$F[0]}{$gene_start}{$gene_end}{$gene_id} = 'gene';
	    $feature_positions->{$F[0]}{$F[3]}{$F[4]}{$transcript} = 'transcript'
		unless exists $feature_positions->{$F[0]}{$F[3]}{$F[4]}{$transcript};
	    
	}
	return ($transcript_line, $feature_positions, $entries_printed);	
    }

    if ($F[2] eq 'CDS') {
	my $cds = pop @F;
	my @gff_lines;
	for my $transcript (keys %{$cds_transcripts->{$cds}}) {
	    if (!exists $entries_printed->{$transcript}) {
		my $transcript_start = $mRNA_details->{$transcript}{'start'};
		my $transcript_end = $mRNA_details->{$transcript}{'end'};
		my $gene_id = $mRNA_details->{$transcript}{'gene'};
		unless (exists $entries_printed->{$gene_id}) {
		    my $gene_start = $gene_positions->{$gene_id}{'start'};
		    my $gene_end = $gene_positions->{$gene_id}{'end'};
		    push @gff_lines, join("\t", $F[0], $F[1], 'gene', $gene_start, $gene_end, $F[5], $F[6], '.',
					  "ID=$gene_id;gene_id=$gene_id;Name=$gene_id;biotype=protein_coding");
		    $feature_positions->{$F[0]}{$gene_start}{$gene_end}{$gene_id} = 'gene';
		    $entries_printed->{$gene_id} = 1;
		}
		push @gff_lines, join("\t", $F[0], $F[1], 'mRNA', $transcript_start, $transcript_end, $F[5],
				      $F[6], '.', "ID=$transcript;transcript_id=$transcript;Name=$transcript;Parent=$gene_id");
		$feature_positions->{$F[0]}{$transcript_start}{$transcript_end}{$transcript} = 'transcript';
		$entries_printed->{$transcript} = 1;
	    }
	    my $attr = "ID=$cds;Name=$cds;Parent=$transcript";
	    if (exists $protein_ids->{$cds} and $protein_ids->{$cds} ne '') { 
		$attr .= ';protein_id=' . $protein_ids->{$cds};
	    }
	    push @gff_lines, join("\t", $F[0], $F[1], $F[2], $F[3], $F[4], $F[5], $F[6], $F[7], $attr);
	}
	return (join("\n", @gff_lines), $feature_positions, $entries_printed);
    }

    if ($F[2] eq 'exon' and $F[-2] eq 'CDS') {
	my $cds = pop @F;
	pop @F;
	my @gff_lines;
	for my $transcript (keys %{$cds_transcripts->{$cds}}) {
	    push @gff_lines, join("\t", @F, "Parent=$transcript");
	}
	return (join("\n", @gff_lines), $feature_positions, $entries_printed);
    }

    my $transcript = pop @F;
    pop @F;
    if ($F[2] =~ /UTR/) {
	my $utr_line = join("\t", @F, "Parent=$transcript");
	$F[2] = 'exon';
	return ($utr_line . "\n" . join("\t", @F, "Parent=$transcript"), $feature_positions, $entries_printed);
    }

    return (join("\t", @F, "Parent=$transcript"), $feature_positions, $entries_printed);
	    
}


=head2 get_all_allele_ids_table_maker

    Title:    get_all_allele_ids_table_maker
    Function: creates TableMaker file and uses it to retrieve all live variations for species
    Args:     Wormbase object
    Returns:  arrayref of variation IDs

=cut

sub get_all_allele_ids_table_maker {
    my ($self, $wb) = @_;

    my $species = "Caenorhabditis " . lc($self->required_param('species'));
    my $query = "Sortcolumn 1\n\n" .
	"Colonne 1\n" .
	"Width 12\n" .
	"Mandatory\n" .
	"Visible\n" .
	"Class\n" .
	"Class Variation\n" .
	"From 1\n" .
	"Condition Flanking_sequences AND Live AND Species = \"$species\"";

    my $def_file = $self->required_param('output_dir') . '/all_allele_ids.def';
    open (DEF, '>', $def_file) or die "Could not open $def_file for writing";
    print DEF $query;
    close (DEF);
    
    my @list;
    
    my $tmfh = $wb->table_maker_query($self->required_param('database'), $def_file);
    while(<$tmfh>) {
        /^\"(\S+)\"$/ and do {
            print $1 . "\n";
            push @list, $1;
        }
    }
    close($tmfh);

    return \@list;
}


=head2 get_mRNA_details_and_genes

    Title:    get_mRNA_details_and_genes
    Function: gets transcript start/end positions and corresponding gene IDs
    Args:     Ace object
    Returns:  hashref of transcript details ({$transcript}{'start'} = $start,
                                             {$transcript}{'end'} = $end,
                                             {$transcript}{'gene'} = $gene)

=cut

sub get_mRNA_details_and_genes {
    my ($self, $ace) = @_;

    my $gff_dir = $self->required_param('gff_dir');
    my %mRNA_details;
    
    my @files = (glob("$gff_dir/*_UTR.gff"), "$gff_dir/UTR.gff");
    for my $file (@files) {
	open (UTR, '<', $file);
	while (<UTR>) {
	    next if $_ =~ /^#/;
	    chomp;
	    s/\"//g;
	    my @F = split;
	    my $transcript = $F[2] =~ /UTR/ ? $F[-1] : $F[-4];
	    $mRNA_details{$transcript}{'start'} = $F[3]
		unless exists $mRNA_details{$transcript}{'start'} and
		$mRNA_details{$transcript}{'start'} < $F[3];
	    $mRNA_details{$transcript}{'end'} = $F[4]
		unless exists $mRNA_details{$transcript}{'start'} and
		$mRNA_details{$transcript}{'end'} > $F[4];
	}
	close (UTR);
    }

    for my $mrna (keys %mRNA_details) {
	my $transcript = $ace->fetch(Transcript => $mrna);
	$mRNA_details{$mrna}{'gene'} = $transcript->Gene->name;
    }
    return \%mRNA_details;
}


=head2 get_noncoding_genes

    Title:    get_noncoding_genes
    Function: maps noncoding transcripts to gene IDs
    Args:     Ace object
    Returns:  hashref mapping noncoding transcripts to genes ({$transcript} = $gene)

=cut

sub get_noncoding_genes {
    my ($self, $ace) = @_;

    my $gff_dir = $self->required_param('gff_dir');
    my %noncoding_genes;
    my @files;
    for my $type ('pre_miRNA', 'ncRNA', 'rRNA', 'scRNA', 'snoRNA', 'snRNA', 'snlRNA', '7kncRNA',
		  'tRNA', 'stRNA', 'piRNA', 'asRNA', 'lincRNA', 'Pseudogene', 'circRNA') {
	push @files, (glob("$gff_dir/*_${type}.gff"), "$gff_dir/${type}.gff");
    }

    for my $file (@files) {
	open (GFF, '<', $file);
	while (<GFF>) {
	    next if $_ =~ /^#/;
	    chomp;
	    s/\"//g;
	    my @F = split;
	    my $object = $F[-2];
	    next if $object eq 'Sequence';
	    my $id = $F[-1];
	    next if exists $noncoding_genes{$id};
	    my $ace_obj = $object eq 'Transcript' ? $ace->fetch(Transcript => $id)
		: $ace->fetch(Pseudogene => $id);
	    if (!defined $ace_obj) {
		$self->write_to("Could not find Ace object for $object $id");
		next;
	    }
	    $noncoding_genes{$id} = $ace_obj->Gene->name if $ace_obj->Gene;
	}
	close (GFF);
    }

    return \%noncoding_genes;
}


=head2 get_protein_id

    Title:    get_protein_id
    Function: retrieves the protein ID corresponding to a given CDS
    Args:     Ace object, CDS ID string
    Returns:  protein ID string

=cut

sub get_protein_id {
    my ($self, $ace, $cds) = @_;

    my $agr_cds = $ace->fetch(CDS => $cds);
    if ($agr_cds->Corresponding_protein) {
	return $agr_cds->Corresponding_protein->Peptide->name;
    }

    return '';
}


=head2 prepare_input_files

    Title:    prepare_input_files
    Function: prepares GFF and FASTA files for input to VEP software
    Args:     none
    Returns:  n/a

=cut

sub prepare_input_files{
    my $self = shift;

    my $fasta = $self->required_param('fasta');
    my $species = $self->required_param('species');
    my $output_dir = $self->required_param('output_dir');
    my $gff = $output_dir . '/' . $species . '.gff';
    
    my @cmds = ("sort -k1,1 -k4,4n -k5,5n -t \$'\\t' $gff > $gff.sorted",
		"mv $gff.sorted $gff",
	        "bgzip -i $output_dir/$species.fa",
		"bgzip $gff",
		"tabix -p gff $gff.gz");
    if (-e $fasta) {
	unshift @cmds, "cp $fasta $output_dir/$species.fa";
    }
    elsif (-e $fasta . '.gz') {
	unshift @cmds, "gunzip -c $fasta.gz > $output_dir/$species.fa";
    }
    else {
	die "No compressed/uncompressed version of $fasta\n";
    }
    
    for my $cmd (@cmds) {
	my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
	die "$flat_cmd - exit code: $exit_code: $stderr" unless $exit_code == 0;
    }
    
    return;
}


=head2 split_vars_into_batches
    
    Title:    split_vars_into_batches
    Function: divides allele IDs into batches
    Args:     array ref of variation IDs
    Returns:  hashref mapping batch IDs to strings of allele IDs ({$batch_id} = $allele_ids_str)

=cut
    
sub split_vars_into_batches {
    my ($self, $variations) = @_;

    my %batches;
    my $batch_id_str = 'Batch0000';
    my $batch_id;
    my $batch_nr = 1;
    my $nr_batches = int(scalar @$variations / $self->required_param('batch_size')) + 1;
    for my $variation (@$variations) {
	$batch_nr++;
	$batch_id = substr($batch_id_str, 0, -length($batch_nr)) . $batch_nr;
	push @{$batches{$batch_id}}, $variation;
	$batch_nr = 0 if $batch_nr == $nr_batches;
    }

    return \%batches;
}


sub write_output {
    my $self = shift;
    
    $self->dataflow_output_id($self->param('batch_info'), 2);
}

1;
