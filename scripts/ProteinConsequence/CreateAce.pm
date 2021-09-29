=head1 NAME

ProteinConsequence::CreateAce - creates AceDB file

=cut

=head1 DESCRIPTION

Creates .ace files with features associated with variations and molecular consequences.
Compares variation/gene connections already stored in AceDB with those generated by pipeline.

=cut

package ProteinConsequence::CreateAce;

use strict;

use Wormbase;
use Ace;
use Storable 'retrieve';

use base ('ProteinConsequence::BaseProteinConsequence');

sub run {
    my $self = shift;

    my $batch_id = $self->required_param('batch_id');
    my $working_dir = $self->required_param('output_dir') . "/$batch_id/";
    my $ace_file = $working_dir . $batch_id . '.ace';
    my $vep_output = $working_dir . $batch_id . '.vep.txt';
    my $mapped_alleles_dump = $working_dir . $batch_id . '_mapped_alleles.dmp';
    my $mapped_alleles = retrieve($mapped_alleles_dump);

    my $wb = Wormbase->new(
	-autoace  => $self->required_param('database'),
	-organism => $self->required_param('species'),
	-debug    => $self->required_param('debug'),
	-test     => $self->required_param('test'),
	);

    my $ace = Ace->connect(-path => $wb->autoace) or die "Could not create AcePerl connection";

    $self->write_to("Writing basic position information for $batch_id\n") if $self->required_param('debug');
    my $fh = new IO::File ">$ace_file" || die($!);
    while (my($var, $allele) = each %$mapped_alleles) {
	print $fh "Sequence : \"$allele->{clone}\"\nAllele $var $allele->{clone_start} $allele->{clone_stop}\n\n";
    }

    my ($transcripts, $hgvsg) = $self->get_transcript_consequences($vep_output);

    my ($gene_alleles, $allele_genes, $pseudogenes);
    ($gene_alleles, $allele_genes, $transcripts, $pseudogenes) = 
	$self->get_overlapping_features($mapped_alleles, $transcripts);


    # compare old<->new genes
    $self->write_to("Comparing old mappings to new for $batch_id\n") if $self->required_param('debug');
    $self->compare_gene_connections($mapped_alleles, $allele_genes);
    
    $self->write_to("Writing gene details for $batch_id\n") if $self->required_param('debug');
    for my $var (keys %$allele_genes) {
	print $fh "Variation : \"$var\"\n";
	for my $gene (@{$allele_genes->{$var}}) {
	    print $fh "Gene $gene\n";
	}
	print $fh "\n";

	for my $gene (@{$allele_genes->{$var}}) {
	    print $fh "Gene : $gene\nAllele $var\n\n";
	}
    }


    my %pseudogenes_printed;
    $self->write_to("Writing transcript details for $batch_id\n") if $self->required_param('debug');
    for my $var (keys %$transcripts) {
	next unless $self->add_transcript_consequences_for_variation($var, $mapped_alleles->{$var}{allele});
	print $fh "Variation : \"$var\"\n";
	print $fh 'HGVSg "' . $hgvsg->{$var} . "\"\n" if exists $hgvsg->{$var};
	for my $transcript (keys %{$transcripts->{$var}}) {
	    my $type = 'Transcript';
	    if (exists $pseudogenes->{$var}{$transcript}) {
		$type = 'Pseudogene';
		$pseudogenes_printed{$var}{$transcript}++;
	    }
	    for my $tag (keys %{$transcripts->{$var}{$transcript}}) {
		print $fh "$type $transcript $tag\n";
	    }
	}
	print $fh "\n";
    }

    
    $self->write_to("Writing pseudogene details for $batch_id\n") if $self->required_param('debug');
    for my $var (keys %$pseudogenes) {
	my $var_header = "Variation : \"$var\"\n";
	my $var_header_printed = 0;
	for my $pseudogene (keys %{$pseudogenes->{$var}}) {
	    next if $pseudogenes_printed{$var}{$pseudogene};
	    print $fh $var_header unless $var_header_printed;
	    $var_header_printed = 1;
	    print $fh "Pseudogene $pseudogene\n";
	}
	print $fh "\n" if $var_header_printed;
    }

    close($fh) or $self->log_and_die("Could not close $ace_file after writing\n");
    $ace->close;

}


=head2 add_transcript_consequences_for_variation

    Title:    add_transcript_consequences_for_variation
    Function: determines whether transcript consequences should be added to output
    Args:     variation ID, AceDB allele object
    Returns:  boolean indicating whether or not to add consequences

=cut

sub add_transcript_consequences_for_variation {
    my ($self, $var, $allele) = @_;
    my ($substitution, $deletion, $insertion, %allele_type);
    if (!defined $allele) {
	$self->write_to(sprintf("WARNING: %s - No allele defined",
				$var));
	return 0;
    }

    foreach my $tp ($allele->Type_of_mutation) {
	my @list;
	if ($tp->right) {
	    push @list, lc($tp->right);
	    if ($tp->right->right) {
		push @list, lc($tp->right->right);
	    }
	}
	if ($tp eq 'Substitution') {
	    $allele_type{Substitution} = 1;
	    $substitution = \@list;
	} elsif ($tp eq 'Insertion') {
	    $allele_type{Insertion} = 1;
	    $insertion = $list[0];
	} elsif ($tp eq 'Deletion') {
	    $allele_type{Deletion} = 1;
	    $deletion = 1;
	}
    }

    if (defined $substitution and (defined $insertion or defined $deletion)) {
      $self->write_to(sprintf("WARNING: %s - Substitution defined together with Insertion/Deletion, ".
                             "so will not calculate prot effect (%s ; Remark:%s)\n", 
                             $var, $allele->Public_name, $allele->Remark));
      return 0;

    }

    if (defined $substitution and scalar(@$substitution) != 2) {
      $self->write_to(sprintf("WARNING: %s - substitution has missing/empty FROM and/or TO, ".
                             "so will not calculate prot effect (%s; Remark:%s)\n", 
                             $var, $allele->Public_name, $allele->Remark));
      return 0;
    }

    if (defined $substitution and ($substitution->[0] =~ /\d|\s/ or $substitution->[1] =~ /\d|\s/)) {
      $self->write_to(sprintf("WARNING: %s - small substitution has numbers/spaces in FROM/TO, ".
                             "so will not calculate prot effect (%s; Remark:%s)\n", 
                             $var, $allele->Public_name, $allele->Remark));
      return 0;
    } 

    return 1;
}


=head2 compare_gene_connections

    Title:    compare_gene_connection
    Function: compares variation/gene connections existing in AceDB with those generated by pipeline
    Args:     hashref of mapped alleles, hashref of allele to gene mappings
    Returns:  n/a

=cut

sub compare_gene_connections {
    my ($self, $old, $new) = @_;
    
    my %check;
    $self->write_to(sprintf('Getting old gene connections for %s', 
			 $self->required_param('batch_id'))) if $self->required_param('debug');
    for my $allele (keys %$old) {
	for my $gene ($old->{$allele}->{allele}->Gene) {
	    next if (qq/${\$old->{$allele}->{allele}->at("Affects.Gene.$gene")->col(1)}/ eq 'Genomic_neighbourhood');
	    next if (qq/${\$old->{$allele}->{allele}->at("Affects.Gene.$gene")->col(1)}/ eq 'Regulatory_feature');
	    
	    $check{$allele}{$gene} = 1;
	}
    }
    
    $self->write_to(sprintf('Getting new gene connections for %s', 
			 $self->required_param('batch_id'))) if $self->required_param('debug');
    for my $allele (keys %$new) {
	for my $gene(@{$new->{$allele}}) {
	    $check{$allele}{$gene} += 2;
	}
    }

    
    $self->write_to(sprintf('Comparing old and new gene connections for %s', 
			 $self->required_param('batch_id'))) if $self->required_param('debug');
    for my $allele (keys %check) {
	for my $gene (keys %{$check{$allele}}) {
	    if ($check{$allele}{$gene} == 1) {
		my $remark = $old->{$allele}->{allele}->Remark;
		$self->write_to("ERROR: $allele (${\$old->{$allele}->{allele}->Public_name}) -> $gene connection is only in geneace (Remark: $remark)\n");
	    }
	    elsif ($check{$allele}{$gene} == 2) {
		$self->write_to("INFO: $allele -> $gene connection created by script") if
		    $self->required_param('debug');
	    }
	    else {
		die "Gene comparison failed" unless $check{$allele}{$gene} == 3;
	    }
	}
    }

    return;
}


=head2 get_overlapping_features

    Title:    get_overlapping_features
    Function: retrieve features overlapping variation
    Args:     hashref of mapped alleles, hashref of transcript annotations
    Returns:  hashref mapping genes to variations, hashref mapping variations to genes, 
              hashref of transcript annotations, hashref mapping variations to pseudogenes

=cut

sub get_overlapping_features {
    my ($self, $mapped_alleles, $transcripts) = @_;

    my (%gene_alleles, %allele_genes, %pseudogenes);
    my $positions_file = $self->required_param('output_dir') . '/' . $self->required_param('species') . '_feature_positions.dmp';
    my $positions = retrieve($positions_file);

    while (my($var, $allele) = each %$mapped_alleles) {
	my $chr = $allele->{chromosome};
	for my $fstart (sort {$a<=>$b} keys %{$positions->{$chr}}) {
	    last if $allele->{stop} < $fstart;
	    for my $fend (sort {$a<=>$b} keys %{$positions->{$chr}{$fstart}}) {
		next if $allele->{start} > $fend;
		for my $fid (keys %{$positions->{$chr}{$fstart}{$fend}}) {
		    my $ftype = $positions->{$chr}{$fstart}{$fend}{$fid};
		    if ($ftype eq 'gene') {
			$gene_alleles{$fid}{$var} = 1;
			$allele_genes{$var}{$fid} = 1;
		    }
		    elsif ($ftype eq 'transcript') {
			$transcripts->{$var}{$fid}{''} = 1 unless exists $transcripts->{$var}{$fid};
		    }
		    elsif ($ftype eq 'pseudogene') {
			$pseudogenes{$var}{$fid} = 1;
		    }
		}
	    }
	}
    }

    return (restructure_hash(\%gene_alleles), restructure_hash(\%allele_genes), $transcripts, \%pseudogenes);
}


=head2 get_transcript_consequences

    Title:    get_transcript_consequences
    Function: parses VEP output to retrieve transcript annotations
    Args:     VEP output filename (string)
    Returns:  hashref of transcript annotations, hashref of HGVSg IDs

=cut

sub get_transcript_consequences {
    my ($self, $vep_output) = @_;

    my %severity_ranking = ('intergenic_variant'                 => 1,
			    'feature_truncation'                 => 2,
			    'regulatory_region_variant'          => 3,
			    'feature_elongation'                 => 4,
			    'regulatory_region_amplification'    => 5,
			    'regulatory_region_ablation'         => 6,
			    'TF_binding_site_variant'            => 7,
			    'TFBS_amplification'                 => 8,
			    'TFBS_ablation'                      => 9,
			    'downstream_gene_variant'            => 10,
			    'upstream_gene_variant'              => 11,
			    'non_coding_transcript_variant'      => 12,
			    'NMD_transcript_variant'             => 13,
			    'intron_variant'                     => 14,
			    'non_coding_transcript_exon_variant' => 15,
			    '3_prime_UTR_variant'                => 16,
			    '5_prime_UTR_variant'                => 17,
			    'mature_miRNA_variant'               => 18,
			    'coding_sequence_variant'            => 19,
			    'synonymous_variant'                 => 20,
			    'stop_retained_variant'              => 21,
			    'start_retained_variant'             => 22,
			    'incomplete_terminal_codon_variant'  => 23,
			    'splice_region_variant'              => 24,
			    'protein_altering_variant'           => 25,
			    'missense_variant'                   => 26,
			    'inframe_deletion'                   => 27,
			    'inframe_insertion'                  => 28,
			    'transcript_amplification'           => 29,
			    'start_lost'                         => 30,
			    'stop_lost'                          => 31,
			    'frameshift_variant'                 => 32,
			    'stop_gained'                        => 33,
			    'splice_donor_variant'               => 34,
			    'splice_acceptor_variant'            => 35,
			    'transcript_ablation'                => 36,
	);
    

    my (%transcripts, %hgvsg, %worst_rankings);
    open (VEP, '<', $vep_output);
    while (<VEP>) {
	next if $_ =~ /^#/;
	chomp;
	my ($var, $pos, $allele, $gene, $feature, $feature_type, $consequence, $cdna_pos, $cds_pos,
	    $prot_pos, $aas, $codons, $existing_var, $attr) = split("\t", $_);

	next unless $feature_type eq 'Transcript';

	my @consequences = split(',', $consequence);
	my $worst_consequence_rank = 0;
	for my $single_consequence (@consequences) {
	    $worst_consequence_rank = $severity_ranking{$single_consequence} if
		$severity_ranking{$single_consequence} > $worst_consequence_rank;
	}
	next if exists $worst_rankings{$var}{$feature} and $worst_rankings{$var}{$feature} >= $worst_consequence_rank;
	$worst_rankings{$var}{$feature} = $worst_consequence_rank;
	
	my %attributes = split /[;=]/, $attr;
	$transcripts{$var}{$feature}{"VEP_consequence \"$consequence\""}++;
	
	$transcripts{$var}{$feature}{'VEP_impact "' . $attributes{'IMPACT'} . '"'} = 1 if exists $attributes{'IMPACT'};
	$transcripts{$var}{$feature}{'HGVSc "' . $attributes{'HGVSc'} . '"'} = 1 if exists $attributes{'HGVSc'};
	if (exists $attributes{'HGVSp'}) {
	    my $hgvsp = $attributes{'HGVSp'};
	    $hgvsp =~ s/%3D/=/;
	    $transcripts{$var}{$feature}{'HGVSp "' . $hgvsp . '"'} = 1;
	}
	$transcripts{$var}{$feature}{'Intron_number "' . $attributes{'INTRON'} . '"'} = 1 if exists $attributes{'INTRON'};
	$transcripts{$var}{$feature}{'Exon_number "' . $attributes{'EXON'} . '"'} = 1 if exists $attributes{'EXON'};
	
	$transcripts{$var}{$feature}{"cDNA_position \"$cdna_pos\""} = 1 unless $cdna_pos eq '-';
	$transcripts{$var}{$feature}{"CDS_position \"$cds_pos\""} = 1 unless $cds_pos eq '-';
	$transcripts{$var}{$feature}{"Protein_position \"$prot_pos\""} = 1 unless $prot_pos eq '-';
	$transcripts{$var}{$feature}{"Codon_change \"$codons\""} = 1 unless $codons eq '-';
	$transcripts{$var}{$feature}{"Amino_acid_change \"$aas\""} = 1 unless $aas eq '-';

	$hgvsg{$var} = $attributes{'HGVSg'};

	for my $pp ('SIFT', 'PolyPhen') {
	    next unless exists $attributes{$pp};
	    my ($prediction, $score) = $attributes{$pp} =~ /^(.+)\((.+)\)$/;
	    $transcripts{$var}{$feature}{"$pp $score \"$prediction\""} = 1;
	}

    }
    close (VEP);

    return (\%transcripts, \%hgvsg);
}


=head2 restructure_hash

    Title:    restructure_hash
    Function: converts ref of hash of hashes to ref of hash of arrays
    Args:     ref to hash of hashes
    Returns:  ref to hash of arrays

=cut

sub restructure_hash{
    my $old = shift;

    my %new;
    for my $key (keys %$old) {
	push @{$new{$key}}, keys %{$old->{$key}};
    }

    return \%new;
}


1;
