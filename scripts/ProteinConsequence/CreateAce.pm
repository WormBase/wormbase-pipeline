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

    my $transcripts = $self->get_transcript_consequences($vep_output);

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


    $self->write_to("Writing transcript details for $batch_id\n") if $self->required_param('debug');
    for my $var (keys %$transcripts) {
	next unless $self->add_transcript_consequences_for_variation($var, $mapped_alleles->{$var}{allele});
	print $fh "Variation : \"$var\"\n";
	for my $transcript (keys %{$transcripts->{$var}}) {
	    for my $tag (keys %{$transcripts->{$var}{$transcript}}) {
		print $fh "Transcript $transcript $tag\n";
	    }
	}
	print $fh "\n";
    }

    
    $self->write_to("Writing pseudogene details for $batch_id\n") if $self->required_param('debug');
    for my $var (keys %$pseudogenes) {
	print $fh "Variation : \"$var\"\n";
	for my $pseudogene (@{$pseudogenes->{$var}}) {
	    print $fh "Pseudogene $pseudogene\n";
	}
	print $fh "\n";
    }

    close($fh) or $self->log_and_die("Could not close $ace_file after writing\n");
    $ace->close;

}


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

    return (restructure_hash(\%gene_alleles), restructure_hash(\%allele_genes), $transcripts, restructure_hash(\%pseudogenes));
}


sub get_transcript_consequences {
    my ($self, $vep_output) = @_;

    my %transcripts;
    open (VEP, '<', $vep_output);
    while (<VEP>) {
	next if $_ =~ /^#/;
	chomp;
	my ($var, $pos, $allele, $gene, $feature, $feature_type, $consequence, $cdna_pos, $cds_pos,
	    $prot_pos, $aas, $codons, $existing_var, $attr) = split("\t", $_);

	next unless $feature_type eq 'Transcript';

	my $types = $self->get_types_from_vep_consequence($consequence, $feature, $cds_pos,
							  $prot_pos, $aas, $codons); # To be deprecated after WS282
	for my $type (keys %$types) {                                                # 
	    $transcripts{$var}{$feature}{$type}++;                                   #
	}                                                                            #
	
	my %attributes = split /[;=]/, $attr;
	$transcripts{$var}{$feature}{"VEP_consequence \"$consequence\""}++;
	
	$transcripts{$var}{$feature}{'VEP_impact "' . $attributes{'IMPACT'} . '"'} = 1 if exists $attributes{'IMPACT'};
	$transcripts{$var}{$feature}{'HGVSc "' . $attributes{'HGVSc'} . '"'} = 1 if exists $attributes{'HGVSc'};
	$transcripts{$var}{$feature}{'HGVSp "' . $attributes{'HGVSp'} . '"'} = 1 if exists $attributes{'HGVSp'};
	$transcripts{$var}{$feature}{'Intron_number "' . $attributes{'INTRON'} . '"'} = 1 if exists $attributes{'INTRON'};
	$transcripts{$var}{$feature}{'Exon_number "' . $attributes{'EXON'} . '"'} = 1 if exists $attributes{'EXON'};

	$transcripts{$var}{$feature}{"cDNA_position \"$cdna_pos\""} = 1 unless $cdna_pos eq '-';
	$transcripts{$var}{$feature}{"CDS_position \"$cds_pos\""} = 1 unless $cds_pos eq '-';
	$transcripts{$var}{$feature}{"Protein_position \"$prot_pos\""} = 1 unless $prot_pos eq '-';
	$transcripts{$var}{$feature}{"Codon_change \"$codons\""} = 1 unless $codons eq '-';
	$transcripts{$var}{$feature}{"Amino_acid_change \"$aas\""} = 1 unless $aas eq '-';

	for my $pp ('SIFT', 'PolyPhen') {
	    next unless exists $attributes{$pp};
	    my ($prediction, $score) = $attributes{$pp} =~ /^(.+)\((.+)\)$/;
	    $transcripts{$var}{$feature}{"$pp $score \"$prediction\""} = 1;
	}

    }
    close (VEP);

    return \%transcripts;
}


sub get_types_from_vep_consequence { # To be deprecated after WS282
    my ($self, $consequence_string, $feature, $cds_pos, $prot_pos, $aas, $codons) = @_;

    my %types;
    my $aa_string = $aas;
    $aa_string =~ s/\// to /;
    my @consequences = split(',', $consequence_string);
    for my $consequence (@consequences) {
	$types{'Intron'} = 1 if $consequence eq 'splice_donor_variant' or $consequence eq 'splice_acceptor_variant'
	     or $consequence eq 'non_coding_transcript_variant';
	$types{'Donor'} = 1 if $consequence eq 'splice_donor_variant';
	$types{'Acceptor'} = 1 if $consequence eq 'splice_acceptor_variant';
	$types{"Readthrough \"$aa_string\""} = 1 if $consequence eq 'stop_lost'; 
	$types{'Coding_exon'} = 1 if $consequence eq 'stop_gained' or $consequence eq 'stop_lost' or $consequence eq 'frameshift_variant'
	    or $consequence eq 'start_lost' or $consequence eq 'inframe_insertion' or $consequence eq 'inframe_deletion'
	    or $consequence eq 'protein_altering_variant' or $consequence eq 'coding_sequence_variant';
	$types{"Silent \"$aas ($prot_pos)\""} = 1 if $consequence eq 'synonymous_variant';
	$types{'UTR_5'} = 1 if $consequence eq '5_prime_UTR_variant';
	$types{'UTR_3'} = 1 if $consequence eq '3_prime_UTR_variant';
	$types{'Noncoding_exon'} = 1 if $consequence eq 'non_coding_transcript_exon_variant';

	if ($consequence eq 'missense_variant' ){
	    $prot_pos = substr($prot_pos, 0, index($prot_pos,'-')) if $prot_pos =~ /\-/;
	    $types{"Missense $prot_pos \"$aa_string\""} = 1;
	}
	elsif ($consequence eq 'stop_gained') {
	    my $var_codon = substr($codons, index($codons, '/') + 1);
	    my ($original_aa, $var_aa) = $aas =~ /^([^\.]+)\/(.+)$/;
	    my $stop_ix = index($var_aa, '*');
	    my $stop_codon = substr($var_codon, 3 * $stop_ix, 3);
	    my ($stop_tag, $stop_colour);
	    if (uc($stop_codon) eq 'TAR') {
		$stop_tag = 'Amber_UAG_or_Ochre_UAA';
		$stop_colour = 'amber or ochre';
	    }
	    elsif (uc($stop_codon) eq 'TRA') {
		$stop_tag = 'Ochre_UAA_or_Opal_UGA';
		$stop_colour = 'ochre or opal';
	    }
	    elsif (uc($stop_codon) =~ /[TWYKHDB]AG/ ) {
		$stop_tag = 'Amber_UAG';
		$stop_colour = 'amber';
	    }
	    elsif (uc($stop_codon) =~ /[TWYKHDB]AA/) {
		$stop_tag = 'Ochre_UAA';
		$stop_colour = 'ochre';
	    }
	    elsif (uc($stop_codon) =~ /[TWYKHDB]GA/) {
		$stop_tag = 'Opal_UGA';
		$stop_colour = 'opal';
	    }
	    else {
		$self->write_to("ERROR: $feature $stop_codon is not Amber/Opal/Ochre");
	    }
	    
	    if (defined $stop_tag) {
		my $text = length $original_aa < length $var_aa ? "$stop_colour stop inserted" :
		    "$original_aa to $stop_colour stop";
		$types{"Nonsense $stop_tag \"$text ($prot_pos)\""} = 1;
	    }
	}
	elsif ($consequence eq 'frameshift_variant') {
	    my ($wt, $allele) = split('/', $codons);
	    my $bp_change = length $allele - length $wt;
	    my $indel_type = $bp_change < 0 ? 'deletion' : 'insertion';
	    $types{'Frameshift "' . abs($bp_change) . "bp $indel_type at CDS position $cds_pos\""} = 1;
	}
	elsif ($consequence eq 'start_lost') {
	    $types{"Missense $prot_pos \"$aa_string\""} = 1 unless length $aas > 3;
	}
	
    }
    
    return \%types;
}


sub restructure_hash{
    my $old = shift;

    my %new;
    for my $key (keys %$old) {
	push @{$new{$key}}, keys %{$old->{$key}};
    }

    return \%new;
}


1;