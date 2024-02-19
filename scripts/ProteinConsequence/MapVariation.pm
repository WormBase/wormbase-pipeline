=head1 NAME

ProteinConsequence::MapVariation - map variations and create VCF files

=cut

=head1 DESCRIPTION

Maps variations using their flanking sequence using map_Alleles package.
Filters out unsuitable and strangely mapped variants.
Creates VCF files for input to VEP.
Dumps mapped allele data to file.

=cut

package ProteinConsequence::MapVariation;

use strict;

use File::Path qw(make_path remove_tree);
use Data::Dumper;
use Wormbase;
use Ace;
use map_Alleles;
use Bio::Seq;
use Storable 'store';

use base ('ProteinConsequence::BaseProteinConsequence');


sub run {
    my $self = shift;

    $self->param('mapping_failure',1);
    my $batch_id = $self->required_param('batch_id');
    my @variation_ids = split(',', $self->required_param('variation_id_string'));
			      
    my $working_dir = $self->required_param('output_dir') . "/$batch_id/";
    $self->param('working_dir', $working_dir);
    
    unless (-d $working_dir) {
	my $err;
	make_path($working_dir, {error => \$err});
	die "make_path failed: ".Dumper($err) if $err && @$err;
    }
    
    chdir $working_dir or die "Failed to chdir to $working_dir";
    
    $self->dbc->disconnect_when_inactive(1);

    my $wb = Wormbase->new(
	-autoace  => $self->required_param('database'),
	-organism => $self->required_param('species'),
	-debug    => $self->required_param('debug'),
	-test     => $self->required_param('test'),
	);
    my $ace = Ace->connect(-path => $wb->autoace) or die "Could not create AcePerl connection";

    MapAlleles::setup($self, $wb, $ace);

    $self->write_to("Fetching $batch_id alleles") if $self->required_param('debug');
    my $alleles = $self->get_alleles(\@variation_ids);
    $self->write_to("Filtering $batch_id alleles") if $self->required_param('debug');
    $alleles = MapAlleles::filter_alleles($alleles);
    $self->write_to("Mapping $batch_id alleles") if $self->required_param('debug');
    my $mapped_alleles = MapAlleles::map($alleles);

    $self->write_to("Removing insanely mapped alleles form $batch_id") if $self->required_param('debug');
    MapAlleles::remove_insanely_mapped_alleles($mapped_alleles);

    my $alleles_file = $working_dir . $batch_id . '_mapped_alleles.dmp';
    store $mapped_alleles, $alleles_file;

    $self->write_to("Extracting data for $batch_id VCF") if $self->required_param('debug');
    my $vcf_details = $self->extract_vcf_data($mapped_alleles);

    $self->write_to("Creating $batch_id VCF file") if $self->required_param('debug');
    $self->write_vcf($vcf_details, $batch_id);

    $self->dbc->disconnect_when_inactive(1);

    $self->param('mapping_failure', 0);
    $ace->close;
}


=head2 extract_vcf_data

    Title:    extract_vcf_data
    Function: retrieves data needed to construct VCF files from AceDB Variation objects
    Args:     hashref of mapped allele objects ({$variation_id} = $allele)
    Returns:  hashref of variation data ({$chr}{$variation_position}{$allele}{'ref'} = $ref,
                                         {$chr}{$variation_position}{$allele}{'alt'} = $alt)

=cut

sub extract_vcf_data {
    my ($self, $mapped_alleles) = @_;

    my %vcf_data;
    my $vcf_count = 0;
    while (my($k, $v) = each(%{$mapped_alleles})) {
	my $allele      = $v->{allele};
	my $public_name = $allele->Public_name;
	my $remark      = $allele->Remark;
	my $chr         = $v->{chromosome};
	my $start       = $v->{start};
	my $stop        = $v->{stop};
	my $strand      = $v->{orientation};

	($start, $stop) = ($stop, $start) if $start > $stop;

	if ($stop - $start >= $self->required_param('max_deletion_size')) {
	    $self->write_to(sprintf("WARNING: %s - Deletion bigger than maximum deletion ".
                             "size, skipping VEP analysis (%s ; Remark:%s)\n", 
                             $k,
                             $public_name,
                             $remark));
	    next;
	}

	my ($substitution, $deletion, $insertion);
	foreach my $tp ($allele->Type_of_mutation) {
	    my @list;
	    if ($tp->right) {
		push @list, lc($tp->right);
		if ($tp->right->right) {
		    push @list, lc($tp->right->right);
		}
	    }
	    if ($tp eq 'Substitution') {
		$substitution = \@list;
	    } elsif ($tp eq 'Insertion') {
		$insertion = $list[0];
	    } elsif ($tp eq 'Deletion') {
		$deletion = 1;
	    }
	}

	
	my ($ref, $alt);
	if (defined $substitution) {
	    if(defined $insertion or defined $deletion) {
		$self->write_to(sprintf("WARNING: %s - Substitution defined together with inDel, ".
                             "skipping VEP analysis (%s ; Remark:%s)\n", 
                             $k,
                             $public_name,
                             $remark));
		next;
	    }
	    if (scalar(@$substitution) != 2) {
		$self->write_to(sprintf("WARNING: %s - Substitution has missing/empty FROM and/or TO, ".
                             "skipping VEP analysis (%s; Remark:%s)\n", 
                             $k,
                             $public_name,
                             $remark));
		next;
	    }
	    ($ref, $alt) = @$substitution;
	    if ($strand eq '-') {
		$ref = Bio::Seq->new(-seq => $ref, -alphabet => 'dna')->revcom->seq;
		$alt = Bio::Seq->new(-seq => $alt, -alphabet => 'dna')->revcom->seq;
	    }
	}
	else {
	    if (!defined $deletion and !defined $insertion) {
		$self->write_to(sprintf("WARNING: %s - Mutation type not substitution, ".
                             "deletion, or insertion - Skipping VEP analysis (%s ; Remark:%s)\n", 
                             $k,
                             $public_name,
                             $remark));
		next;
	    }
	    if (defined $insertion and $insertion !~ /^[A-Za-z]+$/) {
		$self->write_to(sprintf("WARNING: %s - Non-standard insertion %s - ".
				       "Skipping VEP analysis (%s; Remark:%s)\n",
				       $k,
				       $insertion,
				       $public_name,
				       $remark));
		next;
	    }
	    my $inserted_sequence = '';
	    my $deleted_sequence = '';
	    if (defined $insertion) {
		$inserted_sequence = $strand eq '+' ? $insertion :
		    Bio::Seq->new(-seq => $insertion, -alphabet => 'dna')->revcom->seq;
	    }
	    if (defined $deletion) {
		$deleted_sequence = $self->fetch_region($chr, $start, $stop);
	    }
	    my ($left_flank_seq, $right_flank_seq);
	    if (defined $allele->Flanking_sequences->name) {
		$left_flank_seq = $allele->Flanking_sequences->name;
	    }
	    if (defined $allele->Flanking_sequences->right->name) {
		$right_flank_seq = $allele->Flanking_sequences->right->name;
	    }
	    if ($strand eq '-') {
		($left_flank_seq, $right_flank_seq) = ($right_flank_seq, $left_flank_seq);
		$left_flank_seq = Bio::Seq->new(
		    -seq => $left_flank_seq,
		    -alphabet => 'dna')->revcom->seq if defined $left_flank_seq;
		$right_flank_seq = Bio::Seq->new(
		    -seq => $right_flank_seq,
		    -alphabet => 'dna')->revcom->seq if defined $right_flank_seq;
	    }
	    if ($start != 1 or (defined $insertion and !defined $deletion)) {
		$start-- unless defined $insertion and !defined $deletion; # Mapping starts and ends for insertions include flanking bases unless sequence is also deleted
		my $left_pad_base = defined $left_flank_seq ? substr($left_flank_seq, -1) :
		    $self->fetch_region($chr, $start, $start);
		$ref = $left_pad_base . $deleted_sequence;
		$alt = $left_pad_base . $inserted_sequence;
	    }
	    else {
		my $right_pad_base = defined $right_flank_seq ? substr($right_flank_seq, 0, 1) :
		    $self->fetch_region($chr, $stop + 1, $stop + 1);
		$ref = $deleted_sequence . $right_pad_base;
		$alt = $inserted_sequence . $right_pad_base;
	    }
	}
	if ($alt =~ /[^ACGT\-acgt]/) {
	    $alt = $self->get_possible_alleles($alt);
	}
	$vcf_count++;
	$vcf_data{$chr}{$start}{$allele->name}{'ref'} = $ref;
	$vcf_data{$chr}{$start}{$allele->name}{'alt'} = $alt;
    }

    $self->write_to("VCF_INFO: $vcf_count entries in " . $self->required_param('batch_id') . "VCF\n") if $self->required_param('debug');

    return \%vcf_data;
}


=head2 fetch_region

    Title:    fetch_region
    Function: retrieve genomic sequence from FASTA
    Args:     chromosome string, start position, end position
    Returns:  sequence string

=cut

sub fetch_region {
    my ($self, $chr, $start, $end) = @_;

    my $fasta = $self->required_param('fasta');
    my $seq = '';
    open (FAIDX, "samtools faidx $fasta $chr:$start-$end|");
    while (<FAIDX>) {
	chomp;
	next if $_ =~ /^>/;
	$seq .= $_;
    }
    close (FAIDX);

    return $seq;
}


=head2 get_possible_alleles

    Title:    get_possible_alleles
    Function: expand all IUPAC codes in alternative allele to give alternative alleles in
              a string as alt1,alt2,alt3...
    Args:     alternative allele string
    Returns:  reformatted alternative allele string

=cut

sub get_possible_alleles {
    my ($self, $iupac_alt) = @_;

    my %iupac_codes = (R => ['A', 'G'],
		       Y => ['C', 'T'],
		       S => ['G', 'C'],
		       W => ['A', 'T'],
		       K => ['G', 'T'],
		       M => ['A', 'C'],
		       B => ['C', 'G', 'T'],
		       D => ['A', 'G', 'T'],
		       H => ['A', 'C', 'T'],
		       V => ['A', 'C', 'G'],
	);

    my @possible_alleles = ('');
    for my $base_ix (0 .. length($iupac_alt) - 1) {
	my $base = uc(substr($iupac_alt, $base_ix, 1));
	if (exists $iupac_codes{$base}) {
	    my @new_possible_alleles = ();
	    for my $pa (@possible_alleles) {
		for my $pb (@{$iupac_codes{$base}}) {
		    push @new_possible_alleles, $pa . $pb;
		}
	    }
	    @possible_alleles = @new_possible_alleles;
	}
	else {
	    for my $ix (0 .. @possible_alleles - 1) {
		$possible_alleles[$ix] .= $base;
	    }
	}
    }

    return join(',', @possible_alleles);
}


=head2 write_vcf

    Title:    write_vcf
    Function: writes VCF file for VEP input
    Args:     hashref of variant data ({$chr}{$variation_position}{$allele}{'ref'} = $ref,
                                       {$chr}{$variation_position}{$allele}{'alt'} = $alt)
    Returns:  n/a
=cut

sub write_vcf {
    my ($self, $vcf_data, $batch_id) = @_;

    open (VCF, '> ' . $batch_id) or die "Could not open VCF for writing";

    for my $chr (sort keys %$vcf_data) {
	for my $pos (sort {$a<=>$b} keys %{$vcf_data->{$chr}}) {
	    for my $var (keys %{$vcf_data->{$chr}{$pos}}) {
		print VCF join("\t", $chr, $pos, $var, 
			       uc($vcf_data->{$chr}{$pos}{$var}{'ref'}),
			       uc($vcf_data->{$chr}{$pos}{$var}{'alt'}),
			       '.', '.', '.') . "\n";
	    }
	}
    }

    close (VCF);

    return;
}



=head2 get_alleles

    Title:    get_alleles
    Function: retrieve AceDB variation objects using map_Alleles package
    Args:     arrayref of variation IDs
    Returns:  arrayref of AceDB variation objects

=cut

sub get_alleles {
    my ($self, $variation_ids) = @_;

    my @alleles;
    for my $variation_id (@$variation_ids) {
	my $allele = MapAlleles::get_allele($variation_id);
	push @alleles, @$allele;
    }

    return \@alleles;
}


sub write_output {
    my $self = shift;

    if ($self->param('mapping_failure')) {
	$self->dataflow_output_id({batch_id => $self->param('batch_id')}, 3);
    }
    else {
	$self->dataflow_output_id({batch_id => $self->param('batch_id')}, 2);
    }
}


1;
