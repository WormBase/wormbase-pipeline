#!/usr/bin/env perl
#
# overload_gff_variation.pl
#
# Overloads Variation lines with extra info (consequence etc)
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2014-11-28 14:30:09 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use Ace;
use strict;
use Const::Fast;
use Bio::Seq;
use Bio::SeqIO;

const my %VCF_INFO => (
    status => {
	id          => 'ST',
	type        => 'String',
	number      => 1,
	description => 'Status',
    },
    polymorphism => {
	id          => 'PL',
	type        => 'Flag',
	number      => 0,
	description => 'Polymorphism',
    },
    production_method => {
	id          => 'PM',
	type        => 'String',
	number      => 1,
	description => 'Production method',
    },
    engineered => {
	id          => 'EN',
	type        => 'Flag',
	number      => 0,
	description => 'Engineereed',
    },
    public_name => {
	id          => 'PN',
	type        => 'String',
	number      => 1,
	description => 'Public name',
    },
    csq_string => {
	id          => 'CSQ',
	type        => 'String',
	number      => '.',
	description => 'Consequence annotations from EnsEMBL VEP. Format: Transcript|Gene|Consequence|Impact|SIFT_prediction|SIFT_score|PolyPhen_prediction|PolyPhen_score|Codon_change|HGVSp|HGVSc|Amino_acid_position|CDS_position|cDNA_position|Exon_nr|Intron_nr',
    },
    other_name => {
	id          => 'ON',
	type        => 'String',
	number      => 1,
	description => 'Other name',
    },
    strain => {
	id          => 'SN',
	type        => 'String',
	number      => 1,
	description => 'Strain',
    },
    rflp => {
	id          => 'RF',
	type        => 'Flag',
	number      => 0,
	description => 'RFLP',
    },
    hgvsg => {
	id          => 'HG',
	type        => 'String',
	number      => 1,
	description => 'HGVSg identifier',
    },
    );
    
######################################
# variables and command-line options # 
######################################
    
my ($debug, $test, $store, $wormbase, $database, $chromosomes);
my ($species, $gff3, $infile, $outfile, %var, $changed_lines, $added_lines, $no_vcf);
my ($fasta_file, $vcf_file, $vcf_only);

GetOptions (
  "debug=s"    => \$debug,
  "test"       => \$test,
  "store:s"    => \$store,
  "species:s"  => \$species,
  "database:s" => \$database,
  "gff3"       => \$gff3,
  "no_vcf"     => \$no_vcf,
  "vcf_only"   => \$vcf_only,
  "infile:s"   => \$infile,
  "outfile:s"  => \$outfile,
  "vcffile:s"  => \$vcf_file,
  "fasta:s"    => \$fasta_file,
    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( '-debug'   => $debug,
                             '-test'    => $test,
                             '-organism' => $species,
			     );
}

$database = $wormbase->autoace if not defined $database;
$species = $wormbase->species;
my $sp_full_name = $wormbase->full_name;

my $log = Log_files->make_build_log($wormbase);

if (not defined $infile or not defined $outfile) { 
  $log->log_and_die("You must define -infile and -outfile\n") unless $vcf_only;
}

my $vcf = $gff3 and !$no_vcf ? 1 : 0;

my $vcf_out_fh;
if ($vcf) {
    $fasta_file = $wormbase->genome_seq() if !$fasta_file;
    if (!$vcf_file) {
	($vcf_file) = $fasta_file =~ /^(.+)\.genome.fa$/;
	$vcf_file .= '.vcf';
    }
    $chromosomes = get_chromosome_sequences();
    open ($vcf_out_fh, ">$vcf_file") or $log->log_and_die("Could not open $vcf_file for writing\n");
    print_vcf_header($vcf_out_fh, $chromosomes, $sp_full_name);
}

my $transcript_parents = get_transcript_parents($infile);
my $var_consequences = get_molecular_consequences($transcript_parents);

my $db = Ace->connect(-path => $database);

my $gff_in_fh;
if ($vcf_only and $infile =~ /\.gz$/) {
    # May want to run VCF only mode on final gzipped GFF, otherwise will always be in unzipped form
    open ($gff_in_fh, "gunzip -c $infile |") or $log->log_and_die("Could not open $infile for reading\n");
}
else {
    open($gff_in_fh, $infile) or $log->log_and_die("Could not open $infile for reading\n");
}

my $gff_out_fh;
if (!$vcf_only) {
    open($gff_out_fh, ">$outfile") or $log->log_and_die("Could not open $outfile for writing\n");  
}

my %vcf_lines_printed;
while (<$gff_in_fh>) {
  if (/Variation \"(WBVar\d+)\"/ or 
      /Variation:(WBVar\d+)/ or 
      /variation=(WBVar\d+)/) {
    my $allele = $1; 
    
    my $is_putative_change_of_function_allele = 0;

    my @new_els;
    my %allele_vcf_values;
    my @current_els = split(/\t/, $_);
    pop @current_els;
    push @new_els, ['Variation', $allele];

    my $variation = $db->fetch(Variation => $allele);
    
    my @var_types = $variation->at('Variation_type');
    my @other_names = $variation->at('Name.Other_name');
    my ($public_name) = $variation->at('Name.Public_name');
    my ($prodmethod) = $variation->at('Origin.Production_method');
    my $natural_variant = 0;

    if ($public_name) {
      push @new_els, ['Public_name', $public_name];
      $allele_vcf_values{'public_name'} = $public_name;
    }
    if (@other_names) {
      push @new_els, ['Other_name', \@other_names];
      $allele_vcf_values{'other_name'} = join(',', @other_names);
    }
    
    if ($variation->Strain) {
      my @strains = map {$_->Public_name} $variation->Strain;
      push @new_els, ['Strain', \@strains];
      $allele_vcf_values{'strain'} = join(',', @strains);
    }
    if ($prodmethod) {
	push @new_els, ['Production_method', $prodmethod];
	$allele_vcf_values{'production_method'} = $prodmethod;
    }
    
    foreach my $tp (@var_types) {
	if ($tp =~ /confirmed/i) {
	    push @new_els, ['Status', 'Confirmed'];
	    $allele_vcf_values{'status'} = 'Confirmed';
	}
	if ($tp =~ /RFLP/) {
	    push @new_els, ['RFLP'];
	    $allele_vcf_values{'rflp'} = 1;
	}
	if ($tp =~ /Engineered/) {
	    push @new_els, ['Engineered'];
	    $allele_vcf_values{'engineered'} = 1;
	}

      if ($tp eq 'Natural_variant') {
        $natural_variant = 1;
        if ($current_els[2] ne 'transposable_element_insertion_site' and 
            $current_els[1] ne 'Polymorphism') {
          $current_els[1] .= "_Polymorphism";
        }
        push @new_els, ['Polymorphism'];
	$allele_vcf_values{'polymorphism'} = 1;
      }
    }
    
    my @types = $variation->at('Sequence_details.Type_of_mutation');
    my @subs = ();
    foreach my $type (@types) {
      if ($type eq 'Substitution' and defined $variation->Substitution) {
        @subs = $variation->Substitution->row;
        #print NEW " ; Substitution \"$substitution[0]/$substitution[1]\"";
        push @new_els, ['Substitution', join("/", uc($subs[0]), uc($subs[1]))];
      } elsif ($type eq 'Insertion' and defined $variation->Insertion) {
        my ($insert) = $variation->Insertion->row;
        push @new_els, ['Insertion', uc($insert)];
      }
    }

    if ($current_els[2] eq 'sequence_alteration' or $current_els[2] eq 'substitution') {
      # attempt to assign a more specific term here. 
      # "sequence_alteration" is used in the raw dumps for vars from a single project (e.g. MMP) 
      # which can be one of several different types; 
      # "substition" is used specfifically for substitution alleles, but we want to use 
      # the more specfic "point_mutation" or "SNP" where appropriate

      my $new_term = "sequence_alteration";
      if (scalar(@types) == 1) {
        my ($tp) = @types;
        
        if ($tp eq 'Insertion' or $tp eq 'Deletion' or $tp eq 'Inversion') {
          $new_term = lc($tp);
          $new_term .= "_site" if $new_term eq 'insertion';
        } elsif ($tp eq 'Substitution') {
          $new_term = "substitution";
          
          if ($current_els[3] == $current_els[4] and 
              scalar(@subs) >= 2 and
              length($subs[0]) == 1 and 
              length($subs[1]) == 1) {
            # single nucleotide
            if ($natural_variant) {
              $new_term = 'SNP';
            } else {
              $new_term = 'point_mutation';
            }
          }
        } elsif ($tp eq 'Tandem_duplication') {
          $new_term = 'tandem_duplication';
        }
      } elsif (scalar(@types) > 1) {
        if (grep { $_ eq 'Tandem_duplication' } @types) {
          # some tandem duplications are assoiated with micro-insertions/deletions at the flanks, 
          # so will have multiple Type_of_mutation. We still want to classify them as tandem_duplications
          $new_term = "tandem_duplication";
        } else {
          # more than one type, so just put complex_substitution
          $new_term = "complex_substitution";
        }
      }
      
      $current_els[2] = $new_term;
    }
    
    if (exists $var_consequences->{$allele}) {
	for my $attribute ('Consequence', 'VEP_impact', 'AAchange', 'Codon_change', 'HGVSg', 'HGVSc',
			   'SIFT', 'PolyPhen', 'cDNA_position', 'CDS_position', 'AA_position',
			   'Intron_nr', 'Exon_nr') {
	    push @new_els, [$attribute, $var_consequences->{$allele}{$attribute}]
		if exists $var_consequences->{$allele}{$attribute};
	}
	$allele_vcf_values{'csq_string'} = join(',', @{$var_consequences->{$allele}{'csq'}})
	    if exists $var_consequences->{$allele}{'csq'};
	
	if ($var_consequences->{$allele}{'severity'} >= 23) {
	    if ($current_els[2] ne 'transposable_element_insertion_site' and 
		$current_els[2] ne 'tandem_duplication') {
		$is_putative_change_of_function_allele = 1;
	    }
	}
	
    }
    
    my @new_el_strings;
    foreach my $el (@new_els) {
      if (scalar(@$el) == 1) {
        if ($gff3) {
          push @new_el_strings, sprintf("%s=1", lc($el->[0]));
        } else {
          push @new_el_strings, $el->[0];
        }
      } else {
        my ($k, $v) = @$el;
        if ($gff3) {
          if (ref($v) eq 'ARRAY') {
            $v = join(",", @$v);
          } 
          push @new_el_strings, join("=", lc($k), $v);
        } else {
          if (ref($v) eq 'ARRAY') {
            foreach my $vv (@$v) {
              push @new_el_strings, sprintf("%s \"%s\"", $k, $vv);
            }
          } else {
            push @new_el_strings, sprintf("%s \"%s\"", $k, $v);
          }
        }
      }
    }
    
    my $join_str = ($gff3) ? ";" : " ; ";
    my $group = join($join_str, @new_el_strings); 
    print $gff_out_fh join("\t", @current_els, $group), "\n" unless $vcf_only;

    if ($vcf and $current_els[2] ne 'tandem_duplication'
	and $current_els[2] ne 'transposable_element_insertion_site'
	and $current_els[2] ne 'sequence_alteration') {
	next if exists $vcf_lines_printed{$allele};
	my ($vcf_pos, $vcf_ref, $vcf_alt) = vcf_values($allele, \@current_els, $group, $chromosomes);
	if ($vcf_pos) {
	    my $vcf_info_str = vcf_info(\%allele_vcf_values);
	    $vcf_out_fh->print("\n" . join("\t", $current_els[0], $vcf_pos, $allele, $vcf_ref,
					   $vcf_alt, '.', '.', $vcf_info_str));
	    $vcf_lines_printed{$allele} = 1;
	}
    }

    $changed_lines++;
        
    # duplicate putative change-of-function alleles, with a different source, so that they
    # can be drawn in a separate track
    if ($is_putative_change_of_function_allele) {
      $current_els[1] = "PCoF_" . $current_els[1];
      print $gff_out_fh join("\t", @current_els, $group), "\n" unless $vcf_only;
      $added_lines++;
    }
    
    if ($changed_lines % 50000 == 0) {
      $db->close();
      $db = Ace->connect(-path => $database);
    }
  } else {
    print $gff_out_fh "$_" unless $vcf_only;
  }
}
close($gff_out_fh) or $log->log_and_die("Could not close $outfile after writing\n") unless $vcf_only;
if ($vcf) {
    close($vcf_out_fh) or $log->log_and_die("Could not close $vcf_file after writing\n");
}
$db->close();

$log->write_to("Finished processing : $changed_lines lines modified, $added_lines added\n") unless $vcf_only;
$log->mail();
exit(0);


##############################################################
#
# Subroutines
#
##############################################################

sub get_chromosome_sequences {
    my %chromosomes;
    my $chr;
    my $fasta = Bio::SeqIO->new(-file => $fasta_file, -format => 'fasta')
	or $log->log_and_die("Could not read FASTA file $fasta_file");
    while (my $seq = $fasta->next_seq) {
	my $chr = $seq->display_id;
	$chr =~ s/CHROMOSOME_//;
	$chr =~ s/^chr//;
	$chromosomes{$chr} = $seq->seq;
    }

    return \%chromosomes;
}


sub get_molecular_consequences {
    my $transcript_parents = shift;
    
    # ranking based on EnsEMBL VEP order of severity
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
    
    my $table = $wormbase->table_maker_query($database, &write_mol_change_def_file);
  
    my %var_consequences;
    while(<$table>) {
	chomp;
	s/\"//g; 
	s/\\//g;
	next if (! defined $_);
	next if (/acedb/ or /\/\//);
	
	my ($var_name, $transcript, $consequence_string, $vep_impact, $aa_change, $codon_change,
	    $sift_score, $sift_prediction, $polyphen_score, $polyphen_prediction, $hgvsg, $hgvsc,
	    $hgvsp, $cdna_pos, $cds_pos, $prot_pos, $intron_nr, $exon_nr) = split(/\t/, $_);
	my @consequences = split(',', $consequence_string);
	
	for my $consequence (@consequences) {
	    my $severity = $severity_ranking{$consequence};
	    push @{$var_consequences{$var_name}{'csq'}}, join('|', $transcript,
							      $transcript_parents->{$transcript},
							      $consequence, $vep_impact, $sift_prediction,
							      $sift_score, $polyphen_prediction,
							      $polyphen_score, $codon_change, $hgvsp,
							      $hgvsc, $prot_pos, $cds_pos, 
							      $cdna_pos, $exon_nr, $intron_nr);

	    next if exists $var_consequences{$var_name} and $var_consequences{$var_name}{'severity'} > $severity;
	    $var_consequences{$var_name}{'Consequence'} = $consequence;
	    $var_consequences{$var_name}{'severity'} = $severity;
	    $var_consequences{$var_name}{'VEP_impact'} = $vep_impact;
	    $var_consequences{$var_name}{'AAchange'} = $aa_change if $aa_change =~ /\//;
	    $var_consequences{$var_name}{'Codon_change'} = $codon_change if $codon_change =~ /\//;
	    $var_consequences{$var_name}{'HGVSg'} = $hgvsg; 
	    $var_consequences{$var_name}{'HGVSc'} = $hgvsc if $hgvsc;
	    $var_consequences{$var_name}{'HGVSp'} = $hgvsp if $hgvsp;
	    $var_consequences{$var_name}{'SIFT'} = $sift_prediction . '(' . $sift_score . ')' if $sift_score;
	    $var_consequences{$var_name}{'PolyPhen'} = $polyphen_prediction . '(' . $polyphen_score . ')' if $polyphen_score;
	    $var_consequences{$var_name}{'cDNA_position'} = $cdna_pos if $cdna_pos;
	    $var_consequences{$var_name}{'CDS_position'} = $cds_pos if $cds_pos;
	    $var_consequences{$var_name}{'AA_position'} = $prot_pos if $prot_pos;
	    $var_consequences{$var_name}{'Intron_nr'} = $intron_nr if $intron_nr;
	    $var_consequences{$var_name}{'Exon_nr'} = $exon_nr if $exon_nr;
	    
	}
    }
    close($table);

    return \%var_consequences;
}


sub get_transcript_parents {
    my $file = shift;

    my %transcript_parents;
    open (GFF, "<$file") or $log->log_and_die("Can't open $file for reading\n");
    while (<GFF>) {
	next unless $_ =~ /Parent=Gene:([^;]+)/;
	my $parent = $1;
	next unless $_ =~ /ID=Transcript:([^;]+)/;
	$transcript_parents{$1} = $parent;
    }
    close (GFF);

    return \%transcript_parents;
}


sub print_vcf_header {
    my ($out_fh, $chromosomes, $species) = @_;

    $out_fh->print('##fileformat=VCFv4.3');
    for my $chr (keys %$chromosomes) {
	$out_fh->print("\n##contig=<ID=${chr},length=" . length($chromosomes->{$chr}) . ",species=\"${species}\">");
    }
    for my $field (keys %VCF_INFO) {
	$out_fh->print("\n##INFO=<ID=" . $VCF_INFO{$field}{'id'} . ',Number=' .
		       $VCF_INFO{$field}{'number'} . ',Type=' . $VCF_INFO{$field}{'type'} .
		       ',Description="' . $VCF_INFO{$field}{'description'} . '">');
    }
    $out_fh->print("\n" . join("\t", qw(#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT)));

    return;
}

sub reverse_complement {
    my $seq = shift;

    my $seq_obj = Bio::Seq->new(-seq => $seq, -alphabet => 'dna');
    return $seq_obj->revcom->seq;
}

sub vcf_info {
    my $values = shift;

    my @entries;
    for my $attr (keys %VCF_INFO) {
	next unless exists $values->{$attr};
	my $entry = $VCF_INFO{$attr}{'id'};
	unless ($VCF_INFO{$attr}{'type'} eq 'Flag') {
	    $entry .= '=' . $values->{$attr};
	}
	push @entries, $entry;
    }

    return join(';', @entries);
}


sub vcf_values {
    my ($allele, $gff_columns, $attributes, $chromosomes) = @_;

    my $chr = $gff_columns->[0];
    my $start = $gff_columns->[3];
    my $end = $gff_columns->[4];
    my $type = $gff_columns->[2];

    $chr =~ s/CHROMOSOME_//;
    $chr =~ s/^chr//;

    my ($ref, $alt, $pos, $substitution);
    if ($attributes =~ /insertion=([^;]+)/) {
	$alt = $1;
    }
    elsif ($attributes =~ /substitution=([^;]+)/) {
	$substitution = $1;
    }

    unless ($type eq 'insertion_site') {
	$ref = substr($chromosomes->{$chr}, $start - 1, ($end - $start) + 1);
    }

    if ($type eq 'insertion_site' or $type eq 'deletion' or $type eq 'complex_substitution') {
	my $padding_base = substr($chromosomes->{$chr}, $start - 2, 1);
	$pos = $start - 1;
	$ref = $padding_base . $ref;
	$alt = $padding_base . $alt;
    }
    else {
	if (!$substitution) {
	    $log->write_to("INFO: No substitution found for $type $allele, not including in VCF\n");
	    return;
	}
	
	$pos = $start;
	if ($substitution =~ /^$ref\/(.+)$/) {
	    $alt = $1;
	}
	else {
	    my $revref = reverse_complement($ref);
	    if ($substitution =~ /^$revref\/(.+)$/) {
		$alt = reverse_complement($1);
	    }
	    else {
		$log->write_to("INFO: Reference sequence at position $chr:$start-$end ($ref) " .
			       "doesn't match subsitution in GFF ($substitution), not " . 
			       "including $allele in VCF\n");
		return;
	    }
	}
    }

    return ($pos, $ref, $alt);
}


sub write_mol_change_def_file {
  my $def = '/tmp/overload_SNP_GFF_mol_chng.def';
  open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $txt = <<END;

Sortcolumn 1

Colonne 1
Width 12
Optional
Visible
Class
Class Variation
From 1

Colonne 2
Width 12
Mandatory
Hidden
Show_Tag
From 1
Tag Sequence

Colonne 3
Width 12
Mandatory
Visible
Class
Class Transcript
From 1
Tag Transcript

Colonne 4
Width 12
Optional
Hidden
Show_Tag
Right_of 3
Tag Molecular_change # VEP_consequence                                                                                                                                                                                                                              
Colonne 5
Width 12
Optional
Visible
Text
Right_of 4
Tag  HERE

Colonne 6
Width 12
Optional
Hidden
Show_Tag
Right_of 3
Tag Molecular_change # VEP_impact                                                                                                                                                                                                                            
Colonne 7
Width 12
Optional
Visible
Text
Right_of 6
Tag  HERE

Colonne 8
Width 12
Optional
Hidden
Show_Tag
Right_of 3
Tag Molecular_change # Amino_acid_change                                                                   

Colonne 9
Width 12
Optional
Visible
Text
Right_of 8
Tag  HERE

Colonne 10
Width 12
Optional
Hidden
Show_Tag
Right_of 3
Tag Molecular_change # Codon_change                                                          
                                                                                                                                                                                                    
Colonne 11
Width 12
Optional
Visible
Text
Right_of 10
Tag  HERE

Colonne 12
Width 12
Optional
Hidden
Show_Tag
Right_of 3
Tag Molecular_change # SIFT 
                                                                                                                                                               
Colonne 13
Width 12
Optional
Visible
Float
Right_of 12
Tag  HERE

Colonne 14
Width 12
Optional
Visible
Text
Right_of 13
Tag  HERE

Colonne 15 
Width 12
Optional
Hidden
Show_Tag
Right_of 3
Tag Molecular_change # PolyPhen
                                                                                                                                                                                                    
Colonne 16
Width 12
Optional
Visible
Float
Right_of 15
Tag  HERE
                                                                                                                                                                                                    
Colonne 17
Width 12
Optional
Visible
Text
Right_of 16
Tag  HERE

Colonne 18
Width 12
Mandatory
Visible
Text
From 1
Tag HGVSg

Colonne 19
Width 12
Optional
Hidden
Show_Tag
Right_of 3
Tag Molecular_change # HGVSc 
                                                                                                                                                                                                    
Colonne 20
Width 12
Optional
Visible
Text
Right_of 19
Tag  HERE

Colonne 21
Width 12
Optional
Hidden
Show_Tag
Right_of 3
Tag Molecular_change # HGVSp
                                                                                                                                                                                                    
Colonne 22
Width 12
Optional
Visible
Text
Right_of 21
Tag  HERE

Colonne 23
Width 12
Optional
Hidden
Show_Tag
Right_of 3
Tag Molecular_change # cDNA_position
                                                                                                                                                                                                    
Colonne 24
Width 12
Optional
Visible
Text
Right_of 23
Tag  HERE

Colonne 25
Width 12
Optional
Hidden
Show_Tag
Right_of 3
Tag Molecular_change # CDS_position                                                         
                                                                                                                                                                                                    
Colonne 26
Width 12
Optional
Visible
Text
Right_of 25
Tag  HERE

Colonne 27
Width 12
Optional
Hidden
Show_Tag
Right_of 3
Tag Molecular_change # Protein_position                                                          
                                                                                                                                 
Colonne 28
Width 12
Optional
Visible
Text
Right_of 27
Tag  HERE

Colonne 29
Width 12
Optional
Hidden
Show_Tag
Right_of 3
Tag Molecular_change # Intron_number                                                                                                    

Colonne 30
Width 12
Optional
Visible
Text
Right_of 29
Tag  HERE

Colonne 31
Width 12
Optional
Hidden
Show_Tag
Right_of 3
Tag Molecular_change # Exon_number                                                                      

Colonne 32
Width 12
Optional
Visible
Text
Right_of 31
Tag  HERE

END

  print TMP $txt;
  close TMP;

  return $def;
}
