#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;
use Path::Class;
use Const::Fast;

use lib $ENV{CVS_DIR};
use Wormbase;
use Modules::AGR;

my ($debug, $test, $verbose, $store, $wormbase, $curation_test, $limit, $schema);
my ($outdir, $acedbpath, $ws_version, $bgi_json,$disease_file);

const my $LINKML_SCHEMA => 'v1.10.0';
# TODO: check SO term mappings with Stavros
const my %TERM2SO => (
    'Insertion'          => 'SO:0000667', # changed insertion_site to insertion
    'Deletion'           => 'SO:0000159',
    'Point_mutation'     => 'SO:1000008',
    'Substitution'       => 'SO:1000002', # to hack our wrong SO-terms in => DELIN
    'Delins'             => 'SO:1000032', # to hack our wrong SO-terms in => DELIN
    'SNP'                => 'SO:0000694',
    'Tandem_duplication' => 'SO:1000173',
    'Natural_variant'    => 'SO:0001147',
    'RFLP'               => 'SO:0000412',
    'Engineered_allele'  => 'SO:0000783',
    'Transposon_insertion' => 'SO:0001054',
    'Allele'             => 'SO:0001023'
);

GetOptions (
    'debug=s'       => \$debug,
    'test'          => \$test,
    'verbose'       => \$verbose,
    'store:s'       => \$store,
    'database:s'    => \$acedbpath,
    'outdir:s'     => \$outdir,
    'curationtest' => \$curation_test,  
    'wsversion=s'   => \$ws_version,
    'limit:i'       => \$limit
)||die();

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
      );
}


my $tace = $wormbase->tace;
my $date = AGR::get_rfc_date();
my $alt_date = join("/", $date =~ /^(\d{4})(\d{2})(\d{2})/);
my $full_name = $wormbase->full_name;

$acedbpath  ||= $wormbase->autoace;
$ws_version ||= $wormbase->get_wormbase_version_name;
if (!$outdir) {
    $outdir = '.';
}
my $suffix = $curation_test ? '.test_data.json' : '.json';
my $allele_outfile = $outdir . "/wormbase.alleles.${ws_version}.${LINKML_SCHEMA}" . $suffix;
my $assoc_outfile = $outdir . "/wormbase.allele_associations.${ws_version}.${LINKML_SCHEMA}" . $suffix;
my $gene_assoc_outfile = $outdir . "/wormbase.allele_associations.${ws_version}.${LINKML_SCHEMA}" . $suffix . '.tmp';
my $variant_outfile = $outdir . "/wormbase.variants.${ws_version}.${LINKML_SCHEMA}" . $suffix;

my $db = Ace->connect(-path => $acedbpath, -program => $tace) or die("Connection failure: ". Ace->error);

my $allele_out_fh  = file($allele_outfile)->openw;
my $assoc_out_fh = file($assoc_outfile)->openw;
my $gene_assoc_out_fh = file($gene_assoc_outfile)->openw;
my $variant_out_fh = file($variant_outfile)->openw;
$allele_out_fh->print("{\n   \"linkml_version\" : \"" . $LINKML_SCHEMA . "\",\n    \"alliance_member_release_version\" : \"" . $ws_version . "\"\n    \"allele_ingest_set\" : [");
$assoc_out_fh->print("{\n   \"linkml_version\" : \"" . $LINKML_SCHEMA . "\",\n    \"alliance_member_release_version\" : \"" . $ws_version . "\"\n");
$gene_assoc_out_fh->print("    \"allele_gene_association_ingest_set\" : [");
$variant_out_fh->print("{\n    \"linkml_version\" : \"" . $LINKML_SCHEMA . "\",\n    \"alliance_member_release_version\" : \"" . $ws_version . "\"\n    \"variant_ingest_set\" : [");

my $gene_assoc_count = 0;
my $alleles = process_variations();

$alleles = process_transgenes();
$allele_out_fh->print("\n   ]\n}");
$variant_out_fh->print("\n   ]\n}");
$gene_assoc_out_fh->print("\n    ]");

my $gene_assoc_in_fh = file($gene_assoc_outfile)->openr;
while (my $line = $gene_assoc_in_fh->getline()) {
    chomp $line;
    $assoc_out_fh->print($line . "\n");
}
$assoc_out_fh->print("\n}\n");
$db->close;
system("rm $gene_assoc_outfile");

sub process_variations {
    my @alleles = $db->fetch(-query => "Find Variation WHERE Species = \"Caenorhabditis elegans\"");
    
    my $var_count = 0;
    for my $obj (@alleles) {
	$var_count++;
	next unless $obj->isObject();
	unless ($obj->Species) {
	    print "No species for $obj - skipping\n";
	    next;
	}

	my $dp_xref_dto_json = {
	    referenced_curie => 'WB:' . $obj->name,
	    page_area => 'allele',
	    display_name => $obj->name,
	    prefix => 'WB',
	    internal => JSON::false,
	    obsolete => JSON::false
	};

	my $data_provider_dto_json = {
	    source_organization_abbreviation => 'WB',
	    cross_reference_dto => $dp_xref_dto_json,
	    internal => JSON::false,
	    obsolete => JSON::false
	}; 

	if ($obj->Gene) {
	    for my $gene($obj->Gene) {
		$assoc_out_fh->print(",") if $gene_assoc_count > 0;
		$assoc_out_fh->print("{\"allele_curie\": \"WB:" . $obj->name . "\", \"relation_name\": \"is_allele_of\", \"gene_curie\": \"WB:" . $gene->name . "\"}");
		$gene_assoc_count++;
	    }
	}

	my $is_obsolete = ($obj->Status && $obj->Status->name eq 'Dead') ? JSON::true : JSON::false;
	    
	my $allele_obj = {
	    curie                      => "WB:" . $obj->name,
	    allele_symbol_dto          => get_symbol_dto($obj),
	    taxon_curie                => "NCBITaxon:" . $obj->Species->NCBITaxonomyID->name,
	    internal                   => JSON::false,
	    obsolete                   => $is_obsolete,
	    created_by_curie           => 'WB:curator',
	    updated_by_curie           => 'WB:curator',
	    data_provider_dto          => $data_provider_dto_json
	};

	my $var_obj = {
	    curie                      => "WBVar:" . $obj->name,
	    taxon_curie                => "NCBITaxon:" . $obj->Species->NCBITaxonomyID->name,
	    internal                   => JSON::false,
	    obsbolete                  => $is_obsolete,
	    created_by_curie           => 'WB:curator',
	    updated_by_curie           => 'WB:curator',
	    data_provider_dto          => $data_provider_dto_json
	};
	my ($variant_type, $source_general_consequence) = get_variant_type($obj);
	$var_obj->{variant_type_curie} = $variant_type if defined $variant_type;
	$var_obj->{source_general_consequence} = $source_general_consequence if defined $source_general_consequence;
	$var_obj->{variant_status_name} = 'public' unless $is_obsolete;
	
	if ($obj->Method) {
	    my $collection = $obj->Method->name;
	    unless ($collection eq 'Allele'
		    || $collection eq 'CGH_allele'
		    || $collection eq 'Deletion_allele'
		    || $collection eq 'Deletion_and_insertion_allele'
		    || $collection eq 'Deletion_polymorphism'
		    || $collection eq 'Engineered_allele'
		    || $collection eq 'Insertion_allele'
		    || $collection eq 'Insertion_polymorphism'
		    || $collection eq 'Mos_insertion'
		    || $collection eq 'SNP'
		    || $collection eq 'Substitution_allele'
		    || $collection eq 'Transposon_insertion'
		){
		$collection =~ s/ /_/g;
		$collection .= '_project' if $collection eq 'Million_mutation';
		$allele_obj->{in_collection_name} = $collection;
	    }
	}

	my $mutation_types = get_mutation_types($obj);
	$allele_obj->{allele_mutation_type_dtos} = $mutation_types if scalar @$mutation_types > 0;

	my @synonym_dtos;
	if ($obj->Other_name) {
	    for my $other_name ($obj->Other_name) {
		next if $other_name->name =~ /^[^:]+:[pcg]\./;
		my $synonym_dto = {
		    display_text   => $other_name->name,
		    format_text    => $other_name->name,
		    name_type_name => 'unspecified'
		};
		my $synonym_evidence = get_evidence_curies($other_name);
		$synonym_dto->{evidence_curies} = $synonym_evidence if @$synonym_evidence;
		push @synonym_dtos, $synonym_dto;
	    }
	}
	$allele_obj->{allele_synonym_dtos} = \@synonym_dtos if @synonym_dtos;
	
	my @secondary_ids;
	if ($obj->Acquires_merge) {
	    foreach my $g ($obj->Acquires_merge) {
		my $sid_json = {
		    secondary_id => 'WB:' . $g->name,
		    internal => JSON::false,
		    obsolete => JSON::false
		};
		push @secondary_ids, $sid_json;
	    }
	}
	$allele_obj->{secondary_id_dtos} = \@secondary_ids if @secondary_ids;
    

	if ($obj->Reference) {
	    $allele_obj->{reference_curies} = get_papers($obj);
	}

	my @note_dtos;
	if ($obj->Flanking_sequences) {
	    my $left = $obj->get('Flanking_sequences', 1);
	    my $right = $obj->get('Flanking_sequences', 2);
	    my $fs_note_json = {
		free_text => $left . '|' . $right,
		note_type_name => 'flanking_sequences',
		internal => JSON::false,
		obsolete => JSON::false
	    };
	    push @note_dtos, $fs_note_json; 
	}
	if ($obj->Deletion_verification) {
	    my $dv_note_json = {
		free_text => $obj->Deletion_verification->name,
		note_type_name => 'indel_verification',
		internal => JSON::false,
		obsolete => JSON::false
	    };
	    my @deletion_verification_paper_evidence;
	    for my $evi ($obj->Deletion_verification->col) {
		next unless $evi->name eq 'Paper_evidence';
		for my $paper ($evi->col) {
		    my $paper_id = get_paper($paper);
		    push @deletion_verification_paper_evidence, $paper_id if defined $paper_id;
		}
	    }
	    $dv_note_json->{evidence_curies} = \@deletion_verification_paper_evidence if @deletion_verification_paper_evidence > 0;
	    push @note_dtos, $dv_note_json;
	}
	if ($obj->Reference_strain_digest) {
	    for my $rsd ($obj->Reference_strain_digest) {
		my $site = $rsd->right;
		if (defined $site) {
		    my $rflp_string = 'Reference strain digest - Site: ' . $site;
		    my $enzyme = $rsd->right->right;
		    if (defined $enzyme) {
			$rflp_string .= '; Enzyme: ' . $enzyme;
			my $band_size = $rsd->right->right->right;
			if (defined $band_size) {
			    $rflp_string .= '; Band size: ' . $band_size;
			}
		    }
		    
		    my $rflp_note_json = {
			free_text => $rflp_string,
			note_type_name => 'restriction_allele',
			internal => JSON::false,
			obsolete => JSON::false
		    };
		    push @note_dtos, $rflp_note_json;
		}
	    }
	}
	if ($obj->Polymorphic_strain_digest) {
	    for my $psd ($obj->Polymorphic_strain_digest) {
		my $site = $psd->right;
		if (defined $site) {
		    my $rflp_string = 'Polymorphic strain digest - Site: ' . $site;
		    my $enzyme = $psd->right->right;
		    if (defined $enzyme) {
			$rflp_string .= '; Enzyme: ' . $enzyme;
			my $band_size = $psd->right->right->right;
			if (defined $band_size) {
			    $rflp_string .= '; Band size: ' . $band_size;
			}
		    }
		    
		    my $rflp_note_json = {
			free_text => $rflp_string,
			note_type_name => 'restriction_allele',
			internal => JSON::false,
			obsolete => JSON::false
		    };
		    push @note_dtos, $rflp_note_json;
		}
	    }
	}
	if ($obj->Mapping_target) {
	    my $mt_note_json = {
		free_text => $obj->Mapping_target->name,
		note_type_name => 'positive_clone',
		internal => JSON::false,
		obsolete => JSON::false
	    };
	    push @note_dtos, $mt_note_json;
	}

	if ($obj->Detection_method) {
	    my $dm_note_json = {
		free_text => $obj->Detection_method->name,
		note_type_name => 'detection_method',
		internal => JSON::false,
		obsolete => JSON::false
	    };
	    push @note_dtos, $dm_note_json;
	}
	if ($obj->Remark) {
	    for my $remark ($obj->Remark) {
		my $remark_str = $remark->name;
		my $comment_note_json = {
		    free_text => "$remark_str",
		    note_type_name => 'comment',
		    internal => JSON::false,
		    obsolete => JSON::false
		};
		my @remark_paper_evidence;
		for my $evi ($remark->col) {
		    next unless $evi->name eq 'Paper_evidence';
		    for my $paper ($evi->col) {
			my $paper_id = get_paper($paper);
			push @remark_paper_evidence, $paper_id if defined $paper_id;
		    }
		}
		$comment_note_json->{evidence_curies} = \@remark_paper_evidence if @remark_paper_evidence > 0;	
		push @note_dtos, $comment_note_json;
	    }
	}

	# TODO: work out which notes belong on variants and which on alleles
	$allele_obj->{note_dtos} = \@note_dtos if @note_dtos;
	$var_obj->{note_dto} = \@note_dtos if @note_dtos;

	# Stick some random data in for curation DB test files
	if ($curation_test) {
	    my @inheritance_modes = qw(dominant recessive semi-dominant unknown);
	    $allele_obj->{is_extinct} = $var_count % 2 == 0 ? JSON::true : JSON::false;
	    $allele_obj->{date_updated} = get_random_datetime(10);
	    $allele_obj->{date_created} = get_random_datetime(1);
	}

	my $allele_json = JSON->new;
	my $allele_string = $allele_json->allow_nonref->convert_blessed->canonical->pretty->encode($allele_obj);
	$allele_out_fh->print(',') unless $var_count == 1;
	$allele_out_fh->print("\n" . '      ' . $allele_string);
	
	my $variant_json = JSON->new;
	my $variant_string = $variant_json->allow_nonref->convert_blessed->canonical->pretty->encode($var_obj);
	$variant_out_fh->print(',') unless $var_count == 1;
	$variant_out_fh->print("\n" . '      ' . $variant_string);

	last if $limit && $limit <= ($var_count * 2);
    }

    return;
}

sub process_transgenes {
    my @transgenes = $db->fetch(-query => 'find Transgene WHERE Species = "Caenorhabditis elegans"');

    my $var_count = 0;
    for my $obj(@transgenes) {
	$var_count++;
	next unless $obj->isObject();
	unless ($obj->Species) {
	    print "No species for $obj - skipping\n";
	    next;
	}
	
	my $dp_xref_dto_json = {
	    referenced_curie => 'WB:' . $obj->name,
	    page_area => 'transgene',
	    display_name => $obj->name,
	    prefix => 'WB',
	    internal => JSON::false,
	    obsolete => JSON::false
	};
	
	my $data_provider_dto_json = {
	    source_organization_abbreviation => 'WB',
	    cross_reference_dto => $dp_xref_dto_json,
	    internal => JSON::false,
	    obsolete => JSON::false
	}; 
	
	my $allele_obj = {
	    curie         => "WB:" . $obj->name, 
	    allele_symbol_dto => get_symbol_dto($obj),
	    taxon_curie   => "NCBITaxon:" . $obj->Species->NCBITaxonomyID->name,
	    internal      => JSON::false,
	    obsolete      => JSON::false,
	    data_provider_dto => $data_provider_dto_json
	};	

	my @synonym_dtos;
	if ($obj->Synonym) {
	    for my $synonym ($obj->Synonym) {
		my $synonym_dto = {
		    display_text   => $synonym->name,
		    format_text    => $synonym->name,
		    name_type_name => 'unspecified',
		    internal      => JSON::false
		};
		push @synonym_dtos, $synonym_dto;
	    }
	}
	$allele_obj->{allele_synonym_dtos} = \@synonym_dtos if @synonym_dtos;

	if ($obj->Reference) {
	    $allele_obj->{reference_curies} = get_papers($obj);
	}
	
	my @note_dtos;

	if ($obj->Summary) {
	    my $cs_note_json = {
		free_text => $obj->Summary->name,
		note_type_name => 'transgene_content_summary',
		internal => JSON::false,
		obsolete => JSON::false
	    };
		my @cs_paper_evidence;
	    for my $evi ($obj->Summary->col) {
		next unless $evi->name eq 'Paper_evidence';
		for my $paper ($evi->col) {
		    my $paper_id = get_paper($paper);
		    push @cs_paper_evidence, $paper_id if defined $paper_id;
		}
	    }
	    $cs_note_json->{evidence_curies} = \@cs_paper_evidence if @cs_paper_evidence > 0;
	    push @note_dtos, $cs_note_json;
	}

	if ($obj->Construction_summary) {
	    my $cs_note_json = {
		free_text => $obj->Construction_summary->name,
		note_type_name => 'transgene_construction_summary',
		internal => JSON::false,
		obsolete => JSON::false
	    };
	    push @note_dtos, $cs_note_json;
	}

	if ($obj->Coinjection_other) {
	    my $co_note_json = {
		free_text => $obj->Coinjection_other->name,
		note_type_name => 'coinjection_other',
		internal => JSON::false,
		obsolete => JSON::false
	    };
	    push @note_dtos, $co_note_json;
	}
	
	if ($obj->Remark) {
	    for my $remark ($obj->Remark) {
		my $remark_str = $remark->name;
		print $remark_str . "\n";
		my $comment_note_json = {
		    free_text => "$remark_str",
		    note_type_name => 'comment',
		    internal => JSON::false,
		    obsolete => JSON::false
		};
		my @remark_paper_evidence;
		for my $evi ($remark->col) {
		    next unless $evi->name eq 'Paper_evidence';
		    for my $paper ($evi->col) {
			my $paper_id = get_paper($paper);
			push @remark_paper_evidence, $paper_id if defined $paper_id;
		    }
		}
		$comment_note_json->{evidence_curies} = \@remark_paper_evidence if @remark_paper_evidence > 0;	
		push @note_dtos, $comment_note_json;
	    }
	}
      	
	$allele_obj->{note_dtos} = \@note_dtos if @note_dtos;
	
	# Stick some random data in for curation DB test files
	if ($curation_test) {
	    $allele_obj->{is_extinct} = $var_count % 2 == 0 ? JSON::true : JSON::false;
	    $allele_obj->{date_updated} = get_random_datetime(10);
	    $allele_obj->{date_created} = get_random_datetime(1);
	}
	
	my $json = JSON->new;
	my $string = $json->allow_nonref->canonical->pretty->encode($allele_obj);
	$allele_out_fh->print(",\n" . '      ' . $string);
	
	last if $limit && $limit <= ($var_count * 2);
    }

    return;
}

sub get_variant_type {
    my $obj = shift;

    my ($variant_type, $source_general_consequence);
    if ($obj->Engineered_allele) {
	$variant_type = $TERM2SO{'Engineered_allele'};
    } elsif ($obj->Allele) {
	$variant_type = $TERM2SO{'Allele'};
    } elsif ($obj->SNP || $obj->Confirmed_SNP || $obj->Predicted_SNP) {
	$variant_type = $TERM2SO{'SNP'};
    } elsif ($obj->RFLP) {
	$variant_type = $TERM2SO{'RFLP'};
    } elsif ($obj->Transposon_insertion) {
	$variant_type = $TERM2SO{'Transposon_insertion'};
    } elsif ($obj->Natural_variant) {
	$variant_type = $TERM2SO{'Natural_variant'};
    }

    if ($obj->Substitution) {
	my @subs = $obj->Substitution->row;
	if (scalar @subs >= 2 && length($subs[0]) == 1 && length($subs[1]) == 1) {
	    if ($obj->Natural_variant) {
		$source_general_consequence = $TERM2SO{'SNP'};
	    } else {
		$source_general_consequence = $TERM2SO{'Point_mutation'};
	    }
	} else {
	    $source_general_consequence = $TERM2SO{'Delins'};
	}
    } elsif ($obj->Insertion) {
	$source_general_consequence = $TERM2SO{'Insertion'};
    } elsif ($obj->Deletion) {
	$source_general_consequence = $TERM2SO{'Deletion'};
    } elsif ($obj->Tandem_duplication) {
	$source_general_consequence = $TERM2SO{'Tandem_duplication'};
    }

    if (!defined $variant_type) {
	$variant_type = $source_general_consequence;
	$source_general_consequence = undef();
    }

    return ($variant_type, $source_general_consequence);
}


sub get_papers {
    my $obj = shift;

    my @references;
    for my $ref ($obj->Reference) {
	next unless defined $ref;
	my $publication_id = get_paper($ref);
	next unless defined $publication_id;
	push @references, $publication_id;
    }

    return \@references;
}

sub get_paper {
    my $ref = shift;

    if (!$ref->Status) {
	return;
    }

    my $level = 1;
    while ($ref->Merged_into && $level < 6) {
	$level++;
	$ref = $ref->Merged_into;
    }
    return if $ref->Status->name eq 'Invalid';

    my $pmid;
    for my $db ($ref->Database) {
	if ($db->name eq 'MEDLINE') {
	    $pmid = $db->right->right->name;
	    last;
	}
    }
    
    my $publication_id = $pmid ? "PMID:$pmid" : 'WB:' . $ref->name;
    if ($publication_id eq 'WB:WBPaper000045183') {
	$publication_id = 'WB:WBPaper00045183';
    }
    if ($publication_id eq 'WB:WBPaper000042571') {
	$publication_id = 'WB:WBPaper00042571';
    }

    return $publication_id;
}

sub get_random_datetime {
    my $add = shift;
    
    my $year = int(rand(10)) + $add;
    my $month = int(rand(11)) + 1;
    my $day = int(rand(27)) + 1;

    my $year_string = $year > 9 ? $year : '0' . $year;
    my $month_string = $month > 9 ? $month : '0' . $month;
    my $day_string = $day > 9 ? $day : '0' . $day;
    
    return '20' . $year_string . '-' . $month_string . '-' . $day_string . 'T00:00:00+00:00';
}

sub get_mutation_types {
    my $obj = shift;

    my @mutation_types;
    if ($obj->Substitution) {
	my $substitution = {
	   mutation_type_curies => [$TERM2SO{'Substitution'}],
	   internal => JSON::false
	};
	if ($obj->Substitution->right) {
	    my @substitution_paper_evidence;
	    for my $evi ($obj->Substitution->right->col) {
		next unless $evi->name eq 'Paper_evidence';
		for my $paper ($evi->col) {
		    my $paper_id = get_paper($paper);
		    push @substitution_paper_evidence, $paper_id if defined $paper_id;
		}
	    }
	    $substitution->{evidence_curies} = \@substitution_paper_evidence if @substitution_paper_evidence > 0;
	}
	push @mutation_types, $substitution;
    }

    if ($obj->Insertion) {
	my $insertion = {
	    mutation_type_curies => [$TERM2SO{'Insertion'}],
	    internal => JSON::false
	};
	if ($obj->Insertion->right) {
	    my @insertion_paper_evidence;
	    for my $evi ($obj->Insertion->col) {
		next unless $evi->name eq 'Paper_evidence';
		for my $paper ($evi->col) {
		    my $paper_id = get_paper($paper);
		    push @insertion_paper_evidence, $paper_id if defined $paper_id;
		}
	    }
	    $insertion->{evidence_curies} = \@insertion_paper_evidence if @insertion_paper_evidence > 0;
	}
	push @mutation_types, $insertion;
    }

    if ($obj->Deletion) {
	my $deletion = {
	    mutation_type_curies => [$TERM2SO{'Deletion'}],
	    internal => JSON::false
	};
	if ($obj->Deletion->right) {
	    my @deletion_paper_evidence;
	    for my $evi ($obj->Deletion->col) {
		next unless $evi->name eq 'Paper_evidence';
		for my $paper ($evi->col) {
		    my $paper_id = get_paper($paper);
		    push @deletion_paper_evidence, $paper_id if defined $paper_id;
		}
	    }
	    $deletion->{evidence_curies} = \@deletion_paper_evidence if @deletion_paper_evidence > 0;
	}
	push @mutation_types, $deletion;
    }

    if ($obj->Tandem_duplication) {
	my $duplication = {
	    mutation_type_curies => [$TERM2SO{'Tandem_duplication'}],
	    internal => JSON::false
	};
	if ($obj->Tandem_duplication->right) {
	    my @duplication_paper_evidence;
	    for my $evi ($obj->Tandem_duplication->col) {
		next unless $evi->name eq 'Paper_evidence';
		for my $paper ($evi->col) {
		    my $paper_id = get_paper($paper);
		    push @duplication_paper_evidence, $paper_id if defined $paper_id;
		}
	    }
	    $duplication->{evidence_curies} = \@duplication_paper_evidence if @duplication_paper_evidence > 0;
	}
	push @mutation_types, $duplication
    };

    return \@mutation_types;
}

sub get_symbol_dto {
    my $obj = shift;

    my $symbol_dto;
    if ($obj->Public_name) {
	$symbol_dto = {
	    display_text       => $obj->Public_name->name,
	    format_text        => $obj->Public_name->name,
	    name_type_name     => 'nomenclature_symbol',
	    internal      => JSON::false
	};
	if ($obj->Public_name->right) {
	    my $symbol_evidence = get_evidence_curies($obj->Public_name->right);
	    $symbol_dto->{evidence_curies} = $symbol_evidence if @$symbol_evidence;
	}
    } else {
	$symbol_dto = {
	    display_text       => $obj->name,
	    format_text        => $obj->name,
	    name_type_name     => 'systematic_name',
	    internal      => JSON::false
	};
    }

    return $symbol_dto;
}

sub get_evidence_curies {
    my $name_obj = shift;

    my @paper_evidence;
    for my $evi ($name_obj->col) {
	next unless $evi->name eq 'Paper_evidence';
	for my $paper ($evi->col) {
	    my $paper_id = get_paper($paper);
	    push @paper_evidence, $paper_id if defined $paper_id;
	}
    }

    return \@paper_evidence;
}

    
