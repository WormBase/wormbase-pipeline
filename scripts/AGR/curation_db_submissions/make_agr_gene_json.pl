#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;
use Const::Fast;
use Wormbase;

my ($help, $debug, $test, $verbose, $store, $wormbase, $schema);
my ($outfile, $acedbpath, $ws_version, $out_fh, $gtf_file);

const my $LINKML_SCHEMA => 'v2.9.1';

const my %XREF_MAP => (
    "NCBI"       => {
	gene => "NCBI_Gene",
    },
    "SwissProt"  => {
	UniProtAcc => "UniProtKB",
    },
    "TrEMBL"     => {
	UniProtAcc => "UniProtKB",
    },
    "RNAcentral" => {
	URSid => "RNAcentral",
    },
    "Panther" => {
	family => "PANTHER",
    },
    "RefSeq" => {
	protein => "RefSeq",
    },
    "UniProt_GCRP" => {
	UniProtAcc => "UniProtKB",
    },
);



GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "test"        => \$test,
	    "verbose"     => \$verbose,
	    "database:s"  => \$acedbpath,
	    "outfile:s"   => \$outfile,
            "wsversion=s" => \$ws_version,
	    "gtf=s"       => \$gtf_file
	    );

$wormbase = Wormbase->new(-debug => $debug,-test => $test);

my $tace = $wormbase->tace;

$acedbpath = $wormbase->autoace unless $acedbpath;

print "Connecting to AceDB\n" if $verbose;
my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or die("Connection failure: ". Ace->error);

my $query = 'FIND Gene WHERE Species = "Caenorhabditis elegans"';

$outfile = "./wormbase.genes.${ws_version}.${LINKML_SCHEMA}.json" unless defined $outfile;

my @genes;

print "Getting location data\n" if $verbose;
my $locs = get_location_data($db, $gtf_file);


my %unmapped_xrefs;
print "Querying AceDB\n" if $verbose;
my $it = $db->fetch_many(-query => $query);
while (my $obj = $it->next) {
    next unless $obj->isObject();
    next unless $obj->name =~ /^WBGene/;
    unless ($obj->Species) {
	print "No species for $obj - skipping\n";
	next;
    }
    print "Processing $obj\n" if $verbose;
    my ($symbol, $full_name, $systematic_name, $synonyms) = get_name_slot_annotations($obj);
   
    if (!defined $symbol && $obj->Status && $obj->Status->name eq 'Live') {
	print "Live gene with no symbol $obj - skipping\n";
	next;
    }

    # Secondary ids
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
    
    my $dp_xref_dto_json = {
	referenced_curie => 'WB:' . $obj->name,
	page_area => 'gene',
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

    my @xrefs;
    my $gcrp_xref;
    # Cross references
    foreach my $dblink ($obj->Database) {
	if (exists $XREF_MAP{$dblink}) {
	    foreach my $field ($dblink->col) {
		if (exists $XREF_MAP{$dblink}->{$field}) {
		    my $prefix = $XREF_MAP{$dblink}->{$field};
		    if ($dblink eq 'UniProt_GCRP') {
			my $suffix = $field->right->name;
			$gcrp_xref = {
			    referenced_curie => "$prefix:$suffix",
			    page_area => "default",
			    display_name => "$prefix:$suffix",
			    prefix => $prefix
			};
		    } else {
			my @ids = $field->col;
			foreach my $id(@ids){
			    my $suffix = $id->name; 
			    push @xrefs, {
				referenced_curie => "$prefix:$suffix",
				page_area => "default",
				display_name => "$prefix:$suffix",
				prefix => $prefix
			    };
			}
		    }
		}
	    }
	} else {
	    foreach my $field ($dblink->col) {
		$unmapped_xrefs{$dblink}{$field}++;
	    }
	}
    }
    if (defined $gcrp_xref && @xrefs) {
	my @filtered_xrefs;
	for my $xref (@xrefs) {
	    if ($xref->{referenced_curie} ne $gcrp_xref->{referenced_curie}) {
		push @filtered_xrefs, $xref;
	    }
	}
	@xrefs = @filtered_xrefs;
    }

    # Back-to-mod xrefs
    my @baseXrefs = qw(gene gene/expression gene/references);
    push @baseXrefs, 'gene/phenotypes' if ($obj->Allele and $obj->Allele->Phenotype) || ($obj->RNAi_result and $obj->RNAi_result->Phenotype);
    push @baseXrefs, 'gene/expression_images' if ($obj->Expr_pattern && grep {$_->Picture} $obj->Expr_pattern);

    for my $baseXref (@baseXrefs) {
	push @xrefs, {
	    referenced_curie => 'WB:' . $obj->name,
	    page_area => $baseXref,
	    display_name => 'WB:' . $obj->name,
	    prefix => "WB"
	};
    }
    # for SPELL, we should technically only include the ref if this gene exists in SPELL, 
    # which (as defined by Wen Chen) is: 
    if ($obj->RNASeq_FPKM or $obj->Microarray_results) {
	push @xrefs, {
	    referenced_curie => "WB:" . $symbol->{display_text},
	    page_area => "gene/spell",
	    display_name => $symbol->{display_text},
	    prefix => "WB"
	};
    }

    # Ensembl xrefs, only for genomically placed genes
    if (defined $locs and exists $locs->{$obj->name}) {
	push @xrefs, {
	    referenced_curie => 'ENSEMBL:' . $obj->name,
	    page_area => 'default',
	    display_name => 'ENSEMBL:' . $obj->name,
	    prefix => 'ENSEMBL'
	};
    }
    


    my $is_obsolete = $obj->Status && $obj->Status->name eq 'Dead' ? JSON::true : JSON::false;
    my $gene_type = $obj->Biotype ? $obj->Biotype->name : 'SO:0000704';
    my $gene = {
	primary_external_id => 'WB:' . $obj->name,
	gene_type_curie  => $gene_type,
	gene_symbol_dto   => $symbol,
	taxon_curie    => 'NCBITaxon:' . $obj->Species->NCBITaxonomyID,
	obsolete => $is_obsolete,
	internal => JSON::false,
	created_by_curie => 'WB:curator',
	updated_by_curie => 'WB:curator',
	data_provider_dto => $data_provider_dto_json
    };


    $gene->{gene_full_name_dto} = $full_name if $full_name;
    $gene->{gene_systematic_name_dto} = $systematic_name if $systematic_name;
    $gene->{gene_synonym_dtos} = $synonyms if @$synonyms;
    $gene->{gene_secondary_id_dtos} = \@secondary_ids if @secondary_ids;
    $gene->{cross_reference_dtos} = \@xrefs if @xrefs;
    $gene->{gcrp_cross_reference_dto} = $gcrp_xref if $gcrp_xref;
    
    push @genes, $gene;
}

my $data = {
    alliance_member_release_version => $ws_version,
    linkml_version => $LINKML_SCHEMA,
    gene_ingest_set => \@genes,
};

open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";

my $json_obj = JSON->new;
my $string = $json_obj->allow_nonref->canonical->pretty->encode($data);
print $out_fh $string;

$db->close;

open XREF, ">unmapped_gene_xrefs.txt" or die "Could not open unmapped_xrefs.txt for writing\n";
for my $xref_db (keys %unmapped_xrefs) {
    for my $xref_type (keys %{$unmapped_xrefs{$xref_db}}) {
	print XREF $xref_db . ' - ' . $xref_type . ' (' . $unmapped_xrefs{$xref_db}{$xref_type} . ')' . "\n";
    }
}
exit(0);

sub get_name_slot_annotations {
    my $obj = shift;

    my ($symbol_obj, $systematic_name_obj);
    my $symbol_type = "nomenclature_symbol";
    my @synonym_objs = $obj->Other_name;
    if ($obj->CGC_name) {
	$symbol_obj = $obj->CGC_name;
    }
    if ($obj->Sequence_name) {
	unless (defined $symbol_obj) {
	    $symbol_obj = $obj->Sequence_name;
	    $symbol_type = "systematic_name";
	}
	$systematic_name_obj = $obj->Sequence_name;
    }
    if ($obj->Public_name) {
	if (defined $symbol_obj) {
	    push @synonym_objs, $obj->Public_name unless $symbol_obj->name eq $obj->Public_name->name;
	} else {
	    $symbol_obj = $obj->Public_name;
	}
    }

    my $symbol;
    if ($symbol_obj) {
	my $symbol_evidence = get_evidence_curies($symbol_obj);
	$symbol = {
	    display_text => $symbol_obj->name,
	    format_text  => $symbol_obj->name,
	    name_type_name => $symbol_type,
	    internal => JSON::false
	};
	$symbol->{evidence_curies} = $symbol_evidence if @$symbol_evidence;
    } else {
	$symbol = {
	    display_text   => $obj->name,
	    format_text    => $obj->name,
	    name_type_name => 'systematic_name',
	    internal => JSON::false
	};
    }

    my $systematic_name;
    if ($systematic_name_obj) {
	my $systematic_name_evidence = get_evidence_curies($systematic_name_obj);
	$systematic_name = {
	    display_text => $systematic_name_obj->name,
	    format_text => $systematic_name_obj->name,
	    name_type_name => "systematic_name",
	    internal => JSON::false
	};
	$systematic_name->{evidence_curies} = $systematic_name_evidence if @$systematic_name_evidence;
    }

    my @synonyms;
    if (@synonym_objs) {
	for my $synonym_obj (@synonym_objs) {
	    my $symbol_evidence = get_evidence_curies($synonym_obj);
	    my $synonym = {
		format_text => $synonym_obj->name,
		display_text => $synonym_obj->name,
		name_type_name => "unspecified",
		internal => JSON::false
	    };
	    $synonym->{evidence_curies} = $symbol_evidence if @$symbol_evidence;
	    push @synonyms, $synonym;
	}
    }

    my ($full_name, $full_name_evidence);
        
    if ($obj->Gene_class and $obj->Gene_class->Description) {
	my $suffix = '';
	if ($symbol_obj) {
	    if ($symbol_obj->name =~ /-(\S+)$/) {
		$suffix = ' ' . $1;
	    }
	}
	$full_name_evidence = get_evidence_curies($obj->Gene_class->Description);
	$full_name = {
	    display_text => $obj->Gene_class->Description->name . $suffix,
	    format_text => $obj->Gene_class->Description->name . $suffix,
	    name_type_name => "full_name",
	    internal => JSON::false
	};
	$full_name->{evidence_curies} = $full_name_evidence if @$full_name_evidence;
    } elsif ($obj->Corresponding_CDS && $obj->Corresponding_CDS->Brief_identification) {
	$full_name_evidence = get_evidence_curies($obj->Corresponding_CDS->Brief_identification);
	$full_name = {
	    display_text => $obj->Corresponding_CDS->Brief_identification->name,
	    format_text => $obj->Corresponding_CDS->Brief_identification->name,
	    name_type_name => "full_name",
	    internal => JSON::false
	};
	$full_name->{evidence_curies} = $full_name_evidence if @$full_name_evidence;
    }
	
    return ($symbol, $full_name, $systematic_name, \@synonyms);
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



sub get_location_data {
  my ($acedb, $gtf) = @_;

  # get the assembly name for the canonical bioproject
  my ($assembly_name, $fh, %locs);
  
  my $species_obj = $acedb->fetch(-class => 'Species', -name => $wormbase->full_name);
  my @seq_col = $species_obj->at('Assembly');
        
  foreach my $seq_col_name (@seq_col) {
    my ($bioproj);

    my $seq_col = $seq_col_name->fetch;
    my $this_assembly_name = $seq_col->Name;

    my @db = $seq_col->at('Origin.DB_info.Database');
    foreach my $db (@db) {
      if ($db->name eq 'NCBI_BioProject') {
        $bioproj = $db->right->right->name;
      }
    }
    
    if (defined $bioproj and $wormbase->ncbi_bioproject eq $bioproj) {
      $assembly_name = $this_assembly_name->name;
      last;
    }
  }
    
  if (not defined $assembly_name) {
    die "Could not find name of current assembly for " . $wormbase->species() . "\n";
  }

  # parse the GTF
  if ($gtf =~ /\.gz$/) {
    open( $fh, "gunzip -c $gtf |") 
        or die "Could not open gunzip stream GTF file $gtf\n";
  } else {
    open( $fh, $gtf) 
        or die "Could not open GTF file $gtf\n";
  }

  while(<$fh>) {
    /^\#/ and next;
    
    my @l = split(/\t/, $_);
    next if $l[2] ne "gene";

    my ($wbgid) = $l[8] =~ /gene_id\s+\"(\S+)\"/; 

    $locs{$wbgid} = {
      assembly       => $assembly_name,
      chromosome     => $l[0],
      startPosition  => $l[3] + 0, 
      endPosition    => $l[4] + 0,
      strand         => $l[6],
    };
  }
      
  return \%locs;
}

1;
