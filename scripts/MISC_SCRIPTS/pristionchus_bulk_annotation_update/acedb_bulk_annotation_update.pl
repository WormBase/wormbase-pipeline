#!/usr/bin/perl
use strict;
use warnings;

use lib $ENV{'CVS_DIR'};
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";

use Getopt::Long;
use Path::Class;
use Wormbase;
use NameDB_handler;
use Log_files;
use Ace;
use Const::Fast;
use Curate;

const my @ISOFORM_SUFFIXES => qw(a b c d e f g h i j k l m n o p q r s t u v w x y z);
    
my ($create_file, $delete_file, $split_file, $merge_file, $update_file, $wb_gff_file, $gff_file, $geneace_out_file, $curation_out_file, $test, $log_file, $species, $debug, $genedb, $curationdb, $nameserver, $person, $verbose, $start_codons_file, $processed_file);
my @references;

GetOptions(
    "create:s"     => \$create_file,
    "delete:s"     => \$delete_file,
    "split:s"      => \$split_file,
    "merge:s"      => \$merge_file,
    "update:s"     => \$update_file,
    "gff=s"        => \$gff_file,
    "wbgff=s"      => \$wb_gff_file,
    "curation=s"   => \$curation_out_file,
    "geneace=s"    => \$geneace_out_file,
    "genedb=s"     => \$genedb,
    "curationdb=s" => \$curationdb,
    "references=s" => \@references,
    "species:s"    => \$species,
    "test"         => \$test,
    "verbose"      => \$verbose,
    "log=s"        => \$log_file,
    "person=s"     => \$person,
    "processed=s"  => \$processed_file,
    "debug:s"      => \$debug,
    "nameserver"   => \$nameserver
    ) or die ($@);

@references = split(/,/, join(',', @references));

$curation_out_file = './bulk_update.curation.ace' unless defined $curation_out_file;
$geneace_out_file = './bulk_update.geneace.ace' unless defined $geneace_out_file;
$processed_file = './genes_processed_in_names_service.txt' unless defined $processed_file;
if (!-e $processed_file) {
    system ("touch $processed_file");
}
my $curation_out_fh = file($curation_out_file)->openw;
my $geneace_out_fh = file($geneace_out_file)->openw;
$genedb = '/nfs/production/flicek/wormbase/wb/CURATION_DATABASES/geneace_mqt' unless $genedb;
$curationdb = '/nfs/production/flicek/wormbase/wb/CURATION_DATABASES/pristionchus_mqt' unless $curationdb;

my $log = Log_files->make_log($log_file, $debug);
my $wb = Wormbase->new(-organism => $species, -debug => $debug, -test => $test);
my $gdb = Ace->connect(-path => $genedb);
my $cdb = Ace->connect(-path => $curationdb);
my $nh = NameDB_handler->new($wb, $test);
my $curation = Curate->new($geneace_out_fh, $cdb, $wb, $log);

if ($species eq 'elegans' || $species eq 'sratti' || $species eq 'tmuris') {
    $log->log_and_die("This script is not intended for $species updates - see comments in code\n");
}

my $test_count = 0;

my ($day, $mon, $yr)  = (localtime)[3,4,5];
my $date = sprintf("%02d/%02d/%04d", $day, $mon + 1, $yr + 1900);

my @changing_wb_genes;
my %wb_cds_transcript_map;
my %history_created;

my $wormpep_prefix = lc($wb->wormpep_prefix);
my $new_gene_models = parse_gff($gff_file, 0);
my $wb_gene_models = parse_gff($wb_gff_file, 1);

my $processed_genes = parse_single_column_file($processed_file);
my %ns_processed = map {$_ => 1} @$processed_genes;
my $nsp_out_fh = file($processed_file)->opena;

my $to_update = parse_mapping_file($update_file, 0);
my $to_merge = parse_mapping_file($merge_file, 1);
my $to_split = parse_mapping_file($split_file, 1);

check_for_merge_split_clashes();

my %cgc_names;
my %wb_ids_created;

$log->write_to("Creating genes\n");
create_genes() if defined $create_file;

$log->write_to("Deleting genes\n");
delete_genes() if defined $delete_file;

$log->write_to("Splitting genes\n");
split_genes() if scalar keys %$to_split > 0;

$log->write_to("Merging genes\n");
merge_genes() if scalar keys %$to_merge > 0;

$log->write_to("Updating genes\n");
update_genes() if scalar keys %$to_update > 0;

$log->end;
$log->mail;
exit(0);

sub check_for_merge_split_clashes {
    my $conflicts = 0;
    for my $merge_external (keys %$to_merge) {
	for my $wb_merge (@{$to_merge->{$merge_external}}) {
	    if (exists $to_split->{$wb_merge}) {
		$log->error("$wb_merge involved in both merge and split operations - this script cannot currently resolve this situation\n");
		$conflicts++;
	    }
	    if (exists $to_update->{$wb_merge}) {
		$log->error("$wb_merge involved in bother merge and update operations - this script cannot currently resolve this situation\n");
		$conflicts++;
	    }
	}
    }
    if ($conflicts) {
	$log->log_and_die("Cannot continue as $conflicts split/merge/update conflicts\n");
    }
}

sub create_genes {		      
    my $ids_to_create = parse_single_column_file($create_file);
    for my $id (@$ids_to_create) {
	create_gene($id);
    }
    
    return;
}

sub create_gene {
    my ($non_wb_id, $split_from) = @_;
    
    # If this code is ever to be used for elegans, sratti, or tmuris updates then we would need to supply a clone ID to Next_CDS_ID here
    my $seq_name = $curation->Next_CDS_ID();
    my $so_term = 'SO:0001217'; # Are we really only creating coding genes or could we be creating Pseudogenes?
    
    my $wb_id = $nh->idGetByTypedName('Sequence' => $seq_name);
    if ($wb_id) {
	$log->write_to("$wb_id already exists in Names Service, creating using existing ID\n");
    } else {
	if ($split_from) {
	    $log->write_to("Splitting $split_from and creating $seq_name\n") if $verbose;
	    $wb_id = $nh->idSplit($split_from, 'cds', $seq_name, 'cds');
	    $log->write_to("NS->idSplit $wb_id from $split_from\n") if $verbose;
	    if (!defined $wb_id) {
		$log->error("Could not split gene $split_from in Names Services\n");
		return;
	    }
	} else {
	    my $new_gene = $nh->new_gene($seq_name, 'Sequence', $wb->full_name);
	    if (defined $new_gene) {
		$wb_id = $new_gene->{'id'};
	    } else {
		$log->error("Could not retrieve new gene ID from Names Service for $seq_name\n");
		return;
	    }
	}
    }

    $log->write_to("Creating new Gene $wb_id for $seq_name for $non_wb_id\n") if $verbose;
    
    my $event = defined $split_from ? 'Split_from ' . $split_from : 'Created';
    
    # This code comes from GENEACE/newgene.pl - additional code required if this script is ever to be used for elegans
    $geneace_out_fh->print("Gene : $wb_id\n");
    $geneace_out_fh->print("Live\n");
    $geneace_out_fh->print("Version 1\n");
    $geneace_out_fh->print("Biotype \"$so_term\"\n");
    $geneace_out_fh->print("Sequence_name \"$seq_name\"\n");
    $geneace_out_fh->print("Species \"${\$wb->full_name}\"\n");
    $geneace_out_fh->print("History Version_change 1 now $person Event $event\n");
    $geneace_out_fh->print("Split_from $split_from\n") if defined $split_from;
    $geneace_out_fh->print("Method Gene\n");
    $geneace_out_fh->print("Public_name \"$seq_name\"\n\n");
    
    my $ix = 0;
    for my $external_transcript_id(keys %{$new_gene_models->{$non_wb_id}{'transcripts'}}) {
	$ix++;
	my $wb_cds_name = $seq_name;
	if (scalar keys %{$new_gene_models->{$non_wb_id}{'transcripts'}} > 1) {
	    $wb_cds_name .= $ISOFORM_SUFFIXES[$ix];
	}
	create_cds($non_wb_id, $external_transcript_id, $wb_cds_name, $wb_id);
    }
    $wb_ids_created{$non_wb_id} = $wb_id;
    return $wb_id;
}


sub create_cds {
    my ($external_gene_id, $external_transcript_id, $wb_cds_name, $wb_gene_id) = @_;
    
    my ($cds_starts, $cds_ends) = get_relative_positions($external_gene_id, $external_transcript_id, 'CDS');
    
    my @starts_on_sequence = @{$new_gene_models->{$external_gene_id}{'transcripts'}{$external_transcript_id}{'CDS'}{'starts'}};
    my @ends_on_sequence = @{$new_gene_models->{$external_gene_id}{'transcripts'}{$external_transcript_id}{'CDS'}{'ends'}};
    my $remark = 'Created as part of bulk annotation update on ' . $date;
    my $sequence = $new_gene_models->{$external_gene_id}{'gene'}{'chromosome'};
    my $method = 'curated';
    
    $curation_out_fh->print("CDS : \"$wb_cds_name\"\n");
    $curation_out_fh->print("Gene $wb_gene_id\n");
    $curation_out_fh->print("CDS\n");
    $curation_out_fh->print("Sequence $sequence\n");
    $curation_out_fh->print("Species \"${\$wb->full_name}\"\n");
    $curation_out_fh->print("Method $method\n");
    if (scalar @references > 0) {
	for my $ref (@references) {
	    $curation_out_fh->print("Remark \"$remark\" Paper_evidence $ref\n");
	}
    } else {
	$curation_out_fh->print("Remark \"$remark\"\n");
    }
    for my $ix(0 .. scalar @$cds_starts - 1) {
	$curation_out_fh->print("Source_exons $cds_starts->[$ix] $cds_ends->[$ix]\n");
    }

    $curation_out_fh->print("\nGene : $wb_gene_id\n");
    $curation_out_fh->print("Corresponding_CDS $wb_cds_name\n");
    
    $curation_out_fh->print("\nSequence : \"$sequence\"\n" );
    my $start_on_sequence = $new_gene_models->{$external_gene_id}{'gene'}{'strand'} eq '+' ? $starts_on_sequence[0] : $ends_on_sequence[-1];
    my $end_on_sequence = $new_gene_models->{$external_gene_id}{'gene'}{'strand'} eq '+' ? $ends_on_sequence[-1] : $starts_on_sequence[0];
    $curation_out_fh->print("CDS_child \"$wb_cds_name\" $start_on_sequence $end_on_sequence\n\n");
    
}

sub delete_genes {
    my $wb_genes_to_delete = parse_single_column_file($delete_file);
    for my $wb_id (@{$wb_genes_to_delete}) {
	delete_gene($wb_id);
    }
    
    return;
}

sub delete_gene {
    my ($wb_id, $merged_into) = @_;
    
    my $remark = 'Removed as part of bulk annotation update on ' . $date;
    
    my $camace_gene = $cdb->fetch(Gene => $wb_id);
    my $geneace_gene = $gdb->fetch(Gene => $wb_id);
    if ($camace_gene && $geneace_gene) {

	if ($camace_gene->Corresponding_CDS) {
	    for my $cds ($camace_gene->Corresponding_CDS) {
		make_history_cds($cds->name, 'CDS deprecated as part of bulk annotation update on ' . $date);
	    }
	}
	
	my $version = $geneace_gene->Version->name;
	$log->error("Lower gene version for $wb_id in geneace than in curation DB\n") if $camace_gene->Version && $camace_gene->Version->name > $version;

	if ($geneace_gene->CGC_name) {
	    remove_cgc_name($geneace_gene->name);
	    $version++;
	}


	$log->write_to("Deleting gene $wb_id\n") if $verbose;
	if ($nameserver && !exists $ns_processed{$wb_id}) {
	    if ($merged_into) {
		$log->write_to("NS->idMerge $wb_id into $merged_into\n") if $verbose;
		$nh->idMerge($wb_id, $merged_into);
	    } else {
		$log->write_to("NS->kill $wb_id\n") if $verbose;
		$nh->kill_gene($wb_id);
	    }
	    $nsp_out_fh->print($wb_id . "\n");
	}
	
	my $event = defined $merged_into ? 'Merged_into ' . $merged_into : 'Killed';
	
	
	$version++;
	$geneace_out_fh->print("Gene : $wb_id\n");
	$geneace_out_fh->print("Version $version\n");
	$geneace_out_fh->print("History Version_change $version now $person Event $event\n");
	$geneace_out_fh->print("Merged_into $merged_into\n") if defined $merged_into;
	$geneace_out_fh->print("Dead\n");
	if (scalar @references > 0) {
	    for my $ref (@references) {
		$geneace_out_fh->print("Remark \"$remark\" Paper_evidence $ref\n");
	    }
	    $geneace_out_fh->print("\n");
	} else {
	    $geneace_out_fh->print("Remark \"$remark\"\n\n");
	}

	$geneace_out_fh->print("Gene : $wb_id\n");
	$geneace_out_fh->print("-D Sequence_name\n");
	$geneace_out_fh->print("-D Method\n");
	$geneace_out_fh->print("-D Map_info\n");
	$geneace_out_fh->print("-D Other_name\n");
	$geneace_out_fh->print("-D Allele\n");
	$geneace_out_fh->print("-D Reference\n");
	$geneace_out_fh->print("-D Ortholog\n");
	$geneace_out_fh->print("-D Paralog\n\n");
    } else {
	$log->error("Gene $wb_id not found in curation database - cannot delete\n") if !$camace_gene;
	$log->error("Gene $wb_id not found in geneace - cannot delete\n") if !$geneace_gene;
    }
    
}

sub merge_genes {
    for my $external_id (keys %$to_merge) {
	my @wb_ids = @{$to_merge->{$external_id}};
	my @genes_to_merge;
	
	for my $gene_to_merge(@wb_ids) {
	    my $gene = $gdb->fetch(Gene => $gene_to_merge);
	    if (!$gene) {
		$log->error("Could not find Live gene $gene_to_merge in database to merge\n");
	    } else {
		push @genes_to_merge, $gene;
	    }
	}
	
	my $valid_merge = 1;
	my (%gene_starts, %gene_ends);
	for my $ix (0 .. scalar @genes_to_merge - 1) {
	    if ($wb_gene_models->{$genes_to_merge[0]->name}{'gene'}{'chromosome'} ne $wb_gene_models->{$genes_to_merge[$ix]->name}{'gene'}{'chromosome'}){
		$log->error("Attempting to merge genes from different chromosomes for $external_id\n");
		$valid_merge = 0;
	    }
	    
	    if ($wb_gene_models->{$genes_to_merge[0]->name}{'gene'}{'strand'} ne $wb_gene_models->{$genes_to_merge[$ix]->name}{'gene'}{'strand'}){
		$log->error("Attempting to merge genes from different strands for $external_id\n");
		$valid_merge = 0;
	    }
	    $gene_starts{$wb_gene_models->{$genes_to_merge[$ix]->name}{'gene'}{'start'}} = $genes_to_merge[$ix]->name;
	    $gene_ends{$wb_gene_models->{$genes_to_merge[$ix]->name}{'gene'}{'end'}} = $genes_to_merge[$ix]->name;
	}
	next unless $valid_merge;
	
	my $upstream_gene_name;
	if ($wb_gene_models->{$genes_to_merge[0]->name}{'gene'}{'strand'} eq '+') {
	    my @sorted_starts = sort {$a<=>$b} keys %gene_starts;
	    $upstream_gene_name = $gene_starts{$sorted_starts[0]};
	} else {
	    my @sorted_ends = sort {$a<=>$b} keys %gene_ends;
	    $upstream_gene_name = $gene_ends{$sorted_ends[-1]};
	}
	
	my $upstream_gene;
	my @downstream_genes;
	for my $gene_to_merge (@genes_to_merge) {
	    if ($gene_to_merge->name eq $upstream_gene_name) {
		$upstream_gene = $gene_to_merge;
	    } else {
		push @downstream_genes, $gene_to_merge;
	    }
	}
	
	gene_merge($upstream_gene, \@downstream_genes, $external_id);
    }
}

sub gene_merge {
    my ($upstream_gene, $downstream_genes, $external_id) = @_;
    
    my @acquired_merges;

    for my $downstream_gene (@$downstream_genes) {
	$log->write_to("Merging $downstream_gene into $upstream_gene\n") if $verbose;
	delete_gene($downstream_gene->name, $upstream_gene->name);
	
	# transfer operon connections.
	my $operons = $downstream_gene->at('Gene_info.Contained_in_operon');
	if (defined $operons) {
	    foreach my $operon ($operons->col) {
		$geneace_out_fh->print("Gene : $upstream_gene\nContained_in_operon $operon\n\n");
	    }
	}
	
	# transfer the Other_names
	my $other_names = $downstream_gene->at('Identity.Name.Other_name');
	if (defined $other_names) {
	    foreach my $other_name ($other_names->col) {
		$geneace_out_fh->print("Gene : $upstream_gene\nOther_name $other_name\n\n");
	    }
	}	
	
	# transfer alleles
	for my $allele ($downstream_gene->at('Gene_info.Allele')) {
	    if (defined $allele) {
		$geneace_out_fh->print("Gene : $upstream_gene\nAllele $allele\n\n");
	    }
	}
	
	# transfer references
	for my $reference ($downstream_gene->at('Reference')) {
	    if (defined $reference) {
		$geneace_out_fh->print("Gene : $upstream_gene\nReference $reference\n\n");
	    }
	}
	# transfer Features
	for my $feature ($downstream_gene->at('Associated_feature')) {
	    if (defined $feature) {
		$geneace_out_fh->print("Gene : $upstream_gene\nAssociated_feature $feature\n\n");
	    }
	}
	
	# transfer the Ortholog tags
	for my $ortholog ($downstream_gene->at('Gene_info.Ortholog')) {
	    if (defined $ortholog) {
		my @row = $ortholog->row;
		my $row0 = $row[0]->name;
		$log->log_and_die("\n\n!!!!!!!Missing ortholog species data in gene $downstream_gene\n\n") unless defined $row[1];
		my $row1 = $row[1]->name;
		$log->log_and_die("\n\n!!!!!!!Missing ortholog Evidence data in gene $downstream_gene - $row0($row1)\n\n") unless defined $row[2];
		my $row2 = $row[2]->name;
		my @col = $downstream_gene->at("Gene_info.Ortholog.$row0.$row1.$row2")->col;
		$row1 = '"' . $row1 . '"'; # add quotes to the species name
		foreach my $col (@col) {
		    $geneace_out_fh->print("Gene : $upstream_gene\nOrtholog $row0 $row1 $row2 $col\n\n");
		}
	    }
	}
	push @acquired_merges, $downstream_gene->name;
    }
    update_gene($upstream_gene->name, $external_id, undef, \@acquired_merges);
}


sub remove_cgc_name {
    my $wb_id = shift;

#    for my $wb_id (@changing_wb_genes) {
	my $wb_gene = $gdb->fetch(Gene => $wb_id);
	if (!defined $wb_gene) {
	    $log->error("Could not retrieve gene $wb_id to remove CGC name\n");
	    next;
	}
	next unless $wb_gene->CGC_name;
	my $cgc_name = $wb_gene->CGC_name->name;
	$cgc_names{$wb_id} = $cgc_name;
	$log->write_to("Removing CGC name $cgc_name for $wb_id\n") if $verbose;
	my $version = $wb_gene->Version->name;
	$version++;
	$geneace_out_fh->print("Gene : $wb_id\n");
	$geneace_out_fh->print('Public_name ' . $wb_gene->Sequence_name->name . "\n");
	$geneace_out_fh->print("Version $version\n");
	$geneace_out_fh->print("History Version_change $version now $person Name_change CGC_name \"Removed\"\n\n");
	if ($nameserver && !exists $ns_processed{$wb_id}) {
	    $log->write_to("NS->delName $cgc_name\n") if $verbose;
	    $nh->delName($cgc_name);
	    $nsp_out_fh->print($wb_id . "\n");
	}
	$geneace_out_fh->print("Gene : $wb_id\n");
	$geneace_out_fh->print("-D CGC_name\n\n");
#    }
}


sub split_genes {
    for my $wb_id (keys %$to_split) {
	my @new_ids = @{$to_split->{$wb_id}};
	
	my $valid_split = 1;
	my (%gene_starts, %gene_ends);
	for my $ix (0 .. scalar @new_ids - 1) {
	    if ($new_gene_models->{$new_ids[0]}{'gene'}{'chromosome'} ne $new_gene_models->{$new_ids[$ix]}{'gene'}{'chromosome'}) {
		$log->error("Attempting to split $wb_id into genes on different chromosomes\n");
		$valid_split = 0;
	    }
	    if ($new_gene_models->{$new_ids[0]}{'gene'}{'strand'} ne $new_gene_models->{$new_ids[$ix]}{'gene'}{'strand'}) {
		$log->error("Attempting to split $wb_id into genes on different strands\n");
		$valid_split = 0;
	    }
	    $gene_starts{$new_gene_models->{$new_ids[$ix]}{'gene'}{'start'}} = $new_ids[$ix];
	    $gene_ends{$new_gene_models->{$new_ids[$ix]}{'gene'}{'end'}} = $new_ids[$ix];
	    
	}
	next unless $valid_split;
	
	my $upstream_gene_name;
	if ($new_gene_models->{$new_ids[0]}{'gene'}{'strand'} eq '+') {
	    my @sorted_starts = sort {$a<=>$b} keys %gene_starts;
	    $upstream_gene_name = $gene_starts{$sorted_starts[0]};
	} else {
	    my @sorted_ends = sort {$a<=>$b} keys %gene_ends;
	    $upstream_gene_name = $gene_ends{$sorted_ends[-1]};
	}
	
	my $upstream_gene;
	my @downstream_genes;
	for my $split_gene (@new_ids) {
	    if ($split_gene eq $upstream_gene_name) {
		$upstream_gene = $split_gene;
	    } else {
		push @downstream_genes, $split_gene;
	    }
	}
	
	split_gene($wb_id, $upstream_gene, \@downstream_genes);
    }
}

sub split_gene {
    my ($wb_id, $upstream_id, $downstream_ids) = @_;
    
    my $wb_gene = $gdb->fetch(Gene => $wb_id);
    if ($wb_gene) {
	$log->write_to("Updating $wb_id by splitting\n") if $verbose;
	my @split_into;
	for my $downstream_id (@$downstream_ids) {
	    push @split_into, create_gene($downstream_id, $wb_id);
	}
	update_gene($wb_id, $upstream_id, \@split_into);
    } else {
	$log->error("Cannot split $wb_id as gene doesn't exist in database\n");
    }
}

sub update_genes {
    for my $wb_id (keys %$to_update) {
	update_gene($wb_id, $to_update->{$wb_id});
    }
}

sub update_gene {
    my ($wb_id, $external_id, $split_into, $acquired_merges) = @_;

    my $geneace_gene = $gdb->fetch(Gene => $wb_id);
    my $camace_gene = $cdb->fetch(Gene => $wb_id);
    
    $log->error("Cannot update $wb_id as gene doesn't exist in geneace\n") if !$geneace_gene;
    $log->error("Cannot update $wb_id as gene doesn't exist in curation database\n") if !$camace_gene;
    if ($geneace_gene && $camace_gene) {
	$log->write_to("Updating $wb_id\n") if $verbose;
	my $wb_name = $geneace_gene->Sequence_name->name;
	
	# Some convoluted logic to avoid deleting and recreating unchanged isoforms and to resuse / create transcripts and CDS names

	# Create map of CDS coords to ID for new models
	my (%new_coords);
	for my $nt (keys %{$new_gene_models->{$external_id}{'transcripts'}}) {
	    my $cds_coord_summary = coord_summary($external_id, 'CDS', $nt, $new_gene_models);
	    $new_coords{$cds_coord_summary} = $nt;
	}

	# Create map of CDS coords to ID for existing models
	my (%existing_coords);
	for my $ccds ($camace_gene->Corresponding_CDS) {
	    for my $transcript_id (keys %{$wb_cds_transcript_map{$ccds->name}}) {
		my $cds_coord_summary = coord_summary($wb_id, 'CDS', $transcript_id, $wb_gene_models);
		$existing_coords{$cds_coord_summary} = $ccds->name;
	    }
	}

	my (%cds_to_keep, %cds_to_update);
	my (@cds_to_delete, @cds_to_create);

	# Iterate through existing models
	# If new model with matching CDS coords in new set then add ID to %cds_to_keep
	# If no matching model then add ID to @cds_to_delete
	for my $existing_cds (keys %existing_coords) {
	    if (exists $new_coords{$existing_cds}) {
		$cds_to_keep{$existing_coords{$existing_cds}} = $new_coords{$existing_cds};
	    } else {
		push @cds_to_delete, $existing_coords{$existing_cds};
	    }
	}

	# Iterate through new models
	# If no matching old model then attempt to take an ID from @cds_to_delete list and add to %cds_to_update
	# If no entries in @cds_to_delete then add to @cds_to_create
	for my $new_cds (keys %new_coords) {
	    if (!exists $existing_coords{$new_cds}) {
		if (scalar @cds_to_delete > 0) {
		    my $to_update = shift @cds_to_delete;
		    $cds_to_update{$to_update} = $new_coords{$new_cds};
		} else {
		   push @cds_to_create, $new_coords{$new_cds};
		}
	    }
	}

	if (scalar @cds_to_create > 0 && scalar @cds_to_delete > 0) {
	    $log->log_and_die("Something has gone wrong - trying to create and delete CDS when should be updating existing\n");
	}

	my @histories_to_create = @cds_to_delete;
	push @histories_to_create, keys %cds_to_update;

	for my $history_cds (@histories_to_create) {
	    make_history_cds($history_cds, 'CDS deprecated as part of bulk annotation update on ' . $date);
	}
	
	my $version = $geneace_gene->Version->name;
	$log->error("Lower gene version for $wb_id in geneace than in curation DB\n") if $camace_gene->Version && $version < $camace_gene->Version->name;
	
	$version++;
	$version++ if exists $cgc_names{$wb_id};
        if ($split_into || $acquired_merges) {
	    
	    $geneace_out_fh->print("Gene : $wb_id\n");
	    $geneace_out_fh->print("Version $version\n");
	    if ($split_into) {
		$geneace_out_fh->print("History Version_change $version now $person Event Split_into $split_into->[0]\n");
		if (scalar @$split_into > 1) {
		    for my $ix (1 .. scalar @$split_into - 1) {
			$geneace_out_fh->print("Split_into $split_into->[$ix]\n");
		    }
		}
	    } else {
		$geneace_out_fh->print("History Version_change $version now $person Event Acquires_merge $acquired_merges->[0]\n");
		if (scalar @$acquired_merges > 1) {
		    for my $ix (1 .. scalar @$acquired_merges - 1) {
			$geneace_out_fh->print("Acquires_merge $acquired_merges->[$ix]\n");
		    }
		}
	    }
	    $geneace_out_fh->print("\n");
	}

	my $multiple_isoforms = scalar @cds_to_create + scalar (keys %cds_to_update) + scalar (keys %cds_to_keep) > 1 ? 1 : 0;
	my %cds_ids_created;

	if ($multiple_isoforms) {
	    for my $existing_cds (keys %cds_to_keep) {
		if ($existing_cds =~ /\d$/) {
		    my $wb_cds_name = $existing_cds;
		    make_history_cds($existing_cds, 'CDS deprecated as part of bulk annotation update on ' . $date);
		    for my $suffix (@ISOFORM_SUFFIXES) {
			next if exists $cds_ids_created{$wb_cds_name . $suffix};
			$wb_cds_name .= $suffix;
			$cds_ids_created{$wb_cds_name}++;
			last;
		    }
		    create_cds($external_id, $cds_to_keep{$existing_cds}, $wb_cds_name, $wb_id);
		}
	    }
	}
	
	for my $new_cds_id (@cds_to_create) {
	    my $wb_cds_name = $geneace_gene->Sequence_name->name;
	    if ($multiple_isoforms) {
		for my $suffix (@ISOFORM_SUFFIXES) {
		    next if exists $cds_to_update{$wb_cds_name . $suffix};
		    next if exists $cds_ids_created{$wb_cds_name . $suffix};
		    $wb_cds_name .= $suffix;
		    $cds_ids_created{$wb_cds_name}++;
		    last;
		}
	    }
	    create_cds($external_id, $new_cds_id, $wb_cds_name, $wb_id);
	}

	for my $existing_cds (keys %cds_to_update) {
	    if ($multiple_isoforms && $existing_cds =~ /\d$/) {
		my $wb_cds_name = $existing_cds;
		for my $suffix (@ISOFORM_SUFFIXES) {
		    next if exists $cds_ids_created{$wb_cds_name . $suffix};
		    $wb_cds_name .= $suffix;
		    $cds_ids_created{$wb_cds_name}++;
		    last;
		}
		create_cds($external_id, $cds_to_update{$existing_cds}, $wb_cds_name, $wb_id);
	    } else {
		create_cds($external_id, $cds_to_update{$existing_cds}, $existing_cds, $wb_id);
	    }
	}
    }
}


sub make_history_cds {
    my ($cds, $remark) = @_;
    
    my $version = $curation->get_history_version($wb->database('current'));
    my $wormpep_prefix = lc($wb->wormpep_prefix);
    my $history = "$cds:${wormpep_prefix}$version";

    $log->write_to("Creating history CDS $history\n") if $verbose;
    my $obj = $cdb->fetch(CDS => "$cds");
    
    if (!defined $obj) {
	$log->error("$cds is not a valid CDS/Pseudogene/Transcript name\n");
	return;
    }
    if ($cdb->fetch(CDS => $history)) {
	$log->warn("History object $history already exists\n");
	return;
    }

    my $gene = $obj->Gene->name;
    my $seq = $obj->Sequence->name;
    my $clone = $obj->Sequence;
    my @clone_entries = $clone->CDS_child;

    my ($start, $end);
    for my $clone_entry (@clone_entries) {
	next unless $clone_entry->name eq $cds;
	$start = $clone_entry->right->name;
	$end = $clone_entry->right->right->name;
	last;
    }

    $curation_out_fh->print("Sequence : $seq\n");
    $curation_out_fh->print("CDS_child \"$history\" $start $end\n\n");
    
    # Update existing CDS
    $curation_out_fh->print("CDS : $cds\n");
    $curation_out_fh->print("Evidence\n");
    $curation_out_fh->print("Method history\n");
    $curation_out_fh->print("Gene_history $gene\n");
    if (scalar @references > 0) {
	for my $ref (@references) {
	    $curation_out_fh->print("Remark \"$remark\" Paper_evidence $ref\n");
	}
    } else {
	 $curation_out_fh->print("Remark \"$remark\"\n");
    }
    # Remove Gene tag
    $curation_out_fh->print("\nCDS : $cds\n");
    $curation_out_fh->print("-D Gene\n\n");
    # Rename CDS
    $curation_out_fh->print("-R CDS $cds $history\n\n");
    return;
}


sub parse_gff {
    my ($file, $is_wormbase_gff) = @_;

    $log->write_to("Parsing GFF file $file\n");
    my %gene_models;
    my %transcript_gene_map;
    my $fh = file($file)->openr;
    while (my $line = $fh->getline()) {
	next if $line =~ /^#/;
	chomp $line;
	my @col = split("\t", $line);
	next if $col[1] eq 'history';
	my %attr = split(/[=;]/, $col[8]);
	if ($col[2] eq 'gene') {
	    my $gene_id = $attr{'ID'};
	    $gene_id =~ s/^Gene://;
	    $gene_models{$gene_id}{'gene'}{'chromosome'} = $col[0];
	    $gene_models{$gene_id}{'gene'}{'start'} = $col[3];
	    $gene_models{$gene_id}{'gene'}{'end'} = $col[4];
	    $gene_models{$gene_id}{'gene'}{'strand'} = $col[6];
	} elsif ($col[2] eq 'mRNA' || $col[2] eq 'pseudogenic_transcript' || $col[2] eq 'ncRNA') {
	    my $parent_gene_id = $attr{'Parent'};
	    $parent_gene_id =~ s/^Gene://;
	    my $transcript_id = $attr{'ID'};
	    $transcript_id =~ s/^Transcript://;
	    $transcript_id =~ s/^Pseudogene://;
	    $gene_models{$parent_gene_id}{'transcripts'}{$transcript_id}{'start'} = $col[3];
	    $gene_models{$parent_gene_id}{'transcripts'}{$transcript_id}{'end'} = $col[4];
	    $transcript_gene_map{$transcript_id} = $parent_gene_id;
	} elsif ($col[2] eq 'exon' || $col[2] eq 'CDS') {
	    my $parent_transcript_id = $attr{'Parent'};
	    $parent_transcript_id =~ s/Transcript://g;
	    $parent_transcript_id =~ s/Pseudogene://g;
	    my @parent_transcript_ids = split(",", $parent_transcript_id);
	    for my $pt_id (@parent_transcript_ids) {
		my $gene_id = $transcript_gene_map{$pt_id};
		if ($is_wormbase_gff) {
		    my ($cds_id) = $pt_id =~ /^(.+)\.\d+$/;
		    $cds_id = $pt_id if !defined $cds_id;
		    $wb_cds_transcript_map{$cds_id}{$pt_id}++;
		}
		$log->log_and_die("Can't find gene ID for $col[2] with parent $pt_id\n") unless defined $gene_id;
		push @{$gene_models{$gene_id}{'transcripts'}{$pt_id}{$col[2]}{'starts'}}, $col[3];
		push @{$gene_models{$gene_id}{'transcripts'}{$pt_id}{$col[2]}{'ends'}}, $col[4];
	    }
	} else {
	    # do nothing
	}
    }
    
    return \%gene_models;
}

sub get_relative_positions {
    my ($gene_id, $transcript_id, $feature_type) = @_;

    my @starts = sort {$a<=>$b} @{$new_gene_models->{$gene_id}{'transcripts'}{$transcript_id}{$feature_type}{'starts'}};
    my @ends = sort {$a<=>$b} @{$new_gene_models->{$gene_id}{'transcripts'}{$transcript_id}{$feature_type}{'ends'}};
    
    my (@relative_starts, @relative_ends);
    for my $ix (0 .. scalar @starts - 1) {
	my ($relative_start, $relative_end);
	if ($new_gene_models->{$gene_id}{'gene'}{'strand'} eq '-') {
	    my $reverse_ix = scalar(@starts) - ($ix + 1);
	    $relative_start = ($ends[scalar @starts - 1] - $ends[$reverse_ix]) + 1;
	    $relative_end = ($ends[scalar @starts - 1] - $starts[$reverse_ix]) + 1;
	} else {
	    $relative_start = ($starts[$ix] - $starts[0]) + 1;
	    $relative_end = ($ends[$ix] - $starts[0]) + 1;
	}
	push @relative_starts, $relative_start;
	push @relative_ends, $relative_end;
    }
    
    return (\@relative_starts, \@relative_ends);
}

sub coord_summary {
    my ($gene_id, $feature_type, $transcript_id, $gene_models) = @_;
    
    my $coord_type = $feature_type eq 'transcript' ? 'exon' : 'CDS';
    my @starts = sort {$a<=>$b} @{$gene_models->{$gene_id}{'transcripts'}{$transcript_id}{$coord_type}{'starts'}};
    my @ends = sort {$a<=>$b} @{$gene_models->{$gene_id}{'transcripts'}{$transcript_id}{$coord_type}{'ends'}};
    my @coords;
    for my $ix (0 .. @starts - 1) {
	push @coords, $starts[$ix] . '-' . $ends[$ix];
    }
    
    return join("|", @coords);
}


sub parse_mapping_file {
    my ($file, $list_values) = @_;

    $log->write_to("Parsing mapping file $file\n");
    # Mapping file should contain two tab-delimited columns
    # The first column should be a single gene
    # The second should be multiple genes separated by commas
    my %map;
    my $fh = file($file)->openr;
    while (my $line = $fh->getline()) {
	chomp $line;
	next if $line =~ /^#/;
#	my @genes_involved;
	my @col = split("\t", $line);
#	push @genes_involved, $col[0];
#	push @genes_involved, split(',', $col[1]);
	if ($list_values) {
	    push @{$map{$col[0]}}, split(',', $col[1]);
	} else {
	    $map{$col[0]} = $col[1];
	}

#	for my $gene_involved (@genes_involved) {
#	    push @changing_wb_genes, $gene_involved if $gene_involved =~ /^WBGene\d+$/;
#	}
    }

    return \%map;
}

sub parse_single_column_file {
    my $file = shift;
    
    $log->write_to("Parsing file $file\n");
    my @lines;
    my $fh = file($file)->openr;
    while (my $line = $fh->getline()) {
	chomp $line;
	next if $line =~ /^#/;
	push @lines, $line;
    }
    
    return \@lines;
}

