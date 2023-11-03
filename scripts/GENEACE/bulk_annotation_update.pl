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
    
my ($create_file, $change_file, $delete_file, $wb_gff_file, $gff_file, $out_file, $test, $log_file, $species, $debug, $database, $nameserver, $person, $verbose);

GetOptions(
    "create:s"   => \$create_file,
    "delete:s"   => \$delete_file,
    "change:s"   => \$change_file,
    "gff=s"      => \$gff_file,
    "wbgff=s"    => \$wb_gff_file,
    "out=s"      => \$out_file,
    "database=s" => \$database,
    "species:s"  => \$species,
    "test"       => \$test,
    "verbose"    => \$verbose,
    "log=s"      => \$log_file,
    "person=s"   => \$person,
    "debug:s"    => \$debug,
    "nameserver" => \$nameserver
    ) or die ($@);

$out_file = file('./bulk_update.ace') unless defined $out_file;
my $out_fh = file($out_file)->openw;

my $log = Log_files->make_log($log_file, $debug);
my $wb = Wormbase->new(-organism => $species, -debug => $debug, -test => $test);
my $db = Ace->connect(-path => $database);
my $nh = NameDB_handler->new($wb, $test);
my $curation = Curate->new($out_fh, $db, $wb, $log);

if ($species eq 'elegans' || $species eq 'sratti' || $species eq 'tmuris') {
    $log->log_and_die("This script is not intended for $species updates - see comments in code\n");
}

my $test_count = 0;

my ($day, $mon, $yr)  = (localtime)[3,4,5];
my $date = sprintf("%02d%02d%02d",$yr-100, $mon+1, $day);

my $wsversion = `grep "NAME WS" $database/wspec/database.wrm`;
chomp($wsversion);
$wsversion =~ s/.*WS//;
$log->log_and_die("No version\n") unless $wsversion =~ /^\d+$/;

my $wormpep_prefix = lc($wb->wormpep_prefix);

my $new_gene_models = parse_gff($gff_file);
my $wb_gene_models = parse_gff($wb_gff_file);

my ($to_update, $to_merge, $to_split, $to_create, $to_delete) = parse_mapping_file($change_file);

check_for_merge_split_clashes();

my %wb_ids_created;

$log->write_to("Creating genes\n");
create_genes() if defined $create_file;

$log->write_to("Deleting genes\n");
delete_genes() if defined $delete_file;

$log->write_to("Creating genes to replace existing genes with CGC names\n");
create_cgc_replacement_genes() if scalar keys %$to_create > 0;

$log->write_to("Deleting genes with CGC names\n");
delete_cgc_genes() if scalar keys %$to_delete > 0;

$log->write_to("Splitting genes\n");
split_genes() if scalar keys %$to_split > 0;

$log->write_to( "Merging genes\n");
merge_genes() if scalar keys %$to_merge > 0;

print "Updating genes\n";
update_genes() if scalar keys %$to_update > 0;

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

sub create_cgc_replacement_genes {
    for my $external_id (keys %$to_create) {
	create_gene($external_id, undef, $to_create->{$external_id});
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
    my ($non_wb_id, $split_from, $replaced_cgc_genes) = @_;
    
    # If this code is ever to be used for elegans, sratti, or tmuris updates then we would need to supply a clone ID to Next_CDS_ID here
    my $seq_name = $curation->Next_CDS_ID();
    my $so_term = 'SO:0001217'; # Are we really only creating coding genes or could we be creating Pseudogenes?
    my $wb_id = $nh->idGetByTypedName('Sequence' => $seq_name);
    if ($wb_id) {
	$log->write_to("$wb_id already exists in Names Service, creating using existing ID\n");
    } else {
	my $new_gene = $nh->new_gene($seq_name, 'Sequence', $wb->full_name);
	if (defined $new_gene) {
	    $wb_id = $new_gene->{'id'};
	} else {
	    $log->error("Could not retrieve new gene ID from Names Service for $seq_name\n");
	    return;
	}
    }

    if ($verbose) {
	if (defined $replaced_cgc_genes) {
	    $log->write_to("Creating new Gene $wb_id for $seq_name to replace $replaced_cgc_genes\n");
	} else {
	    $log->write_to("Creating new Gene $wb_id for $seq_name\n");
	}
    }
    
    my $event = defined $split_from ? 'Split_from ' . $split_from : 'Created';
    
    # This code comes from GENEACE/newgene.pl - additional code required if this script is ever to be used for elegans
    $out_fh->print("Gene : $wb_id\n");
    $out_fh->print("Live\n");
    $out_fh->print("Version 1\n");
    $out_fh->print("Biotype \"$so_term\"\n");
    $out_fh->print("Sequence_name \"$seq_name\"\n");
    $out_fh->print("Species \"${\$wb->full_name}\"\n");
    $out_fh->print("History Version_change 1 now $person Event $event\n");
    $out_fh->print("Split_from $split_from\n") if defined $split_from;
    $out_fh->print("Method Gene\n");
    $out_fh->print("Remark \"Updated model for deprecated gene(s) $replaced_cgc_genes as part of bulk annotation update on $date - new gene created due to existing genes having CGC names\"\n") if defined $replaced_cgc_genes;
    $out_fh->print("Public_name \"$seq_name\"\n\n");
    
    my $ix = 0;
    for my $external_transcript_id(keys %{$new_gene_models->{$non_wb_id}{'transcripts'}}) {
	$ix++;
	my ($wb_transcript_name, $wb_cds_name);
	if (scalar keys %{$new_gene_models->{$non_wb_id}{'transcripts'}} == 1) {
	    $wb_transcript_name = $seq_name . '.1';
	    $wb_cds_name = $seq_name;
	} else {
	    $wb_transcript_name = $seq_name . $ISOFORM_SUFFIXES[$ix] . '.1';
	    $wb_cds_name = $seq_name . $ISOFORM_SUFFIXES[$ix];
	}
	create_transcript($non_wb_id, $external_transcript_id, $wb_transcript_name, $wb_id);
	create_cds($non_wb_id, $external_transcript_id, $wb_cds_name, $wb_transcript_name, $wb_id);
    }
    $wb_ids_created{$non_wb_id} = $wb_id;
    return $wb_id;
}

sub create_transcript {
    my ($external_gene_id, $external_transcript_id, $wb_transcript_name, $wb_gene_id) = @_;
    
    my ($exon_starts, $exon_ends) = get_relative_positions($external_gene_id, $external_transcript_id, 'exon');
    
    my $remark = 'Created as part of bulk annotation update on ' . $date;
    my $sequence = $new_gene_models->{$external_gene_id}{'gene'}{'chromosome'};
    
    my ($transcript_start, $transcript_end);
    if ($new_gene_models->{$external_gene_id}{'gene'}{'strand'} eq '-') {
	$transcript_start = $new_gene_models->{$external_gene_id}{'transcripts'}{$external_transcript_id}{'end'};
	$transcript_end = $new_gene_models->{$external_gene_id}{'transcripts'}{$external_transcript_id}{'start'};
    } else {
	$transcript_start = $new_gene_models->{$external_gene_id}{'transcripts'}{$external_transcript_id}{'start'};
	$transcript_end = $new_gene_models->{$external_gene_id}{'transcripts'}{$external_transcript_id}{'end'};
    }
    
    my $method = 'Coding_transcript'; # Should be valid for all initial use cases, but need to reassess
    
    $out_fh->print("Transcript : \"$wb_transcript_name\"\n");
    $out_fh->print("Sequence $sequence\n");
    $out_fh->print("Species \"${\$wb->full_name}\"\n");
    $out_fh->print("Gene $wb_gene_id\n");
    $out_fh->print("Method $method\n");
    $out_fh->print("Remark \"$remark\"\n");
    for my $ix(0 .. scalar @$exon_starts - 1) {
	$out_fh->print("Source_exons $exon_starts->[$ix] $exon_ends->[$ix]\n");
    }
    $out_fh->print("\nSequence : \"$sequence\"\n" );
    $out_fh->print("Transcript \"$wb_transcript_name\" $transcript_start $transcript_end\n\n");
}

sub create_cds {
    my ($external_gene_id, $external_transcript_id, $wb_cds_name, $wb_transcript_name, $wb_gene_id) = @_;
    
    my ($cds_starts, $cds_ends) = get_relative_positions($external_gene_id, $external_transcript_id, 'CDS');
    
    my @starts_on_sequence = @{$new_gene_models->{$external_gene_id}{'transcripts'}{$external_transcript_id}{'CDS'}{'starts'}};
    my @ends_on_sequence = @{$new_gene_models->{$external_gene_id}{'transcripts'}{$external_transcript_id}{'CDS'}{'ends'}};
    my $remark = 'Created as part of bulk annotation update on ' . $date;
    my $sequence = $new_gene_models->{$external_gene_id}{'gene'}{'chromosome'};
    my $method = 'curated'; # Going with this for now, what's an appropriate method?
    
    $out_fh->print("CDS : \"$wb_cds_name\"\n");
    $out_fh->print("Gene : $wb_gene_id\n");
    $out_fh->print("Sequence $sequence\n");
    $out_fh->print("Species \"${\$wb->full_name}\"\n");
    $out_fh->print("Corresponding_transcript $wb_transcript_name\n");
    $out_fh->print("Method $method\n");
    $out_fh->print("Remark \"$remark\"\n");
    for my $ix(0 .. scalar @$cds_starts - 1) {
	$out_fh->print("Source_exons $cds_starts->[$ix] $cds_ends->[$ix]\n");
    }
    $out_fh->print("\nSequence : \"$sequence\"\n" );
    for my $ix (0 .. scalar @starts_on_sequence - 1) {
	my $reverse_ix = scalar(@starts_on_sequence) - ($ix + 1);
	my ($start_on_sequence, $end_on_sequence);
	if ($new_gene_models->{$external_gene_id}{'gene'}{'strand'} eq '-') {
	    $start_on_sequence = $ends_on_sequence[$reverse_ix];
	    $end_on_sequence = $starts_on_sequence[$reverse_ix];
	} else {
	    $start_on_sequence = $starts_on_sequence[$ix];
	    $end_on_sequence = $ends_on_sequence[$ix];
	}
	
	$out_fh->print("CDS_child \"$wb_cds_name\" $start_on_sequence $end_on_sequence\n");
    }
    $out_fh->print("\n");
}

sub delete_cgc_genes {
    for my $wb_name (keys %$to_delete) {
	delete_gene($wb_name, undef, $to_delete->{$wb_name});
    }

    return;
}

sub delete_genes {
    my $wb_genes_to_delete = parse_single_column_file($delete_file);
    for my $wb_name (@{$wb_genes_to_delete}) {
	delete_gene($wb_name);
    }
    
    return;
}

sub delete_gene {
    my ($wb_name, $merged_into, $replaced_with) = @_;
    
    my $remark = 'Removed as part of bulk annotation update on ' . $date;
    my $replacement_wb_id;
    if ($replaced_with) {
	$replacement_wb_id = $wb_ids_created{$replaced_with};
	$remark .= " - replaced with $replacement_wb_id (update not possible due to CGC name)";
    }
    
    my @genes = $db->fetch(-query => "Find Gene WHERE Sequence_name = \"$wb_name\"");
    my $gene = $genes[0];
    if ($gene) {
	my $wb_id = $gene->name;
	if ($verbose) {
	    if ($replacement_wb_id) {
		$log->write_to("Deleting gene $wb_id (replaced with $replacement_wb_id)\n");
	    } else {
		$log->write_to("Deleting gene $wb_id\n") if $verbose;
	    }
	    $log->write_to("NS->kill $wb_id\n") if $nameserver;
	}
	$nh->kill_gene($wb_id) if $nameserver;
	
	my $event = defined $merged_into ? 'Merged_into ' . $merged_into : 'Killed';
	
	my $version = $gene->Version->name;
	$version++;
	$out_fh->print("Gene : $wb_id\n");
	$out_fh->print("Version $version\n");
	$out_fh->print("History Version_change $version now $person Event $event\n");
	$out_fh->print("Merged_into $merged_into\n") if defined $merged_into;
	$out_fh->print("Dead\n");
	$out_fh->print("Remark \"$remark\" Curator_confirmed $person\n\n");

	$out_fh->print("Gene : $wb_id\n");
	$out_fh->print("-D Sequence_name\n");
	$out_fh->print("-D Method\n");
	$out_fh->print("-D Map_info\n");
	$out_fh->print("-D Other_name\n");
	$out_fh->print("-D Allele\n");
	$out_fh->print("-D Reference\n");
	$out_fh->print("-D Ortholog\n");
	$out_fh->print("-D Paralog\n\n");
    } else {
	$log->error("Gene $wb_name not found in database - cannot delete\n");
    }
    
}

sub merge_genes {
    for my $external_id (keys %$to_merge) {
	my @wb_names = @{$to_merge->{$external_id}};
	my @genes_to_merge;
	
	for my $gene_to_merge(@wb_names) {
	    my @genes = $db->fetch(-query => "Find Gene WHERE Sequence_name = \"$gene_to_merge\" AND Live");
	    if (scalar @genes < 1) {
		$log->error("Could not find Live gene $gene_to_merge in database to merge\n");
	    } else {
		push @genes_to_merge, $genes[0];
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
	delete_gene($downstream_gene->Sequence_name->name, $upstream_gene->name);
	
	# transfer operon connections.
	my $operons = $downstream_gene->at('Gene_info.Contained_in_operon');
	if (defined $operons) {
	    foreach my $operon ($operons->col) {
		$out_fh->print("Gene : $upstream_gene\nContained_in_operon $operon\n\n");
	    }
	}
	
	# transfer the Other_names
	my $other_names = $downstream_gene->at('Identity.Name.Other_name');
	if (defined $other_names) {
	    foreach my $other_name ($other_names->col) {
		$out_fh->print("Gene : $upstream_gene\nOther_name $other_name\n\n");
	    }
	}	
	
	# transfer alleles
	for my $allele ($downstream_gene->at('Gene_info.Allele')) {
	    if (defined $allele) {
		$out_fh->print("Gene : $upstream_gene\nAllele $allele\n\n");
	    }
	}
	
	# transfer references
	for my $reference ($downstream_gene->at('Reference')) {
	    if (defined $reference) {
		$out_fh->print("Gene : $upstream_gene\nReference $reference\n\n");
	    }
	}
	# transfer Features
	for my $feature ($downstream_gene->at('Associated_feature')) {
	    if (defined $feature) {
		$out_fh->print("Gene : $upstream_gene\nAssociated_feature $feature\n\n");
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
		    $out_fh->print("Gene : $upstream_gene\nOrtholog $row0 $row1 $row2 $col\n\n");
		}
	    }
	}
	push @acquired_merges, $downstream_gene->name;
    }
    update_gene($upstream_gene->Sequence_name->name, $external_id, undef, \@acquired_merges);
}


sub split_genes {
    for my $wb_name (keys %$to_split) {
	my @new_ids = @{$to_split->{$wb_name}};
	
	my $valid_split = 1;
	my (%gene_starts, %gene_ends);
	for my $ix (0 .. scalar @new_ids - 1) {
	    if ($new_gene_models->{$new_ids[0]}{'gene'}{'chromosome'} ne $new_gene_models->{$new_ids[$ix]}{'gene'}{'chromosome'}) {
		$log->error("Attempting to split $wb_name into genes on different chromosomes\n");
		$valid_split = 0;
	    }
	    if ($new_gene_models->{$new_ids[0]}{'gene'}{'strand'} ne $new_gene_models->{$new_ids[$ix]}{'gene'}{'strand'}) {
		$log->error("Attempting to split $wb_name into genes on different strands\n");
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
	
	split_gene($wb_name, $upstream_gene, \@downstream_genes);
    }
}

sub split_gene {
    my ($wb_name, $upstream_id, $downstream_ids) = @_;
    
    my @wb_genes = $db->fetch(-query => "Find Gene WHERE Sequence_name = \"$wb_name\"");
    my $wb_gene = $wb_genes[0];
    if ($wb_gene) {
	$log->write_to("Updating $wb_name by splitting\n") if $verbose;
	my @split_into;
	for my $downstream_id (@$downstream_ids) {
	    push @split_into, create_gene($downstream_id, "$wb_gene");
	}
	update_gene($wb_name, $upstream_id, \@split_into);
    } else {
	$log->error("Cannot split $wb_name as gene doesn't exist in database\n");
    }
}

sub update_genes {
    for my $wb_name (keys %$to_update) {
	update_gene($wb_name, $to_update->{$wb_name});
    }
}

sub update_gene {
    my ($wb_name, $external_id, $split_into, $acquired_merges) = @_;
    
    my @wb_genes = $db->fetch(-query => "Find Gene WHERE Sequence_name = \"$wb_name\"");
    my $wb_gene = $wb_genes[0];
    if ($wb_gene) {
	my $wb_id = $wb_gene->name;
	$log->write_to("Updating $wb_name ($wb_id)\n") if $verbose;
	
	# Some convoluted logic to avoid deleting and recreating unchanged isoforms and to resuse / create transcripts and CDS names
	my (%new_coords);
	for my $nt (keys %{$new_gene_models->{$external_id}{'transcripts'}}) {
	    my $exon_coord_summary = coord_summary($external_id, 'transcript', $nt, $new_gene_models);
	    my $cds_coord_summary = coord_summary($external_id, 'CDS', $nt, $new_gene_models);
	    $new_coords{$exon_coord_summary}{'transcript'} = $nt;
	    $new_coords{$exon_coord_summary}{'CDS_coords'} = $cds_coord_summary;
	}
	
	my (%existing_coords);
	for my $ct ($wb_gene->Corresponding_transcript) {
	    my $exon_coord_summary = coord_summary($wb_id, 'transcript', $ct->name, $wb_gene_models);
	    my $cds_coord_summary = coord_summary($wb_id, 'CDS', $ct->name, $wb_gene_models);
	    $existing_coords{$exon_coord_summary}{'transcript'} = $ct->name;
	    $existing_coords{$exon_coord_summary}{'CDS_name'} = $ct->Corresponding_CDS->name;
	    $existing_coords{$exon_coord_summary}{'CDS_coords'} = $cds_coord_summary;
	}
	
	
	my @transcripts_to_create_with_new_cds; # list of new transcript IDs
	my %transcripts_to_create_with_existing_cds; # map of new transcript IDs to existing CDS names
	my %existing_transcripts_with_new_cds; # map of existing transcript names to new transcript IDs
	for my $new_exon_cs (keys %new_coords) {
	    my $new_cds_cs = $new_coords{$new_exon_cs}{'CDS_coords'};
	    if (!exists $existing_coords{$new_exon_cs}) {
		my $cds_exists = 0;
		for my $existing_exon_cs (keys %existing_coords) {
		    if ($existing_coords{$existing_exon_cs}{'CDS_coords'} eq $new_cds_cs) {
			$transcripts_to_create_with_existing_cds{$new_coords{$new_exon_cs}{'transcript'}} = $existing_coords{$existing_exon_cs}{'CDS_name'};
			$cds_exists = 1;
			last;
		    }
		}
		push @transcripts_to_create_with_new_cds, $new_coords{$new_exon_cs}{'transcript'} if !$cds_exists;
	    } else {
		if ($existing_coords{$new_exon_cs}{'CDS_coords'} ne $new_coords{$new_exon_cs}{'CDS_coords'}) {
		    $existing_transcripts_with_new_cds{$existing_coords{$new_exon_cs}{'transcript'}} = $new_coords{$new_exon_cs}{'transcript'};
		}
	    }
	}
	my %existing_cdses_with_new_transcripts = map {$_ => 1} values %transcripts_to_create_with_existing_cds;
	
	my %existing_transcript_names = map {$_->name => 1} $wb_gene->Corresponding_transcript;
	for my $existing_exon_cs (keys %existing_coords) {
	    if (!exists $new_coords{$existing_exon_cs}) {
		make_history($existing_coords{$existing_exon_cs}{'transcript'}, 'Transcript deprecated as part of bulk annotation update on ' . $date);
		$existing_transcript_names{$existing_coords{$existing_exon_cs}{'transcript'}} = 0;
		make_history($existing_coords{$existing_exon_cs}{'CDS_name'}, 'CDS deprecated as part of bulk annotation update on ' . $date) unless exists $existing_cdses_with_new_transcripts{$existing_coords{$existing_exon_cs}{'CDS_name'}};
	    }
	}
	
	if (scalar keys %existing_transcripts_with_new_cds == 0 && scalar @transcripts_to_create_with_new_cds == 0 && scalar keys %transcripts_to_create_with_existing_cds == 0) {
	    $log->write_to("Nothing to update for $wb_gene\n");
	    return;
	}

	my $version = $wb_gene->Version->name;
	$version++;
        if ($split_into || $acquired_merges) {
	    
	    $out_fh->print("Gene : $wb_id\n");
	    $out_fh->print("Version $version\n");
	    if ($split_into) {
		$out_fh->print("History Version_change $version now $person Event Split_into $split_into->[0]\n");
		if (scalar @$split_into > 1) {
		    for my $ix (1 .. scalar @$split_into - 1) {
			$out_fh->print("Split_into $split_into->[$ix]\n");
		    }
		}
	    } else {
		$out_fh->print("History Version_change $version now $person Event Acquires_merge $acquired_merges->[0]\n");
		if (scalar @$acquired_merges > 1) {
		    for my $ix (1 .. scalar @$acquired_merges - 1) {
			$out_fh->print("Acquires_merge $acquired_merges->[$ix]\n");
		    }
		}
	    }
	    $out_fh->print("\n");
	}

	# Find transcript names made available by making history transcripts and that aren't going to be used for existing CDSes that need new transcripts
	my @available_transcript_names;
	for my $etn (keys %existing_transcript_names) {
	    my ($cds_name) = $etn =~ /^(.+)\.\d$/;
	    push @available_transcript_names, $etn if $existing_transcript_names{$etn} == 0 && !exists($existing_cdses_with_new_transcripts{$etn});
	}
	my $next_transcript_suffix_ix = scalar keys %existing_transcript_names;
	
	for my $existing_transcript (keys %existing_transcripts_with_new_cds) {
	    my ($cds_name) = $existing_transcript =~ /^(.+)\.\d$/;
	    create_cds($external_id, $existing_transcripts_with_new_cds{$existing_transcript}, $cds_name, $existing_transcript, $wb_id);
	}
	for my $new_transcript (@transcripts_to_create_with_new_cds) {
	    my $new_transcript_name;
	    if (scalar @available_transcript_names > 0) {
		$new_transcript_name = shift(@available_transcript_names);
	    } elsif ($next_transcript_suffix_ix == 0 && scalar @transcripts_to_create_with_new_cds == 1 && scalar keys %transcripts_to_create_with_existing_cds == 0) {
		$new_transcript_name = $wb_name . '.1';
	    } else {
		$new_transcript_name = $wb_name . $ISOFORM_SUFFIXES[$next_transcript_suffix_ix] . '.1';
		$next_transcript_suffix_ix++;
	    }
	    my ($cds_name) = $new_transcript_name =~ /^(.+)\.\d$/;
	    create_transcript($external_id, $new_transcript, $new_transcript_name, $wb_id);
	    create_cds($external_id, $new_transcript, $cds_name, $new_transcript_name, $wb_id);
	}
	
	for my $new_transcript (keys %transcripts_to_create_with_existing_cds) {
	    create_transcript($external_id, $new_transcript, $transcripts_to_create_with_existing_cds{$new_transcript} . '.1', $wb_id);
	}
	
    } else {
	$log->error("Cannot update $wb_name as gene doesn't exist in database\n");
    }
    
}


sub make_history {
    my ($cds, $remark) = @_;
    
    my $biotype;
    my $clone_tag; # where it lives under the Sequence object tags
    my $new_method; # the new method to give it
    my $transcript_type;
    my $transcript_type_value1;
    my $transcript_type_value2;
    my $pseudogene_type;
    
    my $obj = $db->fetch(CDS => "$cds");
    if ($obj) {
	$biotype = 'CDS';
	$clone_tag = 'CDS_child';
	$new_method = 'history';
	# Original script doesn't do non-curated CDSes - any reason why not?
    } elsif ($obj = $db->fetch(Pseudogene => "$cds")) {
	$biotype = 'Pseudogene';
	$clone_tag = 'Pseudogene';
	$new_method = 'history_pseudogene';
	$pseudogene_type = $obj->Type->name;
    } elsif ($obj = $db->fetch(Transcript => "$cds")) {
	$biotype = 'Transcript';
	$clone_tag = 'Transcript';
	$new_method = 'history_transcript';
	# Original script doesn't do coding transcripts - any reason why not?;
    }
    
    if (!defined $obj) {
	$log->error("$cds is not a valid CDS/Pseudogene/Transcript name\n");
	return;
    }
    
    my $gene = $obj->Gene->name;
    
    # Update existing CDS / Transcript / Pseudogene
    $out_fh->print("$biotype : $cds\n");
    $out_fh->print("Evidence\n");
    $out_fh->print("Method $new_method\n");
    $out_fh->print("Gene_history $gene\n");
    $out_fh->print("Remark \"$remark\"\n");
    # Remove Gene tag
    $out_fh->print("\n$biotype : $cds\n");
    $out_fh->print("-D Gene\n\n");
    # Rename CDS
    $out_fh->print("-R $biotype $cds $cds:${wormpep_prefix}$wsversion\n\n");
    return;
}


sub parse_gff {
    my $file = shift;

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
	} elsif ($col[2] eq 'mRNA' || $col[2] eq 'pseudogenic_transcript') {
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
    my $file = shift;

    $log->write_to("Parsing mapping file $file\n");
    # Mapping file should contain two or three tab-delimited columns
    # The first column should be the external gene identifier
    # The second (and in some cases third column) should be the corresponding WB gene sequence names
    my (%update_map, %update_by_merge_map, %update_by_split_map, %cgc_to_delete, %to_create);
    my $fh = file($file)->openr;
    while (my $line = $fh->getline()) {
	chomp $line;
	next if $line =~ /^\w*$/;
	my @col = split("\t", $line);
	my @filtered_col;
	my %seen;
	for my $ix (0 .. scalar @col -1) {
	    if ($ix == 0) {
		push @filtered_col, $col[$ix];
	    } else {
		my $col = $col[$ix];
		if ($col =~ /^(\D+\d+)\D$/) {
		    $col = $1;
		}
		push @filtered_col, $col unless $seen{$col};
		$seen{$col}++;
	    }
	}
	
	if (scalar @filtered_col == 2) {
	    my @genes = $db->fetch(-query => "Find Gene WHERE Sequence_name = \"$filtered_col[1]\"");
	    my $gene_to_update = $genes[0];
	    if (!defined $gene_to_update) {
		$log->error("Could not find Gene $filtered_col[1] to update with model for $filtered_col[0] - skipping\n");
		next;
	    }
	    if ($gene_to_update->CGC_name) {
		$log->write_to("Cannot update $gene_to_update as it has CGC name - will kill it and create new gene for $filtered_col[0]\n");
		$to_create{$filtered_col[0]} = $gene_to_update->name;
		$cgc_to_delete{$filtered_col[1]} = $filtered_col[0];
		next;
	    }
	    if (exists $update_map{$filtered_col[1]} || exists $update_by_split_map{$filtered_col[1]}) {
		push @{$update_by_split_map{$filtered_col[1]}}, $filtered_col[0];
		if (exists $update_map{$filtered_col[1]}) {
		    push @{$update_by_split_map{$filtered_col[1]}}, $update_map{$filtered_col[1]};
		    delete($update_map{$filtered_col[1]});
		}
	    } else {
		$update_map{$filtered_col[1]} = $filtered_col[0];
	    }
	} else {
	    my $external_id = shift(@filtered_col);
	    my $contains_cgc_gene = 0;
	    my @wb_ids_to_merge;
	    for my $wb_gene (@filtered_col) {
		my @genes = $db->fetch(-query => "Find Gene WHERE Sequence_name = \"$wb_gene\"");
		my $gene_to_merge = $genes[0];
		push @wb_ids_to_merge, $gene_to_merge->name;
		$contains_cgc_gene = 1 if $gene_to_merge->CGC_name;
	    }
	    if ($contains_cgc_gene) {
		my $genes_to_merge = join(", ", @wb_ids_to_merge);
		$log->write_to("Cannot merge $genes_to_merge as one or more has CGC name - will kill them and create new gene for $external_id\n");
		for my $col (@filtered_col) {
		    $cgc_to_delete{$col} = $external_id;
		}
		$to_create{$external_id} = $genes_to_merge;
	    } else {
		push @{$update_by_merge_map{$external_id}}, @filtered_col;
	    }
	}
    }
    
    return (\%update_map, \%update_by_merge_map, \%update_by_split_map, \%to_create, \%cgc_to_delete);
}

sub parse_single_column_file {
    my $file = shift;
    
    $log->write_to("Parsing file $file\n");
    my @lines;
    my $fh = file($file)->openr;
    while (my $line = $fh->getline()) {
	chomp $line;
	next if $line =~ /^\s*$/;
	push @lines, $line;
    }
    
    return \@lines;
}

