#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;
use Log_files;
use Wormbase;
use Storable;
use Path::Class;
use Const::Fast;

my ($current_release, $current_iteration, $previous_release, $previous_iteration);
my ($database, $species, $log_file, $store, $debug);
GetOptions(
    "curr_release|r1=s"   => \$current_release,   
    "prev_release|r2=s"   => \$previous_release,   
    "curr_iter|i1=i"      => \$current_iteration,  # default is highest
    "prev_iter|i2=i"      => \$previous_iteration, # default is highest
    "database|d=s"        => \$database,
    "species|s=s"         => \$species,
    "logfile|l=s"         => \$log_file,
    "store|st=s"          => \$store,
    "debug|d=s"           => \$debug
    );

my $wormbase;
if ($store) {
    $wormbase = retrieve($store) or croak("Can't restore WormBase from store\n");
}
else {
    $wormbase = Wormbase->new(
	-debug    => $debug,
	-organism => $species,
	-autoace  => $database
	);
}

my $log = $log_file ? Log_files->make_log($log_file, $debug) : Log_files->make_build_log($wormbase);

# Default behaviour is to compare current submission (or that corresponding to specified database) against previous relase
# If one submission specified on command line, compare against current (or submission corresponding to specified database)
# If two submissions specified then compare against each other

my $ena_dir = dir($wormbase->submit_repos);
chdir($ena_dir->stringify);

$current_release =~ s/^WS// if $current_release;
$previous_release =~ s/^WS// if $previous_release;

unless ($current_release and $previous_release and $current_iteration and $previous_iteration) {
    my $iterations = get_latest_release_iterations($log);
    
    $current_release = $wormbase->get_wormbase_version unless $current_release;
    $previous_release = $current_release - 1 unless $previous_release;
    $current_iteration = $iterations->{$current_release} unless $current_iteration;
    $previous_iteration = $iterations->{$previous_release} unless $previous_iteration;
}

my $current_tag = 'WS' . $current_release . '.' . $current_iteration;
my $previous_tag = 'WS' . $previous_release . '.' . $previous_iteration;

my $chr_dir = dir("$ena_dir/0");
my $new_genes = [];
my $removed_genes = [];
my ($shared_genes, $current_data, $previous_data);
for my $child ($chr_dir->children(no_hidden => 1)) {
    for my $grandchild ($child->children(no_hidden => 1)) {
	next unless $grandchild->stringify =~ /\.embl$/;
	my $filename = '0/' . $child->basename . '/' . $grandchild->basename;
	
	($new_genes, $removed_genes, $shared_genes,
	 $current_data, $previous_data) = compare_submissions($filename, $current_tag,
							      $previous_tag, $new_genes,
							      $removed_genes, $shared_genes,
							      $log);
	
    }
}

$log->write_to(@$new_genes . " genes added to $current_tag since $previous_tag:\n" .
	       join(', ', @$new_genes) . "\n\n");
$log->write_to(@$removed_genes . " genes previously in $previous_tag " .
	       "no longer found in $current_tag:\n" .
	       join(', ', @$removed_genes) . "\n\n");

my $changes = {};
for my $gene (@$shared_genes) {
    $changes = compare_genes($current_data, $previous_data, $gene, $changes, $log);
}


$log->write_to(scalar (keys %{$changes->{'gene'}}) . 
	       " genes with changes between $current_tag and $previous_tag:\n" .
	       join(', ', keys %{$changes->{'gene'}}) . "\n\n");
for my $ft (keys %$changes) {
    next if $ft eq 'gene';

    $log->write_to("\t$ft:\n\n");
    $log->write_to("\t\t" . scalar (keys %{$changes->{$ft}{'new'}}) . 
		   " genes with new ${ft}s in $current_tag:\n" .
		   "\t\t" . join(', ', keys %{$changes->{$ft}{'new'}}) . "\n\n");

    $log->write_to("\t\t" . scalar (keys %{$changes->{$ft}{'removed'}}) . 
		   " genes with ${ft}s removed since $previous_tag:\n" .
		   "\t\t" . join(', ', keys %{$changes->{$ft}{'removed'}}) . "\n\n");


    $log->write_to("\t\t" . scalar (keys %{$changes->{$ft}{'changed'}}) . 
		   " genes with changed ${ft}s:\n" .
		   "\t\t" . join(', ', keys %{$changes->{$ft}{'changed'}}) . "\n\n");

    my $nr_attr_changes = count_attribute_changes($changes, $ft);
    for my $attr (keys %$nr_attr_changes) {
	if (defined $nr_attr_changes->{$attr}{'new'}{'features'}) {
	    $log->write_to("\t\t\t" . $nr_attr_changes->{$attr}{'new'}{'features'} .
			   " new $ft $attr attributes found in " .
			   $nr_attr_changes->{$attr}{'new'}{'genes'} .
			   " gene(s) in $current_tag\n\n");
	}
	if (defined $nr_attr_changes->{$attr}{'removed'}{'features'}) {
	$log->write_to("\t\t\t" . $nr_attr_changes->{$attr}{'removed'}{'features'} .
		       " $ft $attr attributes removed from " . 
		       $nr_attr_changes->{$attr}{'removed'}{'genes'} .
		       " gene(s) since $previous_tag\n\n");
	}

	if (defined $nr_attr_changes->{$attr}{'changed'}{'features'}) {
	    $log->write_to("\t\t\t" . $nr_attr_changes->{$attr}{'changed'}{'features'} .
			   " $ft $attr attributes changed in " .
			   $nr_attr_changes->{$attr}{'changed'}{'genes'} .
			   " gene(s) between $previous_tag and $current_tag\n\n")
	}
    }
}

$log->mail;

exit(0);


sub compare_genes {
    my ($current_data, $previous_data, $gene, $changes, $log) = @_;

    
    for my $ft ('mRNA', 'CDS', 'ncRNA', 'prim_transcript',
		'rRNA', 'tRNA', 'precursor_RNA') {
	for my $sn (keys %{$current_data->{$gene}{$ft}}) {
	    if (exists $previous_data->{$gene}{$ft}{$sn}) {
		for my $attr (keys %{$current_data->{$gene}{$ft}{$sn}}) {
		    if (exists $previous_data->{$gene}{$ft}{$sn}{$attr}) {
			unless ($current_data->{$gene}{$ft}{$sn}{$attr} eq
			    $previous_data->{$gene}{$ft}{$sn}{$attr}) {
			    $changes->{$ft}{'changed'}{$gene}{$attr}{'changed'}++;
			    $changes->{'gene'}{$gene}++;
			}   
		    }
		    else {
			$changes->{$ft}{'changed'}{$gene}{$attr}{'new'}++;
			$changes->{'gene'}{$gene}++;
		    }
		}

		for my $attr (keys %{$previous_data->{$gene}{$ft}{$sn}}) {
		    unless (exists $current_data->{$gene}{$ft}{$sn}{$attr}) {
			$changes->{$ft}{'changed'}{$gene}{$attr}{'removed'}++;
			$changes->{'gene'}{$gene}++;
		    }
		}
	    }
	    else {
		$changes->{$ft}{'new'}{$gene}++;
		$changes->{'gene'}{$gene}++;
	    }
	}

	for my $sn (keys %{$previous_data->{$gene}{$ft}}) {
	    unless (exists $current_data->{$gene}{$ft}{$sn}) {
		$changes->{$ft}{'removed'}{$gene}++;
		$changes->{'gene'}{$gene}++;
	    }
	}
    }

    return $changes;
}


sub count_attribute_changes {
    my ($changes, $ft) = @_;

    my (%gene_counts, %ft_counts);
    for my $gene (keys %{$changes->{$ft}{'changed'}}) {
	for my $attr (keys %{$changes->{$ft}{'changed'}{$gene}}) {
	    for my $status (keys %{$changes->{$ft}{'changed'}{$gene}{$attr}}) {
		$gene_counts{$attr}{$status}{$gene} = 1;
		$ft_counts{$attr}{$status} += $changes->{$ft}{'changed'}{$gene}{$attr}{$status};
	    }
	}
    }

    my %summary;
    for my $status ('new', 'removed', 'changed') {
	for my $attr (keys %gene_counts) {
	    $summary{$attr}{$status}{'genes'} = scalar keys %{$gene_counts{$attr}{$status}};
	    $summary{$attr}{$status}{'features'} = $ft_counts{$attr}{$status};
	}
    }

    return \%summary;
}


sub compare_submissions {
    my ($file, $current_tag, $previous_tag, $new_genes, $removed_genes, $shared_genes, $log) = @_;

    my $current_data = parse_file($file, $current_tag, $log);
    my $previous_data = parse_file($file, $previous_tag, $log);

    for my $current_gene (keys %$current_data) {
	if (exists $previous_data->{$current_gene}) {
	    push @$shared_genes, $current_gene;
	}
	else {
	    push @$new_genes, $current_gene;
	}
    }

    for my $previous_gene (keys %$previous_data) {
	push @$removed_genes, $previous_gene unless exists $current_data->{$previous_gene};
    }

    return ($new_genes, $removed_genes, $shared_genes, $current_data, $previous_data);
}


sub get_latest_release_iterations {
    my $log = shift;

    my %iterations;

    open (LOG, "git tag|") or $log->log_and_die($!);
    while (<LOG>) {
	chomp;
	next unless $_ =~ /^WS(\d{3})\.(\d+)$/;
	$iterations{$1} = $2 unless exists $iterations{$1} and $iterations{$1} > $2;
    }
    close (LOG);

    return \%iterations;
}


sub parse_file {
    my ($file, $tag, $log) = @_;

    my ($sn, $gene, $entry, $feature, $attr);
    my %data;
    open (SUB, "git show $tag:$file |");
    while (<SUB>) {
	chomp;
	next unless $_ =~ /^FT/;
	if ($_ =~ /^FT\s{3}(\S+)\s+(\S+)$/) {
	    if ($gene) {
		$data{$gene}{$feature}{$sn} = $entry;
	    }
	    $entry = {};
	    undef $gene;
	    undef $sn;
	    $attr = 'position';
	    $entry->{'position'} = $2;
	    $feature = $1;
	}
	elsif ($_ =~ /^FT\s{19}\/(\w+)=(.+)$/) {
	    $attr = $1;
	    if ($attr eq 'gene') {
		$gene = $2;
		$gene =~ s/"//g;
	    }
	    if ($attr eq 'standard_name') {
		$sn = $2;
		$sn =~ s/"//g;
	    }
	    else {
		$entry->{$attr} = $2;
	    }
	}
	else {
	    next unless $_ =~ /^FT\s{19}(.+)$/;
	    $entry->{$attr} .= $1;
	}
    }
    $data{$gene}{$feature}{$sn} = $entry;
    close (SUB);

    return \%data;
}
