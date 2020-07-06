#!/usr/bin/env perl

use strict;
use warnings;
use List::MoreUtils qw(uniq);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Data::Dumper;

my ($host, $user, $db, $pass, $port, $ipcress_outfile, $ipcress_fh, $phenotypes_file, $phenotypes_fh, $out, $out_fh, $intended_targets_file, $intended_targets_fh);
my (%phenotypes, %intended_targets, %ipcress_targets);


GetOptions(
	'host=s'			=> \$host,
	'dbname=s'			=> \$db,
	'pass=s'			=> \$pass,
	'port=i'			=> \$port,
	'user=s'			=> \$user,
	'ipcress_outfile=s'		=> \$ipcress_outfile,
	'phenotypes_file=s'             => \$phenotypes_file,
	'out=s'				=> \$out,
	'intended_targets_file=s'	=> \$intended_targets_file
	
) || die ("check command line params\n");


my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $host,
    -user   => $user,
    -dbname => $db,
    -pass   => $pass,
    -port   => $port,
      );

open($ipcress_fh, '<', $ipcress_outfile) or die "Could not open $ipcress_outfile\n";
open($phenotypes_fh, '<', $phenotypes_file) or die "Could not open $phenotypes_file\n";
if (defined $intended_targets_file){
	open($intended_targets_fh, '<', $intended_targets_file) or die "Could not open $intended_targets_file\n";
}
open($out_fh, '>', $out) or die "Could not create $out for writing\n";

# parse phenotypes file

while (<$phenotypes_fh>){
	chomp;
	my @temp = split(/\t/,$_);
	if (scalar @temp != 2) { die "phenotype file should have 2 fields: primer_id\tphenotype"; }
	my $primer_id = $temp[0];
	my $phenotype = $temp[1]; 	
	push @{ $phenotypes{$primer_id} }, $phenotype;
}

# parse intended targets

if (defined $intended_targets_file ){
	while(<$intended_targets_fh>){
		chomp;
        	my @temp = split(/\t/,$_);
        	if (scalar @temp != 2) { die "targets file should have 2 fields: primer_id\tintended target"; }
		my $primer_id   = $temp[0];
		my $target_gene = $temp[1];
		push @{ $intended_targets{$primer_id} }, $target_gene; 
	}
}


# prepare sql - given transcript ID, retrieve gene ID

my $dbh = $dba->dbc->db_handle;

my $sql = "SELECT gene.stable_id FROM gene 
	LEFT JOIN transcript ON gene.gene_id = transcript.gene_id 
	WHERE transcript.stable_id = ?";

# parse ipcress outfile : ipcress: Smp_000070.1:filter(unmasked) V5.PRIMER.ID:Smp_000070 704 A 296 0 B 980 0 forward

while(<$ipcress_fh>){
	
	next unless ($_ =~ /^ipcress:/); 
	chomp;
	my %data;
	my @temp	 		 = split(/\s/,$_);
	my $transcript			 = $temp[1];
	my $primer_id 			 = $temp[2];
	$data{'product_length'}		 = $temp[3];
	my $description 		 = $temp[10];
	my $type;
	$transcript =~ s/:filter\(unmasked\)//;

	my @row = $dbh->selectrow_array($sql, undef, $transcript);
	unless(@row){ die "transcript $transcript not found in $db"; }
	my $gene = $row[0];

	if ($description eq 'forward' || $description eq 'revcomp'){
		$type = 'good_product';
	}

	# we don't want bad products
	
	elsif ($description eq 'single_A' || $description eq 'single_B'){
		next;
	}

	else{ die "description $description not recognised"; }
	
	$data{'transcript'} = $transcript;
	$data{'gene'}	    = $gene;
	$data{'type'}	    = $type;

	# check if the target matches (at least one of) the targets intended by the author
	
	if (defined $intended_targets_file){
		foreach my $author_gene (@{ $intended_targets{$primer_id}}){
			if ($author_gene eq $gene){
				$data{'match'} = 'match';
				last;
			}
			$data{'match'} = 'no_match'
		}
	}

	push @{ $ipcress_targets{$primer_id}}, \%data;
}

# for each primer set, check if it targets multiple genes

foreach my $primer_id (keys %ipcress_targets){
	my @genes;
	foreach my $ipcress_target ( @ {$ipcress_targets{$primer_id}} ){
		
		push @genes, $ipcress_target->{'gene'};
	}
	if (scalar uniq(@genes) == 1){
		foreach my $ipcress_target ( @ {$ipcress_targets{$primer_id}} ){
			$ipcress_target->{'specificity'} = 'unique';
		}
	}
	else{
                foreach my $ipcress_target ( @ {$ipcress_targets{$primer_id}} ){
                        $ipcress_target->{'specificity'} = 'multi';
                }
	}
}


# print some stats on the predictions

my $no_target      = 0;
my $unique_target  = 0;
my $multi_target   = 0;
my $match 	   = 0;

foreach my $primer_id (keys %phenotypes){
	
	if (! exists $ipcress_targets{$primer_id}  ){
		$no_target++;
		next;	
	}
	
	if ($ipcress_targets{$primer_id}[0]->{'specificity'} eq 'unique'){
		$unique_target++;
	}
	if (($ipcress_targets{$primer_id}[0]->{'specificity'} eq 'unique') and
	    (defined $intended_targets_file) and
	    ($ipcress_targets{$primer_id}[0]->{'match'} eq 'match')){
		$match++;
	}

	elsif ($ipcress_targets{$primer_id}[0]->{'specificity'} eq 'multi'){
		$multi_target++;
	}

	foreach my $ipcress_target ( @ {$ipcress_targets{$primer_id}} ){
		foreach my $phenotype ( @{ $phenotypes{$primer_id}} ){
			my $print = join "\t", ($ipcress_target->{'gene'},
				     $phenotype,
		 		     $ipcress_target->{'transcript'},
				     $ipcress_target->{'product_length'},
				     $ipcress_target->{'type'},
				     $ipcress_target->{'specificity'},
				     $primer_id
			);
			if (defined $intended_targets_file) {
				$print = join "\t", ($print,$ipcress_target->{'match'});
			}
		print {$out_fh} "$print\n";
		}
	}
}

print "ipcress predicted:
Primers with unique targets    : $unique_target
Primers with no target         : $no_target
Primer with multiple targets   : $multi_target
";


if (defined $intended_targets_file){
	print "Of the unique targets, $match matched the author's intended target\n";
}

