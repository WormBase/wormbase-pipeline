#!/usr/bin/perl

use strict;
use warnings;

use lib '../lib';

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use DBI;
use Domain2Interpro;
use Getopt::Long;

my $debug = 1;

# Vega database
my $vega_dbname = '';
my $vega_dbhost = 'ia64b';
my $vega_dbport =  3306;
my $vega_dbuser = 'wormro'; 
my $vega_dbpass = '';


GetOptions(
           'dbhost=s'     => \$vega_dbhost,
           'dbname=s'     => \$vega_dbname,
           'dbuser=s'     => \$vega_dbuser,
           'dbpass=s'     => \$vega_dbpass,
           'dbport=s'     => \$vega_dbport,
          )or die ("Couldn't get options");

my $vega_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
	-dbname  => $vega_dbname,
	-host    => $vega_dbhost,
	-port    => $vega_dbport,
	-user    => $vega_dbuser,
	-pass    => $vega_dbpass,
#	-species => 'vega',
);

my $mapper = Domain2Interpro->new();

print "##gff-version 3\n";

# Get all Vega toplevel slices
my $sa = $vega_db->get_SliceAdaptor();
my $slices = $sa->fetch_all('toplevel');
foreach my $slice (@$slices) {
	my $slice_name = $slice->seq_region_name();
	print STDERR "Doing $slice_name:\n" if $debug;
	
	# Get all the genes on this slice
	my $genes = $slice->get_all_Genes();
	foreach my $gene (@$genes) {
		my $gene_stable_id = 'gene.'.$gene->stable_id();
		print STDERR " Dumping $gene_stable_id\n" if $debug;
		my %gene_to_dump = (
			stable_id => $gene_stable_id,
			name      => $slice_name,
			start     => $gene->start(),
			end       => $gene->end(),
			strand    => $gene->strand(),
#			display   => $gene->external_name(),
			note      => $gene->status() . " " . $gene->biotype(),
			public_name => $gene->stable_id(),
		);
		foreach my $transcript (@{$gene->get_all_Transcripts()}) {
			my $transcript_stable_id = 'transcript.'.$transcript->stable_id();
			print STDERR "  Dumping $transcript_stable_id\n" if $debug;
			push @{$gene_to_dump{'transcript'}}, {
				stable_id => $transcript_stable_id,
				name      => $slice_name,
				start     => $transcript->start(),
				end       => $transcript->end(),
				strand    => $transcript->strand(),
#				display   => $transcript->external_name(),
				note      => $transcript->biotype(),
				info      => get_info($transcript),
				public_name => $transcript->stable_id(),
			};
			foreach my $exon (@{$transcript->get_all_Exons()}) {
				my $exon_stable_id = 'exon.'.$exon->stable_id();
				print STDERR "   Dumping $exon_stable_id\n" if $debug;
				push @{${${$gene_to_dump{'transcript'}}[-1]}{'exon'}}, {
					stable_id => $exon_stable_id,
					name      => $slice_name,
					start     => $exon->start(),
					end       => $exon->end(),
					strand    => $exon->strand(),
				}
			}
			foreach my $cds (@{$transcript->get_all_translateable_Exons()}) {
#				my $cds_stable_id = 'cds.'.$transcript->translation->stable_id();
				my $cds_stable_id = 'coding_exon.'.$cds->stable_id();

				print STDERR "   Dumping CDS $cds_stable_id\n" if $debug;
				push @{${${$gene_to_dump{'transcript'}}[-1]}{'cds'}}, {
					stable_id => $cds_stable_id,
					name      => $slice_name,
					start     => $cds->start(),
					end       => $cds->end(),
					strand    => $cds->strand(),
					phase     => $cds->phase(),
# add blastp hits and domains as notes?
				}
			}
		}
		print dump_gene(\%gene_to_dump);
	}
	my $features = $slice->get_all_ProteinAlignFeatures();
	print STDERR "dumping ${\scalar @$features} Features\n";
	foreach my $feature(@$features){
		my $stripped_feature = {
			hit_id      => $feature->hseqname,
			target_id   => $feature->slice->seq_region_name,
			target_start=> $feature->hstart,
			target_stop => $feature->hend,
			strand      => ($feature->strand > 0?'+':'-'),
			hit_start   => $feature->start,
			hit_stop    => $feature->end,
			score       => $feature->score,
			dbid        => $feature->dbID,
			logic_name  => $feature->analysis->logic_name,
			cigar       => $feature->cigar_string,
		};
		print dump_feature($stripped_feature);
	}
	print '#'x80;
	print "\n";

}


sub dump_feature {
	my $i=shift;
	my %feature=%{$i};
	my $gff_line=
	 "$feature{target_id}\tprotein_alignment\t$feature{logic_name}\t$feature{target_start}\t$feature{target_stop}\t$feature{score}\t$feature{strand}\t.\t"
	."ID=$feature{logic_name}.$feature{dbid};Name=$feature{hit_id};Target=$feature{hit_id} $feature{hit_start} $feature{hit_stop} Gap=$feature{cigar}\n";
	return $gff_line;
}

sub get_info {
	my $transcript= shift;
	my $info='';
	# get all protein_features on the transcript
	my $features=$transcript->translation->get_all_ProteinFeatures();
	# get logic_name and hit_id
	
	my %plain_features;
        my $rest_features;
	map {
		if ($mapper->get_method2database($_->analysis->logic_name())) {
			push @{$plain_features{$_->analysis->logic_name()}},[$_->id(),$_->start(),$_->end(),$_->hstart(),$_->hend(),$_->score(),$_->p_value()]
		}else {
			push @$rest_features,$_
		}
	} @$features;
	my @interpros=$mapper->get_mapping(\%plain_features);
	map {$info.=sprintf( "position:%d-%d method:%s accession:%s description:%s %%0A", $_->[1], $_->[2], 'InterPro', $_->[0] , $_->[7]) if $_->[1]} @interpros;

	while ( my $pfeature = shift @$rest_features ) {
         my $logic_name = $pfeature->analysis()->logic_name();
         $info.=sprintf( "position:%d-%d %s method:%s accession:%s %%0A", $pfeature->start(), $pfeature->end(), $pfeature->p_value(),$logic_name, $pfeature->id());
	}
	return $info;
}

sub dump_gene {
	my ($gene) = @_;
	
	my $output = '';
	
	# Dump gene
	$output .= "# Gene " . $gene->{'stable_id'} . "\n";
	$output .= gff_line(
		$gene->{'name'}, 'gene', $gene->{'start'}, $gene->{'end'},
		$gene->{'strand'}, $gene->{'stable_id'},undef, $gene->{'display'}, $gene->{'note'},undef,$gene->{'public_name'}
	);
	
	# Dump transcripts
	my $parent = $gene->{'stable_id'};
	my %exon_parent;
	foreach my $transcript (@{$gene->{'transcript'}}) {
		$output .= gff_line(
			$transcript->{'name'}, 'mRNA', $transcript->{'start'}, $transcript->{'end'},
			$transcript->{'strand'}, $transcript->{'stable_id'}, undef ,$transcript->{'display'} || undef, undef,$transcript->{info} 
			|| undef,$transcript->{'public_name'},$parent);
		
		# Store the parent of this transcript's exons
		$parent = $transcript->{'stable_id'};
		foreach my $exon (@{$transcript->{'exon'}}) {
			${$exon_parent{$exon->{'stable_id'}}}{$parent} = 1;
		}
	}
	
	# Dump exons
	foreach my $transcript (@{$gene->{'transcript'}}) {
		foreach my $exon (@{$transcript->{'exon'}}) {
			next if !$exon_parent{$exon->{'stable_id'}}; # If there are no parents then we've already dumped this exon
			my @parents = keys %{$exon_parent{$exon->{'stable_id'}}};
			delete $exon_parent{$exon->{'stable_id'}};
			$output .= gff_line(
				$exon->{'name'}, 'exon', $exon->{'start'}, $exon->{'end'},
				$exon->{'strand'}, $exon->{'stable_id'}, undef,undef, undef,undef,undef ,@parents);
		}
	}
	
	# Dump coding_exons
	foreach my $transcript (@{$gene->{'transcript'}}) {
		my $parent = $transcript->{'stable_id'};
		foreach my $cds (@{$transcript->{'cds'}}) {
			$output .= gff_line(
				$cds->{'name'}, 'coding_exon', $cds->{'start'}, $cds->{'end'},
				$cds->{'strand'}, $cds->{'stable_id'}, $cds->{'phase'},undef, undef,undef ,undef, $parent);
		}
	}
	
	$output .= "###\n";
	
	return $output;
}

sub gff_line {
	my ($seqid, $type, $start, $end, $strand, $stable_id, $phase,$name, $note, $info,$public,@parents) = @_;
	
	my $output = '';

	$phase='.' unless defined $phase;
        $strand= $strand>0? '+' : '-';

	$output .= "$seqid\tWormBase\t$type\t$start\t$end\t.\t$strand\t$phase\t";
	my @tags;
	push @tags, "ID=$stable_id" if defined $stable_id;
	if (@parents) {
		my $parent = join(',', @parents);
		push @tags, "Parent=$parent";
	}
	push @tags, "Name=$name" if defined $name;
	push @tags, "Note=$note" if defined $note;
	push @tags, "Info=$info" if defined $info;
	push @tags, "Public_name=$public" if defined $public;
	$output .= join(';', @tags);
	$output .= "\n";
	
	return $output;
}

