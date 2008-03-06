#!/usr/bin/perl
# based on the ZFish GFF3 VEGA dumper from Ian Sealy
# the Blast Filter is lifted from the WormBlast processing pipeline
#
# To lower the memory footprint and converts the rather large
# Ensembl objects into lighter hashes for later use


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

# Get all toplevel slices
my $sa = $vega_db->get_SliceAdaptor();
my $slices = $sa->fetch_all('toplevel');
while( my $slice = shift @$slices) {
	my $slice_name = $slice->seq_region_name();
	my $slice_size = $slice->length;

	print STDERR "Doing $slice_name:\n" if $debug;

	print "##sequence-region $slice_name 1 $slice_size\n";
	
	# Get all the genes on this slice
	my $genes = $slice->get_all_Genes();
	while( my $gene=shift @$genes) {
		my $gene_stable_id = 'gene:'.$gene->stable_id();
		print STDERR " Dumping $gene_stable_id\n" if $debug;

		my %gene_to_dump = (
			stable_id => $gene_stable_id,
			name      => $slice_name,
			start     => $gene->start(),
			end       => $gene->end(),
			strand    => $gene->strand(),
			note      => ($gene->status()||'PREDICTED' ). " " . $gene->biotype(),
			public_name => $gene->stable_id(),
			display     => $geen->stable_id(), # wrong but fixes db's without xref_mapping
		);

		# get all transcripts of the gene
		my $all_transcripts = $gene->get_all_Transcripts();
		while( my $transcript = shift @{$all_transcripts}) {
			my $translation = $transcript->translation;

			my $transcript_stable_id = 'transcript:'.$transcript->stable_id();
			my $translation_id = 'cds:'.$transcript->translation->stable_id();

			print STDERR "  Dumping $transcript_stable_id\n" if $debug;
			print STDERR "  Dumping $translation_id\n" if $debug;

			push @{$gene_to_dump{'transcript'}}, {
				stable_id => $transcript_stable_id,
				translation_stable_id => $translation_id,
				name      => $slice_name,
				start     => $transcript->start(),
				cds_start => $translation->genomic_start(),
				cds_end   => $translation->genomic_end(),
				end       => $transcript->end(),
				strand    => $transcript->strand(),
				note      => $transcript->biotype(),
				info      => get_info($transcript),
				public_name => $transcript->stable_id(),
				display     => $transcript->stable_id(),
			};

			# get all exons and translatable_exons of the transcript
			my $all_exons=$transcript->get_all_Exons();
			my $all_t_exons = $transcript->get_all_translateable_Exons();

			while( my $exon = shift @{$all_exons}) {
				my $exon_stable_id = 'exon:'.$exon->stable_id();
				print STDERR "   Dumping $exon_stable_id\n" if $debug;
				push @{${${$gene_to_dump{'transcript'}}[-1]}{'exon'}}, {
					stable_id => $exon_stable_id,
					name      => $slice_name,
					start     => $exon->start(),
					end       => $exon->end(),
					strand    => $exon->strand(),
				}
			}
			
			$all_exons=undef;

			while (my $cds = shift @{$all_t_exons}) {
				push @{${${$gene_to_dump{'transcript'}}[-1]}{'cds'}}, {
					stable_id => $translation_id,
					name      => $slice_name,
					start     => $cds->start(),
					end       => $cds->end(),
					strand    => $cds->strand(),
					phase     => (3-$cds->phase())%3, # phase/frame conversion to a sane system
				}
			}
		}
		print dump_gene(\%gene_to_dump);
	}

	# get all protein align features on the slice
	my $features = $slice->get_all_ProteinAlignFeatures();
        my %blastx_features;

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
			p_value     => $feature->p_value,
			dbid        => $feature->dbID,
			logic_name  => $feature->analysis->logic_name,
			cigar       => $feature->cigar_string,
			feature_type=> 'protein_match',
		};
		push @{$blastx_features{$stripped_feature->{logic_name}}}, $stripped_feature;
	}
        $features=undef;

	while (my($k,$v)=each %blastx_features){
		  my @filtered_features=filter_features($v,$slice_size);
		  map {print dump_feature($_)} @filtered_features;
	  }


	# get all repeat features on the slice
        my $repeats = $slice->get_all_RepeatFeatures;
        foreach my $feature (@$repeats){
        	my $stripped_feature = {
			target_id   => $feature->slice->seq_region_name,
			strand      => ($feature->strand > 0?'+':'-'),
			hit_start   => $feature->start,
			hit_stop    => $feature->end,
			score       => ($feature->score||'.'),
			dbid        => $feature->dbID,
			logic_name  => $feature->analysis->logic_name,
			feature_type=> 'repeat_region',
		};
		print dump_feature($stripped_feature);
	}

	print '#'x80;
	print "\n";

}

# CIGAR to old GFF3 CIGAR format converter
sub cigar_to_almost_cigar {
	my $i=shift;
	$i=~s/(\d*)(\w)/$2$1 /g;
	return $i;
}

# print the feature using some funky template
sub dump_feature {
	my $i=shift;
	my %feature=%{$i};
	my $gff_line=
	 "$feature{target_id}\t$feature{logic_name}\t$feature{feature_type}\t$feature{hit_start}\t$feature{hit_stop}\t".
	 "$feature{score}\t$feature{strand}\t.".
	 "\tID=$feature{logic_name}.$feature{dbid}".
	 ($feature{cigar}?";Name=$feature{hit_id};Target=$feature{hit_id} $feature{target_start} $feature{target_stop};Gap=$feature{cigar}\n":"\n");
	return $gff_line;
}


# build the info tag including protein features and interpro
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
			push @{$plain_features{$_->analysis->logic_name()}},[$_->id(),$_->start(),
			$_->end(),$_->hstart(),$_->hend(),$_->score(),$_->p_value()]
		}else {
			push @$rest_features,$_
		}
	} @$features;
	my @interpros=$mapper->get_mapping(\%plain_features);
	map {$info.=sprintf( "position:%d-%d method:%s accession:%s description:%s %%0A", $_->[1], $_->[2],
			'InterPro', $_->[0] , $_->[7]) if $_->[1]} @interpros;

	while ( my $pfeature = shift @$rest_features ) {
         my $logic_name = $pfeature->analysis()->logic_name();
         $info.=sprintf( "position:%d-%d %s method:%s accession:%s %%0A", $pfeature->start(), $pfeature->end(), 
		 $pfeature->p_value(),$logic_name, $pfeature->id());
	}
	return $info;
}

# print the gene including transcripts and exons
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
			$transcript->{'strand'}, $transcript->{'stable_id'}, undef ,$transcript->{'display'} || undef, undef, 
			($transcript->{info}||undef),$transcript->{'public_name'},$parent);

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
				$cds->{'name'}, 'CDS', $cds->{'start'}, $cds->{'end'},
				$cds->{'strand'}, $cds->{'stable_id'}, $cds->{'phase'},undef, 
				undef,undef ,undef, $transcript->{'stable_id'});
		}
	}
	
	$output .= "###\n";
	
	return $output;
}

# a template for a GFF line
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
	push @tags, "info=$info" if defined $info;
	push @tags, "public_name=$public" if defined $public;
	$output .= join(';', @tags);
	$output .= "\n";
	
	return $output;
}

# remove < 80% of evalue features from 100bp windows
sub filter_features {
	my ($features,$length)=@_;
        my %f_features;
	
	# 50bp bin size
	my $size=100;

	# should I bin them instead?
	my @bins;

	# put the feature in all bins between start and stop
	foreach my $f(@$features){
		my ($_start,$_end)=sort {$a <=> $b } ($f->{hit_start},$f->{hit_stop});
		for (my $n=int($_start/$size);$n <=int($_end/$size);$n++){
			push @{$bins[$n]},$f;
		}
	}

	# get the best 5 within 25% of the best hsp and add them to a hash
	foreach my $bin(@bins) {
		next unless $bin; # skip empty bins
		my $best=0;
		my $max_hsp=0;
		my @sorted_bin = sort {$a->{p_value} <=> $b->{p_value}} @$bin;
		$best=&p_value($sorted_bin[0]->{p_value});
		map { $f_features{$_->{dbid}}=$_ if (&p_value($_->{p_value}) > $best*0.80 && $max_hsp++ <5)} @$bin; # <= cutoff place, lets try 20% instead of 25%
	}
	
	# flatten hash into an array
	my @_filtered=values %f_features;
	return @_filtered;
}

# p_value shorthand
sub p_value {
        my ($p)=@_;
        my $log = ($p > 0 ? -(log($p))/log(10) : 999.9); # if e-value is zero set it to 999.9
        return $log;
}



